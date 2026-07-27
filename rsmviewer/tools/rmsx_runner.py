#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rmsx_runner.py

Drives the RNAMotifScanX pipeline for one PDB ID.

Pipeline (Zhong & Zhang, RNA 21:333-346, 2015):
  1. Annotate the target PDB with MC-Annotate and RNAVIEW, merge the two
     interaction sets (union; MC-Annotate wins on conflict), and write the
     per-chain .rmsx.in target file
  2. For each configured motif family, run RNAMotifScanX with the family's
     consensus query against the annotated target
  3. Save output as  {output_dir}/{motif_family}_consensus/result_0_100_withbs.log

RSMViewer uses these files when the user runs:
    rmv_db 7 /path/to/rmsx_pipeline_config.json
    rmv_fetch 1S72
    rmv_load_motif

Usage (standalone):
    python rmsx_runner.py --config rmsx_pipeline_config.json --pdb 1S72
    python rmsx_runner.py --config rmsx_pipeline_config.json --pdb 1S72 --chains 0

Platform note:
    The RNAMotifScanX binary (release_v0.0.5_x86-64_rhel) is Linux x86-64 only.
    On macOS, place pre-generated .log files in the output directory manually;
    RSMViewer will load them without needing the executable.
"""

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


# Motif family name → typical consensus file name fragment
# These names mirror the folder structure expected by RSMViewer's RMSX converter.
DEFAULT_FAMILY_FOLDER_MAP = {
    'k-turn':        'k-turn_consensus',
    'c-loop':        'c-loop_consensus',
    'sarcin-ricin':  'sarcin-ricin_consensus',
    'reverse-kturn': 'reverse-kturn_consensus',
    'e-loop':        'e-loop_consensus',
}


# ─────────────────────────────────────────────────────────────────────────────
# CIF download helper (same SSL-fallback as fr3d_loop_extractor)
# ─────────────────────────────────────────────────────────────────────────────

def _download_cif(pdb_id: str, dest_dir: str) -> str:
    import ssl, urllib.request
    pdb_lower = pdb_id.lower()
    url = f'https://files.rcsb.org/download/{pdb_lower}.cif.gz'
    dest = os.path.join(dest_dir, f'{pdb_lower}.cif.gz')
    print(f"[rmsx_runner] Downloading CIF: {url}")
    last_exc = None
    for ctx in [ssl.create_default_context(), ssl._create_unverified_context()]:
        try:
            with urllib.request.urlopen(url, context=ctx, timeout=60) as r:
                data = r.read()
            with open(dest, 'wb') as fh:
                fh.write(data)
            return dest
        except Exception as exc:
            last_exc = exc
    print(f"[rmsx_runner] Download failed: {last_exc}")
    return ''


def _download_pdb(pdb_id: str, dest_dir: str) -> str:
    import ssl
    import urllib.request
    pdb_upper = pdb_id.upper()
    url = f'https://files.rcsb.org/download/{pdb_upper}.pdb'
    dest = os.path.join(dest_dir, f'{pdb_upper}.pdb')
    print(f"[rmsx_runner] Downloading PDB: {url}")
    last_exc = None
    for ctx in [ssl.create_default_context(), ssl._create_unverified_context()]:
        try:
            with urllib.request.urlopen(url, context=ctx, timeout=60) as r:
                data = r.read()
            with open(dest, 'wb') as fh:
                fh.write(data)
            return dest
        except Exception as exc:
            last_exc = exc
    print(f"[rmsx_runner] PDB download failed: {last_exc}")
    return ''


def _extract_prebuilt_targets_from_archive(prebuild_archive: str, pdb_id: str,
                                           output_dir: str, chains: list) -> dict:
    """Extract prebuilt .rmsx.in/.rmsx.nch files for pdb_id from a tgz archive.

    Returns {chain_id: path}. Empty dict when no matching prebuilt target exists.
    """
    archive = os.path.expanduser(str(prebuild_archive or '').strip())
    if not archive or not os.path.isfile(archive):
        return {}

    pdb_upper = str(pdb_id).strip().upper()
    prefix = f"{pdb_upper}_"
    wanted = {str(c).strip() for c in (chains or []) if str(c).strip()}
    accept_any_chain = (not wanted) or ('0' in wanted)

    def _ext_priority(name: str) -> int:
        upper = name.upper()
        if upper.endswith('.RMSX.IN'):
            return 0
        if upper.endswith('.RMSX.NCH'):
            return 1
        return 2

    produced = {}
    stage_dir = os.path.join(output_dir, '_prebuilt_targets')
    os.makedirs(stage_dir, exist_ok=True)

    try:
        with tarfile.open(archive, 'r:*') as tf:
            candidates = []
            for member in tf.getmembers():
                if not member.isfile():
                    continue
                base = os.path.basename(member.name)
                upper = base.upper()
                if not upper.startswith(prefix):
                    continue
                if not (upper.endswith('.RMSX.IN') or upper.endswith('.RMSX.NCH')):
                    continue

                # Example: 4V9F_B.rmsx.in -> chain token "B"
                tail = base[len(prefix):]
                chain_token = tail.split('.', 1)[0]
                if (not accept_any_chain) and (chain_token not in wanted):
                    continue
                candidates.append((member, chain_token))

            candidates.sort(key=lambda item: (_ext_priority(item[0].name), item[0].name))
            for member, chain_token in candidates:
                if chain_token in produced:
                    continue
                dst = os.path.join(stage_dir, os.path.basename(member.name))
                if not os.path.isfile(dst):
                    src = tf.extractfile(member)
                    if src is None:
                        continue
                    with src, open(dst, 'wb') as out_fh:
                        out_fh.write(src.read())
                produced[chain_token] = dst

    except Exception as exc:
        print(f"[rmsx_runner] WARNING: could not read prebuilt target archive '{archive}': {exc}")
        return {}

    if produced:
        print(f"[rmsx_runner] Using prebuilt RMSX target files from archive: {archive}")
        for chain_token, path in produced.items():
            print(f"[rmsx_runner]   prebuilt chain {chain_token}: {path}")
    return produced


# ─────────────────────────────────────────────────────────────────────────────
# Annotation
# ─────────────────────────────────────────────────────────────────────────────

def run_mc_annotate(mc_annotate_exe: str, structure_file: str, output_dir: str,
                    pdb_id: str, force_fresh: bool = False) -> str:
    """Run MC-Annotate on structure_file. Returns path to annotation output, '' on failure."""
    out_file = os.path.join(output_dir, f'{pdb_id}_mc_annotate.out')
    if force_fresh and os.path.isfile(out_file):
        try:
            os.remove(out_file)
        except Exception as exc:
            print(f"[rmsx_runner] WARNING: could not remove old annotation file: {out_file} ({exc})")

    if os.path.isfile(out_file):
        print(f"[rmsx_runner] Using existing annotation: {out_file}")
        return out_file

    if not mc_annotate_exe or not os.path.isfile(mc_annotate_exe):
        print(f"[rmsx_runner] WARNING: MC-Annotate executable not found: {mc_annotate_exe}")
        print("[rmsx_runner] Annotation step skipped — RMSX will fail if annotation is required.")
        return ''

    cmd = [mc_annotate_exe, structure_file]
    print(f"[rmsx_runner] Running MC-Annotate: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"[rmsx_runner] MC-Annotate failed (exit {result.returncode})")
            return ''
        if not (result.stdout or '').strip():
            print("[rmsx_runner] MC-Annotate produced empty output.")
            print("[rmsx_runner]   Input may be unsupported (MC-Annotate typically expects PDB-formatted structures).")
            return ''
        with open(out_file, 'w', encoding='utf-8') as fh:
            fh.write(result.stdout)
        print(f"[rmsx_runner] Annotation saved: {out_file}")
        return out_file
    except Exception as exc:
        print(f"[rmsx_runner] MC-Annotate error: {exc}")
        return ''


# ─────────────────────────────────────────────────────────────────────────────
# RNAVIEW annotation (union with MC-Annotate, per Zhong & Zhang 2015 preprocessing)
#
# The reference preprocessing (RNAMotifScanX_src/StructureAnnotation/PrepareInput.py)
# annotates the target with BOTH MC-Annotate AND RNAVIEW, then merges the two sets
# (union; MC-Annotate wins on conflict) before writing the .rmsx.in target file.
# The functions below reproduce that behaviour without altering the RNAMotifScanX
# executable or its algorithms.
# ─────────────────────────────────────────────────────────────────────────────

def _locate_rnaview(config: dict):
    """Locate the RNAVIEW executable and its base directory (BASEPARS resources).

    Executable resolution order:
      1. config['rnaview_executable']
      2. 'rnaview' found on PATH
    Base-directory resolution order (needed for the RNAVIEW environment variable
    that points RNAVIEW at its BASEPARS resource files):
      1. config['rnaview_dir']
      2. RNAVIEW environment variable
      3. two levels up from the executable (…/RNAVIEW/bin/rnaview → …/RNAVIEW)
    Returns (exe_path, base_dir); either element may be '' when not found.
    """
    exe = os.path.expanduser(str(config.get('rnaview_executable', '') or '').strip())
    if exe and not os.path.isfile(exe):
        exe = ''
    if not exe:
        exe = shutil.which('rnaview') or ''

    base_dir = os.path.expanduser(str(config.get('rnaview_dir', '') or '').strip())
    if not base_dir:
        base_dir = os.environ.get('RNAVIEW', '') or ''
    if not base_dir and exe:
        # …/RNAVIEW/bin/rnaview → …/RNAVIEW
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(exe)))
    return exe, base_dir


def _stage_short_rnaview_dir(rnaview_dir: str):
    """Ensure RNAVIEW's BASEPARS path fits its fixed 80-char internal buffer.

    RNAVIEW composes ``<RNAVIEW>//BASEPARS/Atomic_?.pdb`` into ``char spdb[80]``.
    When the resolved directory is deep enough that this path would overflow the
    buffer, expose the directory through a short path (a symlink where supported,
    otherwise a copy of BASEPARS) so the unmodified binary runs safely.

    Returns ``(effective_dir, cleanup_dir_or_None)``. The caller must remove
    ``cleanup_dir`` after RNAVIEW finishes.
    """
    if not rnaview_dir:
        return rnaview_dir, None
    # Longest path RNAVIEW appends after the directory: "//BASEPARS/Atomic_I.pdb".
    projected = len(rnaview_dir) + len('//BASEPARS/Atomic_I.pdb')
    if projected < 78:
        return rnaview_dir, None

    tmp_root = '/tmp' if (os.name != 'nt' and os.path.isdir('/tmp')) else tempfile.gettempdir()
    try:
        staging = tempfile.mkdtemp(prefix='rv', dir=tmp_root)
    except Exception as exc:
        print(f"[rmsx_runner] WARNING: could not create short RNAVIEW staging dir ({exc}); "
              f"using original path (RNAVIEW may abort if it is too long).")
        return rnaview_dir, None

    link = os.path.join(staging, 'r')
    try:
        os.symlink(rnaview_dir, link, target_is_directory=True)
        return link, staging
    except (OSError, NotImplementedError, AttributeError):
        # Symlinks unavailable (e.g. unprivileged Windows): copy BASEPARS instead.
        real_basepars = os.path.join(rnaview_dir, 'BASEPARS')
        try:
            os.makedirs(link, exist_ok=True)
            shutil.copytree(real_basepars, os.path.join(link, 'BASEPARS'))
            return link, staging
        except Exception as exc:
            print(f"[rmsx_runner] WARNING: could not stage a short BASEPARS path ({exc}); "
                  f"using original path (RNAVIEW may abort if it is too long).")
            shutil.rmtree(staging, ignore_errors=True)
            return rnaview_dir, None


def run_rnaview(rnaview_exe: str, rnaview_dir: str, structure_file: str,
                force_fresh: bool = False) -> str:
    """Run RNAVIEW on a PDB file. Returns path to the .out annotation, '' on failure.

    RNAVIEW writes ``<input>.out`` next to the input structure and requires the
    RNAVIEW environment variable to point at the directory containing BASEPARS.
    """
    struct_dir = os.path.dirname(os.path.abspath(structure_file)) or '.'
    base = os.path.basename(structure_file)
    # Candidate output names produced by RNAVIEW (regular and NMR variants).
    candidates = [
        structure_file + '.out',
        os.path.join(struct_dir, base + '.out'),
        os.path.splitext(structure_file)[0] + '_nmr.pdb.out',
    ]

    if force_fresh:
        for c in candidates:
            if os.path.isfile(c):
                try:
                    os.remove(c)
                except Exception as exc:
                    print(f"[rmsx_runner] WARNING: could not remove old RNAVIEW file: {c} ({exc})")
    else:
        for c in candidates:
            if os.path.isfile(c):
                print(f"[rmsx_runner] Using existing RNAVIEW annotation: {c}")
                return c

    # RNAVIEW builds BASEPARS file paths into a fixed 80-char stack buffer
    # (get_reference_pdb: char spdb[80]; sprintf(spdb, "%sAtomic_%c.pdb", BDIR, ..)).
    # A deep install path overflows that buffer and aborts the process. Expose the
    # resource directory through a short path so the unmodified binary is safe.
    effective_dir, staging_cleanup = _stage_short_rnaview_dir(rnaview_dir)

    # RNAVIEW can also abort on very long absolute input filenames. When the
    # structure path is long, run against a short staged path and move output
    # back to the canonical location afterward.
    effective_structure = structure_file
    input_staging_cleanup = None
    staged_candidates = []
    abs_structure = os.path.abspath(structure_file)
    if len(abs_structure) >= 96:
        stage_root = staging_cleanup
        if not stage_root:
            tmp_root = '/tmp' if (os.name != 'nt' and os.path.isdir('/tmp')) else tempfile.gettempdir()
            try:
                stage_root = tempfile.mkdtemp(prefix='rv', dir=tmp_root)
                input_staging_cleanup = stage_root
            except Exception as exc:
                stage_root = ''
                print(f"[rmsx_runner] WARNING: could not create RNAVIEW input staging dir ({exc}); "
                      f"using original input path (RNAVIEW may abort if it is too long).")

        if stage_root:
            staged_input = os.path.join(stage_root, 'in.pdb')
            try:
                if os.path.lexists(staged_input):
                    os.remove(staged_input)
                os.symlink(abs_structure, staged_input)
            except Exception:
                try:
                    shutil.copy2(abs_structure, staged_input)
                except Exception as exc:
                    print(f"[rmsx_runner] WARNING: could not stage RNAVIEW input ({exc}); "
                          f"using original input path (RNAVIEW may abort if it is too long).")
                    staged_input = ''
            if staged_input:
                effective_structure = staged_input
                staged_candidates = [
                    effective_structure + '.out',
                    os.path.splitext(effective_structure)[0] + '_nmr.pdb.out',
                ]

    env = dict(os.environ)
    if effective_dir:
        env['RNAVIEW'] = effective_dir

    cmd = [rnaview_exe, effective_structure]
    print(f"[rmsx_runner] Running RNAVIEW: {' '.join(cmd)}  (RNAVIEW={effective_dir or '<unset>'})")
    try:
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=300, cwd=os.path.dirname(os.path.abspath(effective_structure)) or '.', env=env)
        except Exception as exc:
            print(f"[rmsx_runner] RNAVIEW error: {exc}")
            return ''
        if result.returncode != 0:
            print(f"[rmsx_runner] RNAVIEW failed (exit {result.returncode}).")
            stderr_lines = (result.stderr or '').strip().splitlines()
            stdout_lines = (result.stdout or '').strip().splitlines()
            if stderr_lines:
                print(f"[rmsx_runner]   stderr tail: {stderr_lines[-1]}")
            if stdout_lines:
                print(f"[rmsx_runner]   stdout tail: {stdout_lines[-1]}")
            return ''

        # If the run used a staged short input path, move output back so callers
        # always observe the canonical <original_input>.out location.
        if effective_structure != structure_file:
            staged_primary = effective_structure + '.out'
            if os.path.isfile(staged_primary):
                try:
                    shutil.move(staged_primary, structure_file + '.out')
                except Exception as exc:
                    print(f"[rmsx_runner] WARNING: could not move staged RNAVIEW output "
                          f"to canonical path ({exc}).")
            staged_nmr = os.path.splitext(effective_structure)[0] + '_nmr.pdb.out'
            if os.path.isfile(staged_nmr):
                try:
                    shutil.move(staged_nmr, os.path.splitext(structure_file)[0] + '_nmr.pdb.out')
                except Exception as exc:
                    print(f"[rmsx_runner] WARNING: could not move staged RNAVIEW NMR output "
                          f"to canonical path ({exc}).")
    finally:
        if input_staging_cleanup:
            shutil.rmtree(input_staging_cleanup, ignore_errors=True)
        if staging_cleanup:
            shutil.rmtree(staging_cleanup, ignore_errors=True)

    for c in candidates:
        if os.path.isfile(c):
            print(f"[rmsx_runner] RNAVIEW annotation saved: {c}")
            return c
    print("[rmsx_runner] RNAVIEW produced no recognizable .out annotation file.")
    return ''


def run_rnaview_if_enabled(config: dict, structure_file: str,
                           force_fresh: bool = False) -> str:
    """Resolve RNAVIEW policy and, when enabled and available, produce its annotation.

    Returns the path to the RNAVIEW .out file, or '' when RNAVIEW is disabled or
    unavailable. Never silently degrades: if RNAVIEW is enabled but unavailable,
    this either aborts (when config['rnaview_required'] is true) or emits a
    prominent fallback diagnostic before returning ''.
    """
    incorporate = bool(config.get('incorporate_rnaview', True))
    required = bool(config.get('rnaview_required', False))
    if not incorporate:
        print("[rmsx_runner] RNAVIEW incorporation disabled (incorporate_rnaview=false); "
              "using MC-Annotate-only annotation by explicit configuration.")
        return ''

    exe, base_dir = _locate_rnaview(config)
    if not exe:
        print("[rmsx_runner] ============================================================")
        print("[rmsx_runner] RNAVIEW executable not found.")
        print("[rmsx_runner] Reference RNAMotifScanX preprocessing (Zhong & Zhang, 2015)")
        print("[rmsx_runner] annotates the target with BOTH MC-Annotate AND RNAVIEW and")
        print("[rmsx_runner] merges them (union; MC-Annotate wins on conflict).")
        print("[rmsx_runner] Set 'rnaview_executable' (and 'rnaview_dir' for BASEPARS) in")
        print("[rmsx_runner] the config, or place 'rnaview' on PATH.")
        print("[rmsx_runner] ============================================================")
        if required:
            raise RuntimeError(
                "RNAVIEW required (rnaview_required=true) but not found; aborting preprocessing."
            )
        print("[rmsx_runner] FALLBACK: proceeding with MC-Annotate-only annotation. "
              "Results may differ from the published pipeline.")
        return ''

    if not base_dir or not os.path.isdir(base_dir):
        print(f"[rmsx_runner] WARNING: RNAVIEW base dir (BASEPARS) not resolved: '{base_dir}'. "
              "RNAVIEW may fail without the RNAVIEW environment variable set.")

    out = run_rnaview(exe, base_dir, structure_file, force_fresh=force_fresh)
    if not out:
        if required:
            raise RuntimeError(
                "RNAVIEW required (rnaview_required=true) but execution failed; aborting preprocessing."
            )
        print("[rmsx_runner] FALLBACK: RNAVIEW unavailable/failed; proceeding with "
              "MC-Annotate-only annotation. Results may differ from the published pipeline.")
    return out


def _parse_rnaview_output(rnaview_file: str):
    """Parse RNAVIEW .out base pairs.

    Faithful port of the reference ParseStructureAnnotation.GetRNAVIEWInteractions.
    Returns a list of [nt1_id, nt2_id, edge, orientation] entries whose residue-id
    format (``<chain><resnum>``, digit chains quoted) matches the MC-Annotate token
    format so the union-merge keys align.
    """
    interactions = []
    start_parsing = False
    with open(rnaview_file, 'r', encoding='utf-8', errors='ignore') as rvw_fh:
        for raw in rvw_fh:
            line = ' ' + raw
            if re.search(r'BEGIN_base-pair', line):
                start_parsing = True
                continue
            if re.search(r'END_base-pair', line):
                break
            if not start_parsing:
                continue
            decom = re.split(r'\s+', line)
            if (len(decom) >= 9
                    and re.search(r'[HWShws+-]/[HWShws+-]', decom[7])
                    and (decom[8] == 'cis' or decom[8] == 'tran')):
                decom[2] = decom[2].rstrip(':')
                decom[6] = decom[6].rstrip(':')
                if not decom[2] == decom[6]:
                    continue
                if re.search(r'\d', decom[2]):
                    decom[2] = "'" + decom[2] + "'"
                if re.search(r'\d', decom[6]):
                    decom[6] = "'" + decom[6] + "'"
                decom[2] = decom[2] + decom[3]
                decom[6] = decom[6] + decom[5]
                if decom[7] == '+/+' or decom[7] == '-/-':
                    decom[7] = 'W/W'
                decom[7] = decom[7].upper()
                if decom[8] == 'tran':
                    decom[8] = 'trans'
                interactions.append([decom[2], decom[6], decom[7], decom[8]])
    return interactions


def _merge_interactions(mca_interactions, rvw_interactions):
    """Union of MC-Annotate and RNAVIEW interactions.

    Faithful port of the reference ParseStructureAnnotation.MergeInteractions:
    all MC-Annotate entries are kept and hashed on the (nt1, nt2) residue pair;
    an RNAVIEW pair is appended only when that pair is absent from the MC-Annotate
    set (MC-Annotate takes precedence on conflict).
    """
    merged_interactions = []
    interaction_hash = {}
    for single_interaction in mca_interactions:
        merged_interactions.append(single_interaction)
        key = single_interaction[0] + '_' + single_interaction[1]
        interaction_hash[key] = 1
    for single_interaction in rvw_interactions:
        key = single_interaction[0] + '_' + single_interaction[1]
        if key not in interaction_hash:
            merged_interactions.append(single_interaction)
    return merged_interactions


def _load_reference_sequence(reference_fasta: str, seq_tag: str) -> str:
    seq_tag = seq_tag.upper()
    if not os.path.isfile(reference_fasta):
        return ''
    with open(reference_fasta, 'r', encoding='utf-8', errors='ignore') as fh:
        lines = fh.readlines()
    for idx, line in enumerate(lines):
        hdr = line.strip().upper()
        if hdr.startswith('>') and seq_tag in hdr:
            if idx + 1 < len(lines):
                return lines[idx + 1].strip().upper()
    return ''


def _parse_mc_annotate_output(mc_annotate_file: str):
    residues = []
    interactions = []
    section = None

    with open(mc_annotate_file, 'r', encoding='utf-8', errors='ignore') as fh:
        for raw in fh:
            line = raw.rstrip('\n')

            if line.startswith('Adjacent stackings'):
                section = 'adjacent'
                continue
            if line.startswith('Non-Adjacent stackings'):
                section = 'non_adjacent'
                continue
            if line.startswith('Base-pairs'):
                section = 'base_pairs'
                continue

            if section is None:
                decom = re.split(r'\s+', line.strip())
                if len(decom) >= 3 and decom[2] in {'A', 'C', 'G', 'U'}:
                    rid = decom[0]
                    m = re.search(r"'?(\w)'?(\d+)", rid)
                    if m:
                        residues.append([rid, m.group(1), m.group(2), decom[2]])
                continue

            if ':' not in line:
                continue
            decom = re.split(r'\s+', line.strip())
            if len(decom) < 4:
                continue
            pair = decom[0].split('-', 1)
            if len(pair) != 2:
                continue
            a, b = pair[0], pair[1]

            if section == 'non_adjacent' and len(decom) >= 3 and decom[2] in {'inward', 'outward', 'upward', 'downward'}:
                interactions.append([a, b, decom[2]])
            elif section == 'adjacent' and len(decom) >= 4 and decom[3] in {'inward', 'outward', 'upward', 'downward'}:
                interactions.append([a, b, decom[3]])
            elif section == 'base_pairs':
                geom = decom[3] if len(decom) >= 4 else ''
                m = re.search(r'([HWS]).*/([HWS]).*', geom)
                if m:
                    edge = f"{m.group(1)}/{m.group(2)}"
                    orient = 'hbond'
                    if 'cis' in line:
                        orient = 'cis'
                    elif 'trans' in line:
                        orient = 'trans'
                    interactions.append([a, b, edge, orient])
                else:
                    interactions.append([a, b, '-/-', 'hbond'])

    return residues, interactions


def _global_align_map(ref_seq: str, rec_seq: str):
    """Needleman-Wunsch mapping equivalent to pairwise2.globalms(...)."""
    match_score = 3
    mismatch_score = -100
    gap_open = -10
    gap_extend = -2

    n, m = len(ref_seq), len(rec_seq)
    neg_inf = -10**12

    M = [[neg_inf] * (m + 1) for _ in range(n + 1)]
    X = [[neg_inf] * (m + 1) for _ in range(n + 1)]
    Y = [[neg_inf] * (m + 1) for _ in range(n + 1)]
    back = [[('M', 'M')] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = gap_open + (i - 1) * gap_extend
    for j in range(1, m + 1):
        Y[0][j] = gap_open + (j - 1) * gap_extend

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match_score if ref_seq[i - 1] == rec_seq[j - 1] else mismatch_score
            prev_vals = [(M[i - 1][j - 1], 'M'), (X[i - 1][j - 1], 'X'), (Y[i - 1][j - 1], 'Y')]
            best_prev = max(prev_vals, key=lambda t: t[0])
            M[i][j] = best_prev[0] + s

            open_x = M[i - 1][j] + gap_open
            ext_x = X[i - 1][j] + gap_extend
            if open_x >= ext_x:
                X[i][j] = open_x
            else:
                X[i][j] = ext_x

            open_y = M[i][j - 1] + gap_open
            ext_y = Y[i][j - 1] + gap_extend
            if open_y >= ext_y:
                Y[i][j] = open_y
            else:
                Y[i][j] = ext_y

            best_state = max([(M[i][j], 'M'), (X[i][j], 'X'), (Y[i][j], 'Y')], key=lambda t: t[0])[1]
            back[i][j] = (best_state, best_prev[1])

    state = max([(M[n][m], 'M'), (X[n][m], 'X'), (Y[n][m], 'Y')], key=lambda t: t[0])[1]
    i, j = n, m
    map_rec_to_ref = {}
    while i > 0 or j > 0:
        if state == 'M':
            if i > 0 and j > 0:
                map_rec_to_ref[j - 1] = i - 1
                prev_state = back[i][j][1]
                i -= 1
                j -= 1
                state = prev_state
            elif i > 0:
                i -= 1
                state = 'X'
            else:
                j -= 1
                state = 'Y'
        elif state == 'X':
            if i > 0 and X[i][j] == X[i - 1][j] + gap_extend:
                i -= 1
                state = 'X'
            else:
                i -= 1
                state = 'M'
        else:
            if j > 0 and Y[i][j] == Y[i][j - 1] + gap_extend:
                j -= 1
                state = 'Y'
            else:
                j -= 1
                state = 'M'

    return map_rec_to_ref


def prepare_rmsx_inputs_from_annotation(annotation_file: str, pdb_id: str,
                                        output_dir: str, chains: list[str],
                                        reference_fasta: str,
                                        rnaview_file: str = '') -> dict[str, str]:
    """Generate per-chain .rmsx.in files from the annotation.

    Reproduces the upstream PrepareInput.py behaviour: when an RNAVIEW annotation
    is supplied, the MC-Annotate and RNAVIEW interaction sets are merged (union,
    MC-Annotate precedence) before the .rmsx.in target file is written.
    """
    residues, interactions = _parse_mc_annotate_output(annotation_file)
    if not residues:
        return {}

    if rnaview_file and os.path.isfile(rnaview_file):
        rvw_interactions = _parse_rnaview_output(rnaview_file)
        before = len(interactions)
        interactions = _merge_interactions(interactions, rvw_interactions)
        added = len(interactions) - before
        print(f"[rmsx_runner] Merged annotations: MC-Annotate={before}, "
              f"RNAVIEW={len(rvw_interactions)}, +{added} unique from RNAVIEW "
              f"(union, MC-Annotate precedence)")
    else:
        print("[rmsx_runner] Annotation source: MC-Annotate only "
              "(no RNAVIEW annotation merged)")

    residues_by_chain = {}
    for rid, chain, idx, nuc in residues:
        residues_by_chain.setdefault(chain, []).append((rid, int(idx), nuc))

    wanted = set(chains or [])
    produced = {}
    pdb_upper = pdb_id.upper()
    for chain, items in residues_by_chain.items():
        if wanted and chain not in wanted:
            continue
        items = sorted(items, key=lambda t: t[1])
        rec_seq = ''.join(nuc for _, _, nuc in items)
        rec_ids = [rid for rid, _, _ in items]

        seq_tag = f"{pdb_upper}_{chain}"
        ref_seq = _load_reference_sequence(reference_fasta, seq_tag)
        if not ref_seq:
            # Fallback keeps pipeline running when seqref entry is unavailable.
            ref_seq = rec_seq
            map_rec_to_ref = {i: i for i in range(len(rec_ids))}
        else:
            map_rec_to_ref = _global_align_map(ref_seq, rec_seq)

        nucleotide_hash = {}
        for rec_i, rid in enumerate(rec_ids):
            if rec_i in map_rec_to_ref:
                nucleotide_hash[rid] = map_rec_to_ref[rec_i]

        out_file = os.path.join(output_dir, f"{seq_tag}.rmsx.in")
        with open(out_file, 'w', encoding='utf-8') as out_fh:
            out_fh.write(f">{seq_tag}\n")
            out_fh.write(f"{ref_seq}\n")
            out_fh.write("#info=basepair\n")
            for it in interactions:
                if len(it) == 4 and it[0] in nucleotide_hash and it[1] in nucleotide_hash:
                    i = nucleotide_hash[it[0]]
                    j = nucleotide_hash[it[1]]
                    out_fh.write(f"{i}-{j},{it[2]},{it[3]},{it[0]}-{it[1]}\n")
            out_fh.write("#info=stacking\n")
            for it in interactions:
                if len(it) == 3 and it[0] in nucleotide_hash and it[1] in nucleotide_hash:
                    i = nucleotide_hash[it[0]]
                    j = nucleotide_hash[it[1]]
                    out_fh.write(f"{i}-{j},{it[2]},{it[0]}-{it[1]}\n")

        produced[chain] = out_file
        print(f"[rmsx_runner] Prepared RNAMotifScanX structure: {out_file}")

    return produced


# ─────────────────────────────────────────────────────────────────────────────
# RMSX execution
# ─────────────────────────────────────────────────────────────────────────────

def run_rmsx_for_family(rmsx_exe: str, query_file: str, annotation_file: str,
                        out_dir: str, pdb_id: str, chains: list,
                        max_strands: int, num_threads: int,
                        force_fresh: bool = False) -> bool:
    """Run RMSX for one motif family. Returns True on success."""
    os.makedirs(out_dir, exist_ok=True)
    out_log = os.path.join(out_dir, 'result_0_100_withbs.log')

    if force_fresh and os.path.isfile(out_log):
        try:
            os.remove(out_log)
            print(f"[rmsx_runner]   Removed old log for fresh run: {out_log}")
        except Exception as exc:
            print(f"[rmsx_runner]   WARNING: could not remove old log: {out_log} ({exc})")

    if os.path.isfile(out_log):
        if _result_file_has_pdb_hits(out_log, pdb_id):
            print(f"[rmsx_runner]   Already exists for {pdb_id}: {out_log}")
            return True
        print(f"[rmsx_runner]   Existing log does not match {pdb_id}; regenerating: {out_log}")
        try:
            os.remove(out_log)
        except Exception as exc:
            print(f"[rmsx_runner]   WARNING: could not remove stale log: {out_log} ({exc})")

    if not os.path.isfile(rmsx_exe):
        print(f"[rmsx_runner]   ERROR: RMSX executable not found: {rmsx_exe}")
        return False

    chain_str = ','.join(chains)

    # Prefer modern scan CLI (boost::program_options long options), and keep
    # legacy short-option invocation as fallback for older wrappers.
    cmd_modern = [
        rmsx_exe,
        '--query_motif', query_file,
        '--structure', annotation_file,
        '--max_num_strands', str(max_strands),
        '--num_threads', str(num_threads),
    ]
    cmd_legacy = [
        rmsx_exe,
        '-q', query_file,
        '-t', annotation_file,
        '-c', chain_str,
        '-s', str(max_strands),
        '-p', str(num_threads),
        '-o', out_dir,
    ]

    env = os.environ.copy()
    runtime_root = os.path.abspath(os.path.join(os.path.dirname(query_file), os.pardir))
    if os.path.isdir(os.path.join(runtime_root, 'mat')):
        env['RNAMOTIFSCANX_PATH'] = runtime_root

    def _run_and_capture(cmd: list[str], write_stdout_log: bool) -> tuple[bool, int, str]:
        print(f"[rmsx_runner]   Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, env=env)
        if result.returncode != 0:
            return False, result.returncode, (result.stderr or result.stdout or '')[-1000:]
        if write_stdout_log:
            with open(out_log, 'w', encoding='utf-8') as fh:
                fh.write(result.stdout or '')
        return True, 0, ''

    try:
        ok, exit_code, err = _run_and_capture(cmd_modern, write_stdout_log=True)
        if not ok:
            print(f"[rmsx_runner]   Modern CLI failed (exit {exit_code}); trying legacy flags...")
            if err:
                print(err)
            ok, exit_code, err = _run_and_capture(cmd_legacy, write_stdout_log=False)
        if not ok:
            print(f"[rmsx_runner]   RMSX failed (exit {exit_code})")
            if err:
                print(err)
            return False
        print(f"[rmsx_runner]   Done → {out_log}")
        return True
    except subprocess.TimeoutExpired:
        print("[rmsx_runner]   ERROR: RMSX timed out (>1 hour)")
        return False
    except Exception as exc:
        print(f"[rmsx_runner]   RMSX error: {exc}")
        return False


# ─────────────────────────────────────────────────────────────────────────────
# Platform check
# ─────────────────────────────────────────────────────────────────────────────

def _check_platform_compatibility(rmsx_exe: str) -> bool:
    """Best-effort compatibility check; allow execution when an executable is present."""
    if not rmsx_exe or not os.path.isfile(rmsx_exe):
        return False
    return True


def _result_file_has_pdb_hits(log_file: str, pdb_id: str) -> bool:
    """Return True if RMSX result file contains at least one hit for pdb_id."""
    pdb_prefix = f"{str(pdb_id).strip().upper()}_"
    try:
        with open(log_file, 'r', encoding='utf-8', errors='ignore') as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith('#'):
                    continue
                # First tab-delimited token is fragment_ID, e.g. 1S72_0:75-85_89-98_58-60
                fragment_id = line.split('\t', 1)[0].strip()
                if fragment_id.upper().startswith(pdb_prefix):
                    return True
    except Exception:
        return False
    return False


# ─────────────────────────────────────────────────────────────────────────────
# High-level API (called by RSMViewer gui.py)
# ─────────────────────────────────────────────────────────────────────────────

def check_results_exist(config: dict, pdb_id: str) -> dict:
    """Return {family: path} where the family log has at least one hit for pdb_id."""
    output_dir = os.path.expanduser(str(config.get('output_dir', '') or ''))
    query_file = str(config.get('query_file', '') or '').strip()
    if query_file:
        families = [Path(query_file).stem]
    else:
        families = config.get('motif_families', list(DEFAULT_FAMILY_FOLDER_MAP.keys()))
    found = {}
    for family in families:
        folder = DEFAULT_FAMILY_FOLDER_MAP.get(family, family if family.endswith('_consensus') else f'{family}_consensus')
        log_file = os.path.join(output_dir, folder, 'result_0_100_withbs.log')
        if os.path.isfile(log_file) and _result_file_has_pdb_hits(log_file, pdb_id):
            found[family] = log_file
    return found


def run_pipeline(config: dict, pdb_id: str, cif_file: str = '', force_fresh: bool = False) -> dict:
    """Run the full RMSX pipeline for a PDB.

    Returns {family: result_log_path} for successfully produced families.
    Empty dict on complete failure.
    """
    pdb_id = pdb_id.strip().upper()
    rmsx_exe = os.path.expanduser(str(config.get('rmsx_executable', '') or ''))
    mc_exe   = os.path.expanduser(str(config.get('mc_annotate_executable', '') or ''))
    query_dir = os.path.expanduser(str(config.get('query_motifs_dir', '') or ''))
    query_file = os.path.expanduser(str(config.get('query_file', '') or ''))
    output_dir = os.path.expanduser(str(config.get('output_dir', '.') or '.'))
    cif_in_dir = os.path.expanduser(str(config.get('cif_input_dir', '') or ''))
    auto_dl    = bool(config.get('auto_download_cif', True))
    auto_dl_pdb = bool(config.get('auto_download_pdb', True))
    families   = config.get('motif_families', list(DEFAULT_FAMILY_FOLDER_MAP.keys()))
    chains     = config.get('target_chains', ['0'])
    max_str    = int(config.get('max_strands', 3))
    threads    = int(config.get('num_threads', 4))
    seq_ref = os.path.expanduser(str(config.get(
        'reference_sequence_fasta',
        os.path.join(os.path.dirname(__file__), 'rmsx_runtime', 'src', 'RNAMotifScanX_src', 'StructureAnnotation', 'pdb_seqres.na.fa')
    ) or ''))

    os.makedirs(output_dir, exist_ok=True)

    # ── Check if results already exist ────────────────────────────────────
    existing = {} if force_fresh else check_results_exist(config, pdb_id)
    if (not force_fresh) and len(existing) == len(families):
        print(f"[rmsx_runner] All {len(existing)} family result files already exist — skipping run")
        return existing
    if force_fresh:
        print(f"[rmsx_runner] Fresh run requested for {pdb_id}: ignoring cached RMSX results")

    # ── Platform + executable check ───────────────────────────────────────
    if not _check_platform_compatibility(rmsx_exe):
        if existing:
            print(f"[rmsx_runner] {len(existing)}/{len(families)} result files found; "
                  f"using available files.")
            return existing
        if not rmsx_exe or not os.path.isfile(rmsx_exe):
            print("[rmsx_runner] No result files found and no runnable executable is configured on this platform.")
            print("[rmsx_runner] To generate results on Linux:")
            print(f"[rmsx_runner]   1. Run: python rmsx_runner.py --config <config> --pdb {pdb_id}")
            print(f"[rmsx_runner]   2. Copy {output_dir}/* to this machine")
            return {}
        print("[rmsx_runner] Attempting execution using configured executable/wrapper on this platform...")

    if not rmsx_exe or not os.path.isfile(rmsx_exe):
        print(f"[rmsx_runner] ERROR: RMSX executable not found: {rmsx_exe}")
        return existing if existing else {}

    prebuild_archive = os.path.expanduser(str(config.get('pdb_prebuild_archive', '') or ''))
    prepared_targets = {}
    annot_file = ''
    pdb_file = ''

    # Prefer prebuilt targets when available (unless fresh rerun is requested).
    if not force_fresh:
        prepared_targets = _extract_prebuilt_targets_from_archive(
            prebuild_archive, pdb_id, output_dir, chains
        )

    if not prepared_targets:
        # ── Find CIF ──────────────────────────────────────────────────────
        if not cif_file or not os.path.isfile(cif_file):
            pdb_lower = pdb_id.lower()
            for sd in [d for d in [cif_in_dir, output_dir, '.'] if d]:
                for name in [f'{pdb_lower}.cif', f'{pdb_id}.cif',
                             f'{pdb_lower}.cif.gz', f'{pdb_id}.cif.gz']:
                    candidate = os.path.join(sd, name)
                    if os.path.isfile(candidate):
                        cif_file = candidate
                        break
                if cif_file:
                    break
        if not cif_file or not os.path.isfile(cif_file):
            if auto_dl:
                cif_file = _download_cif(pdb_id, output_dir)
        if not cif_file or not os.path.isfile(cif_file):
            print(f"[rmsx_runner] ERROR: CIF file not found for {pdb_id}")
            return existing if existing else {}

        print(f"[rmsx_runner] CIF: {cif_file}")

        # ── Prepare PDB for MC-Annotate (original workflow expects PDB-like input) ──
        for sd in [d for d in [cif_in_dir, output_dir, '.'] if d]:
            for name in [f'{pdb_id}.pdb', f'{pdb_id.lower()}.pdb']:
                candidate = os.path.join(sd, name)
                if os.path.isfile(candidate):
                    pdb_file = candidate
                    break
            if pdb_file:
                break
        if not pdb_file and auto_dl_pdb:
            pdb_file = _download_pdb(pdb_id, output_dir)
        if not pdb_file:
            print(f"[rmsx_runner] ERROR: PDB file not found for {pdb_id}; cannot run MC-Annotate preparation step")
            return existing if existing else {}
        print(f"[rmsx_runner] PDB for annotation: {pdb_file}")

        # ── Annotate (MC-Annotate + RNAVIEW → union merge → .rmsx.in) ─────
        annot_file = run_mc_annotate(mc_exe, pdb_file, output_dir, pdb_id, force_fresh=force_fresh)
        if not annot_file:
            print("[rmsx_runner] WARNING: annotation step failed; RMSX search may not work.")

        if annot_file:
            try:
                rnaview_file = run_rnaview_if_enabled(config, pdb_file, force_fresh=force_fresh)
            except RuntimeError as exc:
                print(f"[rmsx_runner] ERROR: {exc}")
                return existing if existing else {}
            prepared_targets = prepare_rmsx_inputs_from_annotation(
                annot_file, pdb_id, output_dir, chains, seq_ref, rnaview_file
            )

    if not prepared_targets and prebuild_archive:
        prepared_targets = _extract_prebuilt_targets_from_archive(
            prebuild_archive, pdb_id, output_dir, chains
        )

    if not prepared_targets:
        print("[rmsx_runner] WARNING: preparation step produced no .rmsx.in target files")

    # ── Resolve query mode ────────────────────────────────────────────────
    if query_file:
        if not os.path.isfile(query_file):
            print(f"[rmsx_runner] ERROR: query file not found: {query_file}")
            return existing if existing else {}
        families = [Path(query_file).stem]

    # ── Run RMSX per family ───────────────────────────────────────────────
    results = dict(existing)
    for family in families:
        if family in results:
            continue  # already exists
        folder = DEFAULT_FAMILY_FOLDER_MAP.get(family, family if family.endswith('_consensus') else f'{family}_consensus')
        family_out_dir = os.path.join(output_dir, folder)

        # Find query consensus file
        family_query_file = query_file
        if not family_query_file and query_dir and os.path.isdir(query_dir):
            for name in [
                f'{family}_consensus.struct', f'{family}.struct', f'{folder}.struct', f'{family}_query.struct',
                f'{family}_consensus.txt', f'{family}.txt', f'{folder}.txt', f'{family}_query.txt'
            ]:
                candidate = os.path.join(query_dir, name)
                if os.path.isfile(candidate):
                    family_query_file = candidate
                    break

        if not family_query_file:
            print(f"[rmsx_runner] WARNING: No query file found for '{family}' in {query_dir}")
            print(f"[rmsx_runner]   Expected: {family}_consensus.struct or {family}_consensus.txt (in query_motifs_dir)")
            continue

        print(f"[rmsx_runner] Running RMSX for motif family: {family}")
        chain_id = str(chains[0]) if chains else ''
        target_for_scan = prepared_targets.get(chain_id)
        if not target_for_scan:
            if prepared_targets:
                target_for_scan = next(iter(prepared_targets.values()))
            else:
                target_for_scan = annot_file or pdb_file or cif_file

        ok = run_rmsx_for_family(
            rmsx_exe, family_query_file, target_for_scan,
            family_out_dir, pdb_id, chains, max_str, threads,
            force_fresh=force_fresh
        )
        if ok:
            log = os.path.join(family_out_dir, 'result_0_100_withbs.log')
            if os.path.isfile(log) and _result_file_has_pdb_hits(log, pdb_id):
                results[family] = log
            elif os.path.isfile(log):
                print(f"[rmsx_runner]   WARNING: {family} log produced but has no hits for {pdb_id}; skipping")

    if results:
        print(f"[rmsx_runner] Done: {len(results)}/{len(families)} families available")
    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Run RNAMotifScanX pipeline for a PDB structure'
    )
    parser.add_argument('--config', required=True,
                        help='Path to rmsx_pipeline_config.json')
    parser.add_argument('--pdb', required=True,
                        help='PDB ID (e.g. 1S72)')
    parser.add_argument('--cif', default='',
                        help='Explicit path to CIF file (optional)')
    parser.add_argument('--check', action='store_true',
                        help='Only check which result files exist, do not run')
    parser.add_argument('--fresh', action='store_true',
                        help='Force fresh execution (ignore/remove cached outputs)')
    args = parser.parse_args()

    with open(args.config, 'r', encoding='utf-8') as fh:
        config = json.load(fh)

    if args.check:
        existing = check_results_exist(config, args.pdb)
        families = config.get('motif_families', list(DEFAULT_FAMILY_FOLDER_MAP.keys()))
        print(f"Result files for {args.pdb.upper()} ({len(existing)}/{len(families)} families):")
        for family in families:
            path = existing.get(family, '(missing)')
            marker = '✓' if family in existing else '✗'
            print(f"  {marker} {family:<20} {path}")
        return 0

    results = run_pipeline(config, args.pdb, args.cif, force_fresh=args.fresh)
    if not results:
        print("[rmsx_runner] Pipeline produced no results.")
        return 1

    print(f"\nResult files ({len(results)}):")
    for family, path in results.items():
        print(f"  {family}: {path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
