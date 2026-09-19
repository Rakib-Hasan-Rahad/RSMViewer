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
    rmv_fetch 1S72
    rmv_db RNAMotifScanX

Usage (standalone):
    python rmsx_runner.py --config config/rmsx_config.json --pdb 1S72
    python rmsx_runner.py --config config/rmsx_config.json --pdb 1S72 --chains 0

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
import ssl
import subprocess
import sys
import tempfile
import tarfile
import time
import uuid
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


def _copy_prebuilt_targets_from_directory(prebuild_dir: str, pdb_id: str,
                                          output_dir: str, chains: list) -> dict:
    """Copy prepared targets from expanded ``rmsx_work_default`` data."""
    root = Path(os.path.expanduser(str(prebuild_dir or ""))).resolve()
    entry = root / str(pdb_id).lower()
    if not entry.is_dir():
        entry = root / str(pdb_id).upper()
    if not entry.is_dir():
        return {}

    wanted = {str(chain).strip() for chain in (chains or []) if str(chain).strip()}
    accept_any = not wanted or "0" in wanted
    candidates = []
    for path in entry.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in (".in", ".nch"):
            continue
        prefix = f"{str(pdb_id).lower()}_"
        if not path.name.lower().startswith(prefix):
            continue
        chain = path.name[len(str(pdb_id)):].split(".", 1)[0].lstrip("_")
        if not accept_any and chain not in wanted:
            continue
        priority = 0 if path.suffix.lower() == ".in" else 1
        candidates.append((priority, str(path), chain))

    produced = {}
    stage_dir = Path(output_dir) / "_prebuilt_targets"
    stage_dir.mkdir(parents=True, exist_ok=True)
    for _priority, source, chain in sorted(candidates):
        if chain in produced:
            continue
        destination = stage_dir / Path(source).name
        destination.write_bytes(Path(source).read_bytes())
        produced[chain] = str(destination)
    if produced:
        print(f"[rmsx_runner] Using expanded prebuilt RMSX targets: {entry}")
    return produced


def copy_preannotated_results(prebuild_path: str, pdb_id: str, output_dir: str) -> dict:
    """Copy preannotated RMSX consensus logs into the normal output cache.

    A single PDB may carry per-chain consensus logs for the same motif family
    (e.g. ``1s72/0/sarcin-ricin_consensus.log`` and
    ``1s72/9/sarcin-ricin_consensus.log``). All chains for a given family are
    concatenated into one ``result_0_100_withbs.log`` so no chain's hits are
    lost to overwriting.
    """
    source = Path(os.path.expanduser(str(prebuild_path or ""))).resolve()
    pdb_lower = str(pdb_id).lower()

    # Enumerating members of a large gzip tar requires decompressing the whole
    # stream, so each PDB's small consensus logs are extracted once and reused
    # on later loads and across PyMOL restarts. The cache is stamped with the
    # source identity so a changed archive/dir invalidates it automatically.
    cache_root = Path(output_dir) / ".preannotated_cache" / pdb_lower
    marker = cache_root / ".source"

    def _source_stamp() -> str:
        try:
            stat = source.stat()
            return f"{source}|{int(stat.st_mtime)}|{stat.st_size}"
        except OSError:
            return str(source)

    def _write_output(family_bytes: dict) -> int:
        written = 0
        for family, blob in family_bytes.items():
            destination = Path(output_dir) / family
            destination.mkdir(parents=True, exist_ok=True)
            (destination / "result_0_100_withbs.log").write_bytes(blob)
            written += 1
        return written

    stamp = _source_stamp()

    if marker.is_file():
        try:
            cache_valid = marker.read_text(encoding="utf-8").strip() == stamp
        except OSError:
            cache_valid = False
        if cache_valid:
            family_bytes = {
                cached.stem: cached.read_bytes()
                for cached in sorted(cache_root.glob("*.log"))
            }
            if family_bytes:
                copied = _write_output(family_bytes)
                return {"copied": copied, "output_dir": str(output_dir), "cached": True}

    family_chunks: dict = {}
    if source.is_file() and tarfile.is_tarfile(source):
        # Single forward pass: reading each matching member while iterating
        # avoids the second full decompression that getmembers() would add.
        with tarfile.open(source, "r:*") as archive:
            for member in archive:
                if not member.isfile():
                    continue
                name_lower = member.name.lower()
                if f"/{pdb_lower}/" not in f"/{name_lower}":
                    continue
                if not name_lower.endswith("_consensus.log"):
                    continue
                extracted = archive.extractfile(member)
                if extracted is not None:
                    family = Path(member.name).stem
                    family_chunks.setdefault(family, []).append(extracted.read())
    else:
        candidates = []
        if source.is_dir():
            entry = source / pdb_lower
            if not entry.is_dir():
                entry = source / str(pdb_id).upper()
            if entry.is_dir():
                candidates = list(entry.glob("*/*_consensus.log"))
        for candidate in candidates:
            family_chunks.setdefault(candidate.stem, []).append(candidate.read_bytes())

    family_bytes = {family: b"\n".join(chunks) for family, chunks in family_chunks.items()}

    if family_bytes:
        try:
            cache_root.mkdir(parents=True, exist_ok=True)
            for stale in cache_root.glob("*.log"):
                stale.unlink()
            for family, blob in family_bytes.items():
                (cache_root / f"{family}.log").write_bytes(blob)
            marker.write_text(stamp, encoding="utf-8")
        except OSError:
            pass

    copied = _write_output(family_bytes)
    return {"copied": copied, "output_dir": str(output_dir), "cached": False}


DEFAULT_PREANNOTATED_BASE_URL = (
    "https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default"
)


def download_preannotated_pdb(base_url: str, pdb_id: str, dest_dir: str,
                              timeout: int = 60) -> dict:
    """Download ``<base_url>/<pdb_lower>.tar.gz`` and extract it under ``dest_dir``.

    The archive is expected to contain a top-level ``<pdb_lower>/`` folder, i.e.
    the same layout as ``rmsx_work_default/<pdb_lower>/``; if it does not, its
    contents are placed under ``dest_dir/<pdb_lower>/`` instead. Extraction is
    staged in a temporary directory and moved into place only when complete, so
    an interrupted download never leaves a partial ``<pdb_lower>/`` folder that
    would later be mistaken for real local data.

    Returns ``{"ok": bool, "url": str, "path": str, "error": str}``; never raises.
    """
    import urllib.error
    import urllib.request

    pdb_lower = str(pdb_id or "").strip().lower()
    base = str(base_url or DEFAULT_PREANNOTATED_BASE_URL).strip().rstrip("/")
    result = {"ok": False, "url": "", "path": "", "error": "", "insecure_tls": False}
    if not re.fullmatch(r"[0-9a-z]{4}", pdb_lower):
        result["error"] = f"not a 4-character PDB ID: {pdb_id!r}"
        return result

    url = f"{base}/{pdb_lower}.tar.gz"
    result["url"] = url
    dest_root = Path(os.path.expanduser(str(dest_dir))).resolve()
    final_dir = dest_root / pdb_lower

    tmp_root = None
    try:
        dest_root.mkdir(parents=True, exist_ok=True)
        tmp_root = Path(tempfile.mkdtemp(prefix=f".dl_{pdb_lower}_", dir=str(dest_root)))
        archive_path = tmp_root / f"{pdb_lower}.tar.gz"

        request = urllib.request.Request(url, headers={"User-Agent": "RSMViewer"})
        # Verified TLS first. Some Python installs (notably python.org builds on
        # macOS) ship no CA bundle and reject valid certificates, so retry with
        # certifi's bundle, and only as a last resort without verification.
        contexts = [("verified", ssl.create_default_context())]
        try:
            import certifi
            contexts.append(("verified", ssl.create_default_context(cafile=certifi.where())))
        except ImportError:
            pass
        unverified = ssl.create_default_context()
        unverified.check_hostname = False
        unverified.verify_mode = ssl.CERT_NONE
        contexts.append(("unverified", unverified))

        for index, (mode, context) in enumerate(contexts):
            try:
                with urllib.request.urlopen(request, timeout=timeout, context=context) as response, \
                        open(archive_path, "wb") as out_fh:
                    shutil.copyfileobj(response, out_fh)
            except urllib.error.URLError as exc:
                is_cert_error = isinstance(exc.reason, ssl.SSLCertVerificationError)
                if is_cert_error and index < len(contexts) - 1:
                    continue
                raise
            result["insecure_tls"] = (mode == "unverified")
            break

        extract_dir = tmp_root / "extracted"
        extract_dir.mkdir()
        with tarfile.open(archive_path, "r:*") as archive:
            members = archive.getmembers()
            for member in members:
                target = (extract_dir / member.name).resolve()
                if extract_dir.resolve() not in target.parents and target != extract_dir.resolve():
                    raise ValueError(f"unsafe path in archive: {member.name}")
                if member.issym() or member.islnk() or member.isdev():
                    raise ValueError(f"unsupported member type in archive: {member.name}")
            archive.extractall(extract_dir, members=members)

        top_level = {Path(m.name).parts[0] for m in members if Path(m.name).parts}
        staged = extract_dir / pdb_lower if top_level == {pdb_lower} else extract_dir
        if not any(staged.rglob("*_consensus.log")):
            raise ValueError("archive contains no *_consensus.log result files")

        if final_dir.exists():
            # Merge rather than replace: an existing folder may hold the user's
            # own prepared .rmsx.in/.nch inputs, which must not be deleted.
            shutil.copytree(staged, final_dir, dirs_exist_ok=True)
        else:
            shutil.move(str(staged), str(final_dir))
        result["ok"] = True
        result["path"] = str(final_dir)
    except urllib.error.HTTPError as exc:
        result["error"] = f"HTTP {exc.code} {exc.reason}"
    except Exception as exc:  # noqa: BLE001 - reported to the caller, never raised
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        if tmp_root is not None:
            shutil.rmtree(tmp_root, ignore_errors=True)
    return result


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
        '--pvalue', '1.0',
    ]
    cmd_legacy = [
        rmsx_exe,
        '-q', query_file,
        '-t', annotation_file,
        '-c', chain_str,
        '-s', str(max_strands),
        '-p', str(num_threads),
        '-v', '1.0',
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
    prebuild_dir = os.path.expanduser(str(config.get('pdb_prebuild_dir', '') or ''))
    prepared_targets = {}
    annot_file = ''
    pdb_file = ''

    # Prefer prebuilt targets when available (unless fresh rerun is requested).
    # --fresh means rerun scan outputs; it does not discard valid prepared
    # .rmsx.in/.nch inputs supplied by the prebuilt-data directory/archive.
    prepared_targets = _copy_prebuilt_targets_from_directory(
        prebuild_dir, pdb_id, output_dir, chains
    )
    if not prepared_targets:
        prepared_targets = _extract_prebuilt_targets_from_archive(
            prebuild_archive, pdb_id, output_dir, chains
        )

    if not prepared_targets:
        # run_from_scratch consumes user-provided inputs only; RSMViewer never
        # runs MC-Annotate/RNAVIEW itself. Generate the inputs externally and
        # place them under pdb_prebuild_dir/<pdb>/<chain>/<pdb>_<chain>.rmsx.{in,nch}.
        expected = os.path.join(prebuild_dir or '<pdb_prebuild_dir>', pdb_id.lower(), '<chain>')
        print(f"[rmsx_runner] ERROR: no prepared RMSX inputs (.rmsx.in/.rmsx.nch) found for {pdb_id}.")
        print("[rmsx_runner] run_from_scratch needs externally generated inputs (RNAVIEW + MC-Annotate).")
        print(f"[rmsx_runner]   Expected: {expected}/{pdb_id.lower()}_<chain>.rmsx.in and .rmsx.nch")
        return existing if existing else {}

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
        scan_targets = list(prepared_targets.items())
        if not scan_targets:
            scan_targets = [(str(chains[0]) if chains else '', annot_file or pdb_file or cif_file)]

        chain_logs = []
        for chain_id, target_for_scan in scan_targets:
            chain_out_dir = os.path.join(family_out_dir, f"_chain_{chain_id or 'default'}")
            ok = run_rmsx_for_family(
                rmsx_exe, family_query_file, target_for_scan,
                chain_out_dir, pdb_id, [str(chain_id)] if chain_id else chains,
                max_str, threads, force_fresh=force_fresh
            )
            chain_log = os.path.join(chain_out_dir, 'result_0_100_withbs.log')
            if ok and os.path.isfile(chain_log):
                chain_logs.append(chain_log)

        log = os.path.join(family_out_dir, 'result_0_100_withbs.log')
        if chain_logs:
            with open(log, 'wb') as output:
                for chain_log in chain_logs:
                    with open(chain_log, 'rb') as source:
                        output.write(source.read())
                        output.write(b'\n')
            for chain_log in chain_logs:
                shutil.rmtree(os.path.dirname(chain_log), ignore_errors=True)
            if _result_file_has_pdb_hits(log, pdb_id):
                results[family] = log
            else:
                print(f"[rmsx_runner]   WARNING: {family} log produced but has no hits for {pdb_id}; skipping")

    if results:
        print(f"[rmsx_runner] Done: {len(results)}/{len(families)} families available")
    return results


# ─────────────────────────────────────────────────────────────────────────────
# scan_prepared mode
#
# Minimal prepared-input workflow: run the RNAMotifScanX ``scan`` executable
# directly against locally prepared ``.rmsx.in``/``.rmsx.nch`` inputs, skipping
# MC-Annotate and RNAVIEW. The command mirrors the one recorded in the header of
# the distributed preannotated logs:
#
#   scan <query.struct> <target.rmsx.in> --map_pdb=<target.rmsx.nch> \
#        --pvalue 1.0 --num_threads <N> --write_alignment
#
# Preannotated results are never read, substituted, or used as a fallback in
# this mode.
# ─────────────────────────────────────────────────────────────────────────────

SCAN_PREPARED_PVALUE = '1.0'


def locate_prepared_inputs(prebuild_dir: str, pdb_id: str, chains=None):
    """Locate prepared ``.rmsx.in``/``.rmsx.nch`` pairs for ``pdb_id``.

    Returns ``(pairs, problems)`` where ``pairs`` maps ``chain -> {'in','nch'}``
    for readable, matched pairs only, and ``problems`` lists every missing or
    unreadable input encountered. Nothing is substituted for a missing pair.
    """
    problems: list[str] = []
    root = Path(os.path.expanduser(str(prebuild_dir or ''))).resolve()
    pdb_lower = str(pdb_id).strip().lower()

    entry = root / pdb_lower
    if not entry.is_dir():
        entry = root / str(pdb_id).strip().upper()
    if not entry.is_dir():
        problems.append(f"No prepared-input directory for {pdb_id} under {root}")
        return {}, problems

    wanted = {str(c).strip() for c in (chains or []) if str(c).strip()}
    accept_any = not wanted  # '0' is a real chain label, never a wildcard here.

    pairs: dict = {}
    # Prepared chains live in per-chain subdirectories (e.g. {pdb}/A, {pdb}/0);
    # ``_prep*`` staging folders are intermediate and are not scanned directly.
    chain_dirs = sorted(
        d for d in entry.iterdir() if d.is_dir() and not d.name.startswith('_')
    )
    for chain_dir in chain_dirs:
        chain = chain_dir.name
        if not accept_any and chain not in wanted:
            continue
        in_files = sorted(chain_dir.glob('*.rmsx.in'))
        nch_files = sorted(chain_dir.glob('*.rmsx.nch'))
        if not in_files:
            problems.append(f"chain {chain}: missing .rmsx.in in {chain_dir}")
            continue
        if not nch_files:
            problems.append(f"chain {chain}: missing .rmsx.nch in {chain_dir}")
            continue
        in_path, nch_path = in_files[0], nch_files[0]
        if not os.access(in_path, os.R_OK):
            problems.append(f"chain {chain}: unreadable {in_path}")
            continue
        if not os.access(nch_path, os.R_OK):
            problems.append(f"chain {chain}: unreadable {nch_path}")
            continue
        pairs[chain] = {'in': str(in_path), 'nch': str(nch_path)}

    if wanted:
        for chain in sorted(wanted):
            if chain not in pairs and not any(f"chain {chain}:" in p for p in problems):
                problems.append(f"chain {chain}: no prepared-input directory under {entry}")
    if not pairs and not problems:
        problems.append(f"No prepared per-chain inputs found under {entry}")
    return pairs, problems


def _executable_runs_natively(exe: str) -> bool:
    """Return True when ``exe`` can be exec'd directly on the current host."""
    try:
        with open(exe, 'rb') as fh:
            magic = fh.read(4)
    except OSError:
        return False
    system = platform.system()
    machine = platform.machine().lower()
    if magic[:2] == b'#!':
        return True
    if magic == b'\x7fELF':
        return system == 'Linux' and machine in ('x86_64', 'amd64')
    if magic[:2] == b'MZ':
        return system == 'Windows'
    if magic in (b'\xcf\xfa\xed\xfe', b'\xce\xfa\xed\xfe',
                 b'\xca\xfe\xba\xbe', b'\xbe\xba\xfe\xca'):
        return system == 'Darwin'
    return False


def _find_scan_docker_wrapper(config: dict) -> str:
    """Return a runnable ``scan_docker_*.sh`` wrapper path, or '' when absent."""
    exe = os.path.expanduser(str(config.get('rmsx_executable', '') or '').strip())
    search_dirs = []
    if exe:
        search_dirs.append(Path(exe).resolve().parent)
    for base in list(search_dirs):
        for name in ('scan_docker_x86_64.sh', 'scan_docker.sh'):
            candidate = base / name
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
    return ''


def resolve_scan_command(config: dict):
    """Resolve how to invoke ``scan``.

    Returns ``(argv_prefix, note)``. ``argv_prefix`` is a list to prepend to the
    scan arguments, or ``None`` when no runnable executable/wrapper exists (with
    ``note`` explaining why).
    """
    exe = os.path.expanduser(str(config.get('rmsx_executable', '') or '').strip())

    # Candidate native binaries: the configured path first, then the compiled
    # binary that ships in the RNAMotifScanX source tree (config points at the
    # unversioned bin/ folder, which may be empty on a fresh checkout).
    native_candidates = []
    if exe:
        native_candidates.append(exe)
        rmsx_root = Path(exe).resolve().parent.parent  # …/external/rmsx/bin/scan → …/external/rmsx
        native_candidates.append(str(rmsx_root / 'RNAMotifScanX_src' / 'scan'))
        native_candidates.append(str(rmsx_root / 'bin' / 'scan'))
    for candidate in native_candidates:
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK) \
                and _executable_runs_natively(candidate):
            return [candidate], f"native scan executable: {candidate}"

    wrapper = _find_scan_docker_wrapper(config)
    if wrapper and shutil.which('docker'):
        return [wrapper], f"docker scan wrapper: {wrapper}"

    reasons = []
    if not exe or not os.path.isfile(exe):
        reasons.append(f"configured rmsx_executable not found: {exe or '(unset)'}")
    elif not os.access(exe, os.X_OK):
        reasons.append(f"scan binary is not executable (chmod +x needed): {exe}")
    elif not _executable_runs_natively(exe):
        reasons.append(
            f"scan binary cannot run natively on {platform.system()}/{platform.machine()}: {exe}"
        )
    if wrapper and not shutil.which('docker'):
        reasons.append(f"docker wrapper present but Docker is not installed: {wrapper}")
    elif not wrapper:
        reasons.append("no scan_docker_*.sh wrapper found next to the executable")
    return None, '; '.join(reasons)


def _candidate_query_dirs(config: dict) -> list:
    """Directories that may hold ``{family}_consensus.struct`` templates."""
    dirs = []
    qd = os.path.expanduser(str(config.get('query_motifs_dir', '') or '').strip())
    if qd:
        dirs.append(qd)
    exe = os.path.expanduser(str(config.get('rmsx_executable', '') or '').strip())
    if exe:
        rmsx_root = Path(exe).resolve().parent.parent  # …/external/rmsx/bin/scan → …/external/rmsx
        dirs.append(str(rmsx_root / 'RNAMotifScanX_src' / 'Queries'))
        dirs.append(str(rmsx_root / 'RNAMotifScanX_src' / 'Queries' / 'reduced'))
        dirs.append(str(rmsx_root / 'queries'))
    seen, unique = set(), []
    for d in dirs:
        if d and d not in seen and os.path.isdir(d):
            seen.add(d)
            unique.append(d)
    return unique


def _resolve_family_query(family: str, query_dirs: list) -> str:
    names = [
        f"{family}_consensus.struct", f"{family}.struct",
        f"{family}_consensus.txt", f"{family}.txt",
    ]
    for directory in query_dirs:
        for name in names:
            candidate = os.path.join(directory, name)
            if os.path.isfile(candidate):
                return candidate
    return ''


def _rmsx_src_root(scan_prefix: list) -> str:
    """Locate the RNAMotifScanX source root that holds the ``mat/`` resources."""
    exe = scan_prefix[0] if scan_prefix else ''
    if not exe:
        return ''
    base = Path(exe).resolve().parent
    for candidate in (base, base.parent, base.parent / 'RNAMotifScanX_src'):
        if (candidate / 'mat').is_dir():
            return str(candidate)
    return ''


def _count_alignment_hits(log_path: str) -> int:
    """Count ``# Aligning`` alignment blocks in a scan output log."""
    count = 0
    try:
        with open(log_path, 'r', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                if re.search(r"^#\s+Aligning\s+", line):
                    count += 1
    except OSError:
        return 0
    return count


def _run_one_scan(scan_prefix, query_file, in_file, nch_file,
                  num_threads, run_dir, timeout=3600, cancel_event=None):
    """Run scan once for one (chain, family). Captures stdout/stderr/exit code.

    Supports cooperative cancellation: when ``cancel_event`` is set mid-run the
    scanner process is terminated and, for the docker wrapper, the container is
    force-removed so no work continues in the background.
    """
    os.makedirs(run_dir, exist_ok=True)
    map_pdb_abs = os.path.abspath(nch_file)
    # The docker wrapper only rewrites bare arguments that live under the repo
    # root into container paths, so pass --map_pdb as a separate token there.
    # For a native binary use the exact "=" form recorded in the validated logs.
    is_wrapper = bool(scan_prefix) and str(scan_prefix[0]).endswith('.sh')
    map_pdb_args = ['--map_pdb', map_pdb_abs] if is_wrapper else [f"--map_pdb={map_pdb_abs}"]
    cmd = list(scan_prefix) + [
        os.path.abspath(query_file),
        os.path.abspath(in_file),
        *map_pdb_args,
        '--pvalue', SCAN_PREPARED_PVALUE,
        '--num_threads', str(num_threads),
        '--write_alignment',
    ]
    stdout_path = os.path.join(run_dir, 'scan.stdout.log')
    stderr_path = os.path.join(run_dir, 'scan.stderr.log')
    cmd_path = os.path.join(run_dir, 'command.txt')
    with open(cmd_path, 'w', encoding='utf-8') as fh:
        fh.write(' '.join(cmd) + '\n')

    env = os.environ.copy()
    src_root = _rmsx_src_root(scan_prefix)
    if src_root:
        env['RNAMOTIFSCANX_PATH'] = src_root
    container_name = ''
    if is_wrapper:
        container_name = f"rmsx_scan_{os.getpid()}_{uuid.uuid4().hex[:8]}"
        env['RMSX_CONTAINER_NAME'] = container_name

    result = {
        'chain': '', 'family': '', 'cmd': cmd, 'cmd_str': ' '.join(cmd),
        'run_dir': run_dir, 'stdout': stdout_path, 'stderr': stderr_path,
        'container': container_name, 'exit_code': None, 'ok': False,
        'error': '', 'hits': 0,
    }

    def _terminate(proc):
        # Kill the scanner and, for docker, force-remove the container so the
        # scan actually stops instead of only skipping later tasks.
        try:
            proc.terminate()
        except Exception:
            pass
        if container_name:
            try:
                subprocess.run(['docker', 'rm', '-f', container_name],
                               capture_output=True, text=True, timeout=30)
            except Exception:
                pass
        try:
            proc.kill()
        except Exception:
            pass

    try:
        with open(stdout_path, 'w', encoding='utf-8') as out_fh, \
                open(stderr_path, 'w', encoding='utf-8') as err_fh:
            proc = subprocess.Popen(cmd, stdout=out_fh, stderr=err_fh, env=env,
                                    start_new_session=True)
            deadline = time.time() + timeout
            while True:
                try:
                    return_code = proc.wait(timeout=0.5)
                    break
                except subprocess.TimeoutExpired:
                    return_code = None
                if cancel_event is not None and cancel_event.is_set():
                    _terminate(proc)
                    result['error'] = 'cancelled'
                    return result
                if time.time() > deadline:
                    _terminate(proc)
                    result['error'] = 'timeout'
                    return result
    except OSError as exc:
        with open(stderr_path, 'a', encoding='utf-8') as fh:
            fh.write(f"OSError: {exc}\n")
        result['error'] = f"OSError: {exc}"
        return result

    result['exit_code'] = return_code
    result['ok'] = (return_code == 0)
    if not result['ok'] and not result['error']:
        result['error'] = 'nonzero_exit'
    if result['ok']:
        result['hits'] = _count_alignment_hits(stdout_path)
    return result


def run_scan_prepared(config: dict, pdb_id: str, output_dir: str, chains=None,
                      families=None, progress_cb=None, cancel_event=None) -> dict:
    """Run RNAMotifScanX against prepared inputs for ``pdb_id``.

    Returns a structured report. Never reads or substitutes preannotated data.
    """
    def emit(message: str):
        if progress_cb:
            try:
                progress_cb(message)
            except Exception:
                pass
        print(f"[scan_prepared] {message}")

    pdb_upper = str(pdb_id).strip().upper()
    prebuild_dir = os.path.expanduser(str(config.get('pdb_prebuild_dir', '') or ''))
    families = list(families or config.get('motif_families', list(DEFAULT_FAMILY_FOLDER_MAP.keys())))
    num_threads = int(config.get('num_threads', 4))
    output_dir = os.path.abspath(os.path.expanduser(output_dir))

    report = {
        'pdb_id': pdb_upper, 'output_dir': output_dir, 'executable': '',
        'runs': [], 'families': {}, 'failed_runs': [], 'problems': [],
        'ok': False, 'cancelled': False, 'total_hits': 0,
    }

    pairs, problems = locate_prepared_inputs(prebuild_dir, pdb_upper, chains)
    report['problems'].extend(problems)
    for problem in problems:
        emit(f"input problem: {problem}")
    if not pairs:
        emit("No readable prepared .rmsx.in/.rmsx.nch pairs found; aborting "
             "(preannotated data is NOT substituted).")
        return report

    emit("Using locally prepared RMSX inputs; MC-Annotate and RNAVIEW are skipped.")

    scan_prefix, note = resolve_scan_command(config)
    report['executable'] = note
    if not scan_prefix:
        emit(f"RMSX scan executable/wrapper is not runnable: {note}")
        report['problems'].append(note)
        return report
    emit(f"Scan command source: {note}")

    query_dirs = _candidate_query_dirs(config)
    if not query_dirs:
        emit("No query-template directory found; cannot supply motif templates.")
        report['problems'].append("no query template directory")
        return report

    os.makedirs(output_dir, exist_ok=True)

    for chain in sorted(pairs):
        if cancel_event is not None and cancel_event.is_set():
            report['cancelled'] = True
            emit(f"Cancelled before chain {chain}.")
            break
        paths = pairs[chain]
        emit(f"Running RNAMotifScanX for PDB {pdb_upper}, chain {chain}.")
        for family in families:
            if cancel_event is not None and cancel_event.is_set():
                report['cancelled'] = True
                emit("Cancelled during scan.")
                break
            query_file = _resolve_family_query(family, query_dirs)
            if not query_file:
                message = f"chain {chain}: no query template for family '{family}'"
                report['problems'].append(message)
                emit(message)
                continue
            run_dir = os.path.join(output_dir, '_runs', f"chain_{chain}", f"{family}_consensus")
            emit(f"  scanning family {family} (chain {chain})...")
            result = _run_one_scan(scan_prefix, query_file, paths['in'], paths['nch'],
                                   num_threads, run_dir, cancel_event=cancel_event)
            result['chain'] = chain
            result['family'] = family
            if result.get('error') == 'cancelled':
                report['cancelled'] = True
                report['runs'].append(result)
                emit(f"  chain {chain} / {family}: CANCELLED (scanner/container terminated)")
                break
            report['runs'].append(result)
            if result['ok']:
                family_dir = os.path.join(output_dir, f"{family}_consensus")
                os.makedirs(family_dir, exist_ok=True)
                aggregated = os.path.join(family_dir, 'result_0_100_withbs.log')
                with open(result['stdout'], 'r', encoding='utf-8', errors='ignore') as src, \
                        open(aggregated, 'a', encoding='utf-8') as dst:
                    dst.write(src.read())
                    dst.write('\n')
                family_entry = report['families'].setdefault(
                    family, {'log': aggregated, 'hits': 0, 'chains': []}
                )
                family_entry['hits'] += result['hits']
                family_entry['chains'].append(chain)
                emit(f"  chain {chain} / {family}: exit 0, {result['hits']} hit(s)")
            else:
                report['failed_runs'].append(result)
                emit(f"  chain {chain} / {family}: FAILED "
                     f"(exit {result.get('exit_code')}, {result.get('error')})")
        if report['cancelled']:
            break

    report['total_hits'] = sum(entry['hits'] for entry in report['families'].values())
    report['ok'] = (
        bool(report['runs'])
        and not report['failed_runs']
        and not report['cancelled']
    )
    return report


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Run RNAMotifScanX pipeline for a PDB structure'
    )
    parser.add_argument('--config', required=True,
                        help='Path to config/rmsx_config.json')
    parser.add_argument('--pdb', required=True,
                        help='PDB ID (e.g. 1S72)')
    parser.add_argument('--cif', default='',
                        help='Explicit path to CIF file (optional)')
    parser.add_argument('--check', action='store_true',
                        help='Only check which result files exist, do not run')
    parser.add_argument('--fresh', action='store_true',
                        help='Force fresh execution (ignore/remove cached outputs)')
    parser.add_argument('--scan-prepared', action='store_true',
                        help='Run scan against locally prepared .rmsx.in/.rmsx.nch inputs '
                             '(skips MC-Annotate/RNAVIEW; never uses preannotated data)')
    parser.add_argument('--chains', default='',
                        help='Comma-separated prepared chains to scan (default: all)')
    parser.add_argument('--out', default='',
                        help='Output directory for scan_prepared results')
    args = parser.parse_args()

    with open(args.config, 'r', encoding='utf-8') as fh:
        config = json.load(fh)

    # Resolve config-relative paths the same way RSMViewer does, so the
    # standalone CLI locates prepared inputs and the executable correctly.
    config_dir = os.path.dirname(os.path.abspath(args.config))
    for key in ('rmsx_executable', 'mc_annotate_executable', 'rnaview_executable',
                'rnaview_dir', 'pdb_prebuild_archive', 'pdb_prebuild_dir',
                'query_motifs_dir', 'cif_input_dir', 'output_dir'):
        value = str(config.get(key, '') or '').strip()
        if value:
            expanded = os.path.expanduser(value)
            if not os.path.isabs(expanded):
                expanded = os.path.abspath(os.path.join(config_dir, expanded))
            config[key] = expanded

    if args.scan_prepared:
        chains = [c.strip() for c in args.chains.split(',') if c.strip()]
        out_dir = args.out or os.path.join(
            str(config.get('output_dir', '.') or '.'), 'scan_prepared', args.pdb.upper()
        )
        report = run_scan_prepared(config, args.pdb, out_dir, chains=chains or None)
        print("\n" + "=" * 70)
        print(f"scan_prepared report for {report['pdb_id']}")
        print("=" * 70)
        print(f"Executable : {report['executable'] or '(none)'}")
        print(f"Output dir : {report['output_dir']}")
        print(f"Runs       : {len(report['runs'])}  "
              f"(failed: {len(report['failed_runs'])})")
        print(f"Total hits : {report['total_hits']}")
        print(f"Cancelled  : {report['cancelled']}")
        print(f"OK         : {report['ok']}")
        if report['problems']:
            print("Problems:")
            for problem in report['problems']:
                print(f"  - {problem}")
        for family, entry in sorted(report['families'].items()):
            print(f"  {family:<20} hits={entry['hits']:<4} chains={entry['chains']} -> {entry['log']}")
        return 0 if report['ok'] else 1

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
