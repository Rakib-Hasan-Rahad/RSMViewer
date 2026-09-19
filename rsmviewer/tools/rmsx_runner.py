#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RNAMotifScanX data handling and scan execution for RSMViewer.

Two ways RSMViewer gets RNAMotifScanX (RMSX) annotations, chosen by
``data_mode`` in ``config/rmsx_config.json``:

  preannotated     Precomputed ``*_consensus.log`` results. The requested PDB's
                   archive is downloaded from the project's results server the
                   first time (so results are always the current published
                   annotations) and extracted into ``pdb_prebuild_dir``.
  run_from_scratch Run the ``scan`` program on the prepared ``.rmsx.in`` /
                   ``.rmsx.nch`` inputs in ``pdb_prebuild_dir/<pdb>/<chain>/``.
                   The scanner runtime (native / WSL / Docker) comes from
                   :mod:`rmsx_runtime`. MC-Annotate and RNAVIEW are not run
                   here: the prepared inputs already contain their output.

Result folder layout read by the provider::

    <output>/<family>_consensus/result_0_100_withbs.log
"""

import os
import re
import shutil
import subprocess
import tarfile
import tempfile
import time
import urllib.error
from pathlib import Path

from .rmsx_runtime import (
    Runtime, build_command, ensure_runtime_bundle, find_query, find_src_root,
    open_url, query_dirs, resolve_runtime, setup,
)

DEFAULT_FAMILY_FOLDER_MAP = {
    'k-turn':        'k-turn_consensus',
    'c-loop':        'c-loop_consensus',
    'sarcin-ricin':  'sarcin-ricin_consensus',
    'reverse-kturn': 'reverse-kturn_consensus',
    'e-loop':        'e-loop_consensus',
}

DEFAULT_PREANNOTATED_BASE_URL = (
    "https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default"
)
DEFAULT_SCAN_TIMEOUT_SECONDS = 6 * 3600


# ─────────────────────────────────────────────────────────────────────────────
# Preannotated data
# ─────────────────────────────────────────────────────────────────────────────

def copy_preannotated_results(prebuild_dir: str, pdb_id: str, output_dir: str) -> dict:
    """Copy a PDB's preannotated consensus logs into the working output folder.

    Reads ``<prebuild_dir>/<pdb>/<chain>/<family>_consensus.log`` and writes one
    ``<output_dir>/<family>_consensus/result_0_100_withbs.log`` per family. A PDB
    may carry the same family for several chains (e.g. ``1s72/0`` and
    ``1s72/9``); all chains of a family are concatenated so no chain's hits are
    lost. Logs are always read fresh from ``prebuild_dir``, so a re-downloaded
    PDB is picked up immediately.

    Returns ``{"copied": <number of families written>, "output_dir": ...}``.
    """
    root = Path(os.path.expanduser(str(prebuild_dir or ""))).resolve()
    entry = root / str(pdb_id).lower()
    if not entry.is_dir():
        entry = root / str(pdb_id).upper()

    family_chunks: dict = {}
    if entry.is_dir():
        for log in sorted(entry.glob("*/*_consensus.log")):
            family_chunks.setdefault(log.stem, []).append(log.read_bytes())

    written = 0
    for family, chunks in family_chunks.items():
        destination = Path(output_dir) / family
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "result_0_100_withbs.log").write_bytes(b"\n".join(chunks))
        written += 1
    return {"copied": written, "output_dir": str(output_dir)}

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


def download_preannotated_pdb(base_url: str, pdb_id: str, dest_dir: str,
                              timeout: int = 60) -> dict:
    """Download ``<base_url>/<pdb_lower>.tar.gz`` and extract it under ``dest_dir``.

    The archive is expected to contain a top-level ``<pdb_lower>/`` folder, i.e.
    the same layout as ``rmsx_work_default/<pdb_lower>/``; if it does not, its
    contents are placed under ``dest_dir/<pdb_lower>/`` instead. Extraction is
    staged in a temporary directory and moved into place only when complete, so
    an interrupted download never leaves a partial ``<pdb_lower>/`` folder that
    would later be mistaken for real local data.

    Returns ``{"ok": bool, "url": str, "path": str, "error": str,
    "insecure_tls": bool}``; never raises.
    """
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

        response, insecure = open_url(url, timeout=timeout)
        with response, open(archive_path, "wb") as out_fh:
            shutil.copyfileobj(response, out_fh)
        result["insecure_tls"] = insecure

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
# Prepared inputs
# ─────────────────────────────────────────────────────────────────────────────

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


def ensure_prepared_inputs(config: dict, pdb_id: str, chains=None, say=None):
    """Prepared ``.rmsx.in``/``.rmsx.nch`` pairs for ``pdb_id``, fetching them if absent.

    Looks in ``pdb_prebuild_dir/<pdb>/<chain>/``. When nothing is there, the
    PDB's archive is downloaded from the results server (it carries the
    prepared inputs alongside the precomputed logs); the precomputed logs are
    never used as scan results. Returns ``(pairs, problems)``.
    """
    say = say or (lambda message: None)
    prebuild_dir = os.path.expanduser(str(config.get('pdb_prebuild_dir', '') or ''))
    pairs, problems = locate_prepared_inputs(prebuild_dir, pdb_id, chains)
    if pairs or not prebuild_dir:
        return pairs, problems
    base = str(config.get('preannotated_base_url') or DEFAULT_PREANNOTATED_BASE_URL)
    say(f"No prepared inputs for {str(pdb_id).upper()} locally; downloading them from {base} ...")
    downloaded = download_preannotated_pdb(base, pdb_id, prebuild_dir)
    if not downloaded['ok']:
        return {}, problems + [f"could not download prepared inputs: {downloaded['error']}"]
    return locate_prepared_inputs(prebuild_dir, pdb_id, chains)


# ─────────────────────────────────────────────────────────────────────────────
# Running the scanner
# ─────────────────────────────────────────────────────────────────────────────

_ALIGNMENT_END = re.compile(r">{5,}\s+ALIGNMENT ENDS\s+>{5,}")


def _complete_alignments(text: str) -> str:
    """The output up to and including the last complete alignment block, or ''.

    ``scan`` can segfault part-way through a large structure (the published
    preannotated logs contain the same ``# scan failed rc=-11``). Everything it
    printed before that is valid; the truncated final block is not.
    """
    matches = list(_ALIGNMENT_END.finditer(text))
    return text[:matches[-1].end()] + "\n" if matches else ""


def _run_one_scan(runtime: Runtime, src_root: str, query_file: str, in_file: str,
                  nch_file: str, num_threads: int, run_dir: str, timeout: float) -> dict:
    """Run ``scan`` once for one (chain, family) and capture stdout/stderr/exit code."""
    os.makedirs(run_dir, exist_ok=True)
    argv, env, container = build_command(
        runtime, src_root, os.path.abspath(query_file), os.path.abspath(in_file),
        os.path.abspath(nch_file), num_threads,
    )
    stdout_path = os.path.join(run_dir, 'scan.stdout.log')
    stderr_path = os.path.join(run_dir, 'scan.stderr.log')
    with open(os.path.join(run_dir, 'command.txt'), 'w', encoding='utf-8') as fh:
        fh.write(' '.join(argv) + '\n')

    result = {
        'cmd': argv, 'run_dir': run_dir, 'stdout': stdout_path, 'stderr': stderr_path,
        'exit_code': None, 'ok': False, 'error': '', 'hits': 0, 'seconds': 0.0,
    }
    started = time.time()
    try:
        with open(stdout_path, 'w', encoding='utf-8') as out_fh, \
                open(stderr_path, 'w', encoding='utf-8') as err_fh:
            proc = subprocess.Popen(argv, stdout=out_fh, stderr=err_fh, env=env)
            try:
                result['exit_code'] = proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
                if container:
                    subprocess.run(['docker', 'rm', '-f', container], capture_output=True, timeout=60)
                result['error'] = f'timeout after {int(timeout)}s'
                return result
            finally:
                result['seconds'] = round(time.time() - started, 1)
    except OSError as exc:
        result['error'] = f'{type(exc).__name__}: {exc}'
        return result

    result['log'] = stdout_path
    if result['exit_code'] == 0:
        result['ok'] = True
        result['hits'] = _count_alignment_hits(stdout_path)
        return result

    result['error'] = 'nonzero_exit'
    try:
        with open(stdout_path, 'r', encoding='utf-8', errors='ignore') as fh:
            kept = _complete_alignments(fh.read())
    except OSError:
        kept = ''
    if kept:
        # Crashed after producing complete alignments: keep those, as the published logs do.
        complete_path = os.path.join(run_dir, 'scan.complete.log')
        with open(complete_path, 'w', encoding='utf-8') as fh:
            fh.write(kept)
            fh.write(f"\n# scan failed rc={result['exit_code']} (output kept up to the last complete alignment)\n")
        result.update(ok=True, partial=True, log=complete_path, hits=_count_alignment_hits(complete_path),
                      error=f"scanner crashed (exit {result['exit_code']}) after {_count_alignment_hits(complete_path)} complete hit(s)")
    return result


def run_scan_prepared(config: dict, pdb_id: str, output_dir: str, chains=None,
                      families=None, progress_cb=None, runtime: Runtime = None) -> dict:
    """Run RNAMotifScanX on the prepared inputs for ``pdb_id`` and write result logs.

    Returns a report dict (``ok``, ``runs``, ``failed_runs``, ``families``,
    ``problems``, ``total_hits``, ``runtime``). Never reads or substitutes
    preannotated results, and never raises for expected failures.
    """
    def emit(message: str):
        if progress_cb:
            try:
                progress_cb(message)
            except Exception:
                pass

    pdb_upper = str(pdb_id).strip().upper()
    runtime_dir = str(config.get('rmsx_runtime_dir') or '')
    families = list(families or config.get('motif_families') or DEFAULT_FAMILY_FOLDER_MAP)
    num_threads = int(config.get('num_threads', 4))
    timeout = float(config.get('scan_timeout_seconds') or DEFAULT_SCAN_TIMEOUT_SECONDS)
    output_dir = os.path.abspath(os.path.expanduser(output_dir))

    report = {
        'pdb_id': pdb_upper, 'output_dir': output_dir, 'runtime': '',
        'runs': [], 'families': {}, 'failed_runs': [], 'partial_runs': [], 'problems': [],
        'ok': False, 'total_hits': 0,
    }

    pairs, problems = ensure_prepared_inputs(config, pdb_upper, chains, emit)
    report['problems'].extend(problems)
    if not pairs:
        report['problems'].append(
            f"no readable prepared .rmsx.in/.rmsx.nch pairs for {pdb_upper} "
            "(preannotated results are not substituted)")
        return report

    # First use: fetch the scanner runtime (source, scoring matrices, query
    # models) from the project's server. Already there -> nothing is downloaded.
    bundle = ensure_runtime_bundle(config, runtime_dir, emit)
    if not bundle['ok']:
        report['problems'].append(
            f"could not download the RNAMotifScanX runtime from {bundle['url']}: {bundle['error']}. "
            "Download it manually and extract it into external/rmsx/ (see external/rmsx_setup.md)")
        return report

    if runtime is None:
        runtime, reasons = resolve_runtime(config, runtime_dir)
        if runtime is None:
            # No scanner runs on this machine yet: do the one-time preparation
            # (build from source / WSL2 / Docker) now instead of stopping.
            emit("No RNAMotifScanX scanner is ready yet; preparing one (first run only)...")
            prepared = setup(config, runtime_dir, emit, install_deps=True)
            runtime = prepared['runtime'] if prepared['ok'] else None
            if runtime is None:
                report['problems'].append(
                    "no RNAMotifScanX scanner is available on this machine: "
                    + "; ".join(reasons + prepared['problems']))
                report['problems'].extend(prepared['next'])
                return report
    report['runtime'] = runtime.note
    emit(f"Scanner: {runtime.note}")

    src_root = find_src_root(runtime_dir)
    if src_root is None:
        report['problems'].append(f"RNAMotifScanX scoring matrices (mat/) not found under {runtime_dir}")
        return report
    qdirs = query_dirs(config, runtime_dir)
    if not qdirs:
        report['problems'].append("no query-model directory found")
        return report

    os.makedirs(output_dir, exist_ok=True)
    total = len(pairs) * len(families)
    done = 0
    for chain in sorted(pairs):
        paths = pairs[chain]
        for family in families:
            done += 1
            query_file = find_query(family, qdirs)
            if not query_file:
                message = f"chain {chain}: no query model for family '{family}'"
                report['problems'].append(message)
                emit(message)
                continue
            emit(f"[{done}/{total}] scanning {pdb_upper} chain {chain} for {family} ...")
            run_dir = os.path.join(output_dir, '_runs', f"chain_{chain}", f"{family}_consensus")
            result = _run_one_scan(runtime, str(src_root), query_file, paths['in'], paths['nch'],
                                   num_threads, run_dir, timeout)
            result['chain'], result['family'] = chain, family
            report['runs'].append(result)
            if not result['ok']:
                report['failed_runs'].append(result)
                detail = _tail(result['stderr'])
                emit(f"    FAILED ({result['error']}, exit {result['exit_code']}) {detail}")
                continue
            family_dir = os.path.join(output_dir, f"{family}_consensus")
            os.makedirs(family_dir, exist_ok=True)
            aggregated = os.path.join(family_dir, 'result_0_100_withbs.log')
            if result.get('partial'):
                report['partial_runs'].append(result)
                emit(f"    WARNING: {result['error']}; later hits are missing")
            with open(result['log'], 'r', encoding='utf-8', errors='ignore') as src, \
                    open(aggregated, 'a', encoding='utf-8') as dst:
                dst.write(src.read())
                dst.write('\n')
            entry = report['families'].setdefault(family, {'log': aggregated, 'hits': 0, 'chains': []})
            entry['hits'] += result['hits']
            entry['chains'].append(chain)
            emit(f"    done in {result['seconds']}s: {result['hits']} hit(s)")

    report['total_hits'] = sum(e['hits'] for e in report['families'].values())
    report['ok'] = bool(report['runs']) and not report['failed_runs']
    return report


def _tail(path: str, n: int = 3) -> str:
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
            lines = [line.strip() for line in fh.read().strip().splitlines() if line.strip()]
        return ' | '.join(lines[-n:])[:300]
    except OSError:
        return ''
