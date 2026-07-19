#!/usr/bin/env python3
"""First-run setup helper for integrated Source-5 FR3D runtime."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def _which_python_candidates() -> list[str]:
    candidates = []
    env_python = os.environ.get('FR3D_PYTHON', '').strip()
    if env_python:
        candidates.append(env_python)

    for value in [sys.executable, '/opt/homebrew/bin/python3', '/usr/local/bin/python3', 'python3']:
        if value and value not in candidates:
            candidates.append(value)
    return candidates


def _looks_like_fr3d_root(path: Path) -> bool:
    return (path / 'fr3d' / '__init__.py').exists()


def _detect_fr3d_root(runtime_dir: Path) -> str:
    env_root = os.environ.get('FR3D_ROOT', '').strip()
    if env_root:
        p = Path(os.path.expanduser(env_root)).resolve()
        if _looks_like_fr3d_root(p):
            return str(p)

    candidates = [
        runtime_dir / 'python' / 'fr3d-python',
        runtime_dir / 'python',
        runtime_dir.parent.parent.parent / 'fr3d-python',
        runtime_dir.parent.parent / 'fr3d-python',
        runtime_dir.parent / 'fr3d-python',
    ]

    for candidate in candidates:
        p = candidate.resolve()
        if _looks_like_fr3d_root(p):
            return str(p)

    return ''


def _python_has_fr3d_deps(python_exe: str) -> bool:
    try:
        proc = subprocess.run(
            [python_exe, '-c', 'import scipy, pdbx'],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
        return proc.returncode == 0
    except Exception:
        return False


def _resolve_executable(candidate: str) -> str:
    if not candidate:
        return ''
    expanded = os.path.abspath(os.path.expanduser(candidate))
    if os.path.isfile(expanded):
        return expanded
    found = shutil.which(candidate)
    return os.path.abspath(found) if found else ''


def _detect_python(preferred: str = '') -> str:
    candidates = []
    if preferred:
        candidates.append(preferred)
    candidates.extend(_which_python_candidates())

    seen = set()
    for candidate in candidates:
        expanded = _resolve_executable(candidate)
        if not expanded:
            continue
        if expanded in seen:
            continue
        seen.add(expanded)
        if _python_has_fr3d_deps(expanded):
            return expanded
    return ''


def _run_cmd(cmd: list[str], cwd: Path | None = None, timeout: int = 900) -> tuple[bool, str]:
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
        msg = ((proc.stdout or '') + '\n' + (proc.stderr or '')).strip()
        return proc.returncode == 0, msg[-4000:]
    except Exception as exc:
        return False, f'{type(exc).__name__}: {exc}'


def _bootstrap_runtime(runtime_dir: Path) -> list[str]:
    """Attempt to provision local FR3D runtime assets for strict Source-5 mode."""
    logs: list[str] = []
    python_dir = runtime_dir / 'python'
    venv_dir = python_dir / 'venv'
    fr3d_repo = python_dir / 'fr3d-python'

    host_python = _resolve_executable(os.environ.get('FR3D_BOOTSTRAP_PYTHON', '').strip())
    if not host_python:
        host_python = _resolve_executable(sys.executable) or _resolve_executable('python3')

    if not host_python:
        logs.append('Could not find host python to create runtime venv.')
        return logs

    python_dir.mkdir(parents=True, exist_ok=True)

    # Create venv once.
    venv_python = venv_dir / 'bin' / 'python'
    if not venv_python.exists():
        ok, msg = _run_cmd([host_python, '-m', 'venv', str(venv_dir)])
        logs.append('Create venv: ok' if ok else f'Create venv failed: {msg}')
        if not ok:
            return logs

    venv_python_str = str(venv_python)

    # Install required dependencies for runtime detection + FR3D execution.
    ok, msg = _run_cmd([venv_python_str, '-m', 'pip', 'install', '--upgrade', 'pip'])
    logs.append('Upgrade pip: ok' if ok else f'Upgrade pip failed: {msg}')

    ok, msg = _run_cmd([venv_python_str, '-m', 'pip', 'install', 'numpy', 'scipy', 'mmcif-pdbx'])
    logs.append('Install deps (numpy/scipy/mmcif-pdbx): ok' if ok else f'Install deps failed: {msg}')

    # Ensure FR3D source checkout exists under runtime/python/fr3d-python.
    if not _looks_like_fr3d_root(fr3d_repo):
        git_exe = shutil.which('git')
        if not git_exe:
            logs.append('Git not found; cannot auto-clone FR3D source.')
            return logs
        if fr3d_repo.exists() and not _looks_like_fr3d_root(fr3d_repo):
            try:
                shutil.rmtree(fr3d_repo)
            except Exception as exc:
                logs.append(f'Could not clean invalid FR3D checkout: {exc}')
                return logs
        ok, msg = _run_cmd([git_exe, 'clone', '--depth', '1', 'https://github.com/BGSU-RNA/fr3d-python.git', str(fr3d_repo)])
        logs.append('Clone FR3D source: ok' if ok else f'Clone FR3D source failed: {msg}')
        if not ok:
            return logs

    # Install FR3D package into venv to maximize import compatibility.
    ok, msg = _run_cmd([venv_python_str, '-m', 'pip', 'install', '-e', str(fr3d_repo)])
    logs.append('Install FR3D package: ok' if ok else f'Install FR3D package failed: {msg}')

    return logs


def _build_report(runtime_dir: Path, do_setup: bool = False) -> dict:
    runtime_dir = runtime_dir.resolve()
    user_ann_root = runtime_dir.parent.parent / 'database' / 'user_annotations' / 'fr3d'
    cache_dir = runtime_dir / 'cache'
    scripts_dir = runtime_dir / 'scripts'
    python_dir = runtime_dir / 'python'
    bin_dir = runtime_dir / 'bin'

    if do_setup:
        for p in [cache_dir, scripts_dir, python_dir, bin_dir, user_ann_root]:
            p.mkdir(parents=True, exist_ok=True)
        bootstrap_logs = _bootstrap_runtime(runtime_dir)
    else:
        bootstrap_logs = []

    fr3d_root = _detect_fr3d_root(runtime_dir)
    preferred_python = str((python_dir / 'venv' / 'bin' / 'python').resolve()) if (python_dir / 'venv' / 'bin' / 'python').exists() else ''
    python_exe = _detect_python(preferred=preferred_python)

    missing = []
    if not python_exe:
        missing.append('python_executable_with_scipy_pdbx')
    if not fr3d_root:
        missing.append('fr3d_python_dir')

    setup_message = ''
    if not fr3d_root:
        setup_message += (
            'FR3D local runtime not found. Source 5 requires local FR3D runtime. '
            'Set FR3D_ROOT to enable local FR3D extraction. '
        )
    if not python_exe:
        setup_message += (
            'No Python interpreter with scipy+pdbx was found. '
            'Install dependencies to enable local FR3D extraction. '
        )
    if bootstrap_logs:
        setup_message += ' '.join(bootstrap_logs)

    report = {
        'ok': bool(fr3d_root and python_exe),
        'runtime_mode': 'local_pipeline' if (fr3d_root and python_exe) else 'not_ready',
        'runtime_dir': str(runtime_dir),
        'fr3d_runtime_dir': str(runtime_dir),
        'fr3d_python_dir': fr3d_root,
        'python_executable': python_exe,
        'output_dir': str(user_ann_root.resolve()),
        'cache_dir': str(cache_dir.resolve()),
        'scripts_dir': str(scripts_dir.resolve()),
        'python_dir': str(python_dir.resolve()),
        'bin_dir': str(bin_dir.resolve()),
        'missing': missing,
        'setup_message': setup_message.strip(),
        'bootstrap_logs': bootstrap_logs,
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description='FR3D runtime setup/doctor helper')
    parser.add_argument('--runtime-dir', default='', help='Path to fr3d_runtime directory')
    parser.add_argument('--build', action='store_true', help='Create runtime cache/scripts/bin directories')
    parser.add_argument('--json', action='store_true', help='Emit JSON only')
    args = parser.parse_args()

    runtime_dir = Path(args.runtime_dir).resolve() if args.runtime_dir else Path(__file__).resolve().parent
    report = _build_report(runtime_dir, do_setup=args.build)

    if args.json:
        print(json.dumps(report))
    else:
        print(json.dumps(report, indent=2))

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
