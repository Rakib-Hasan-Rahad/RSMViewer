from __future__ import annotations

import os
import platform
from pathlib import Path
from typing import List, Optional


def validate_fr3d_checkout(repo_root: str) -> dict:
    """Validate a local FR3D checkout and report missing required paths."""
    root = Path(os.path.abspath(os.path.expanduser(repo_root))).resolve()
    if not root.exists() or not root.is_dir():
        return {'ok': False, 'repo_root': str(root), 'missing': ['checkout-does-not-exist']}

    required = [
        Path('fr3d/__init__.py'),
        Path('fr3d/search/FR3D.py'),
        Path('fr3d/search/query_processing.py'),
        Path('fr3d/classifiers/NA_pairwise_interactions.py'),
    ]
    missing = [str(item.as_posix()) for item in required if not (root / item).is_file()]
    # Newer FR3D checkouts may not ship fr3d_configuration.py in-tree. RSMViewer
    # can synthesize a runtime-compatible module when needed.
    optional_missing = []
    legacy_cfg = Path('fr3d/search/fr3d_configuration.py')
    if not (root / legacy_cfg).is_file():
        optional_missing.append(str(legacy_cfg.as_posix()))
    return {
        'ok': not missing,
        'repo_root': str(root),
        'missing': missing,
        'optional_missing': optional_missing,
    }


def discover_rmsx_runtime_dir(explicit_runtime_dir: str = '', fallback_dirs: Optional[List[str]] = None) -> str:
    """Discover an RMSX runtime directory from explicit config, env, or fallback locations."""
    env_value = os.environ.get('RSMVIEWER_RMSX_RUNTIME_DIR', '').strip()
    candidates = []
    if explicit_runtime_dir:
        candidates.append(explicit_runtime_dir)
    if env_value:
        candidates.append(env_value)
    if fallback_dirs:
        candidates.extend(fallback_dirs)

    for candidate in candidates:
        if not candidate:
            continue
        path = Path(os.path.abspath(os.path.expanduser(candidate))).resolve()
        if path.exists() and path.is_dir():
            return str(path)

    return ''


def get_runtime_platform_dir() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == 'darwin':
        return 'macos-arm64' if machine in ('arm64', 'aarch64') else 'macos-x86_64'
    if system == 'windows':
        return 'windows-x86_64'
    return 'linux-x86_64'
