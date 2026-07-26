#!/usr/bin/env python3
"""First-run setup helper for integrated Source-7 RMSX runtime."""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
from pathlib import Path


def detect_platform_dir() -> str:
    sys_name = platform.system().lower()
    machine = platform.machine().lower()
    if sys_name == 'darwin':
        return 'macos-arm64' if machine in ('arm64', 'aarch64') else 'macos-x86_64'
    if sys_name == 'windows':
        return 'windows-x86_64'
    return 'linux-x86_64'


def ensure_executable(path: Path) -> None:
    try:
        mode = path.stat().st_mode
        path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except Exception:
        pass


def find_first(paths: list[Path]) -> Path | None:
    for p in paths:
        if p.exists() and p.is_file():
            return p
    return None


def find_tool(name_candidates: list[str]) -> str:
    for name in name_candidates:
        loc = shutil.which(name)
        if loc:
            return loc
    return ''


def detect_boost_prefix() -> tuple[str, str]:
    env_prefixes = [
        os.environ.get('BOOST_ROOT', ''),
        os.environ.get('BOOST_HOME', ''),
    ]
    prefixes = [
        p for p in env_prefixes if p
    ] + [
        '/opt/homebrew',
        '/usr/local',
        '/usr',
        'C:/local/boost_1_84_0',
        'C:/vcpkg/installed/x64-windows',
    ]
    for prefix in prefixes:
        inc = Path(prefix) / 'include'
        lib = Path(prefix) / 'lib'
        if 'vcpkg' in prefix.lower():
            inc = Path(prefix) / 'include'
            lib = Path(prefix) / 'lib'
        if (inc / 'boost').exists() and lib.exists():
            return str(inc), str(lib)
    return '', ''


def detect_binary_format(path: Path) -> str:
    """Return a coarse binary format marker by magic bytes."""
    try:
        with path.open('rb') as fh:
            magic = fh.read(4)
    except Exception:
        return 'unknown'

    if magic.startswith(b'\x7fELF'):
        return 'elf'
    if magic[:2] == b'MZ':
        return 'pe'
    if magic in (b'\xfe\xed\xfa\xce', b'\xfe\xed\xfa\xcf', b'\xce\xfa\xed\xfe', b'\xcf\xfa\xed\xfe'):
        return 'macho'
    return 'unknown'


def is_binary_compatible(path: Path) -> bool:
    fmt = detect_binary_format(path)
    sys_name = platform.system().lower()
    if sys_name == 'windows':
        return fmt in ('pe', 'unknown')
    if sys_name == 'darwin':
        return fmt in ('macho', 'unknown')
    return fmt in ('elf', 'unknown')


def run_cmd(cmd: list[str], cwd: Path | None = None) -> tuple[bool, str]:
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            check=True,
            capture_output=True,
            text=True,
        )
        return True, (proc.stdout or '').strip()
    except subprocess.CalledProcessError as exc:
        msg = (exc.stderr or exc.stdout or '').strip()
        return False, msg[-2000:]


def _apply_mccore_compat_patches(repo_dir: Path) -> None:
    """Apply minimal source patches needed for modern libc++/compilers."""
    graph_h = repo_dir / 'lib' / 'Graph.h'
    if graph_h.exists():
        text = graph_h.read_text(encoding='utf-8', errors='ignore')
        if 'V& getHeadVertex () throw (NoSuchElementException)' in text:
            text = text.replace('V& getHeadVertex () throw (NoSuchElementException)', 'V& getHeadVertex (Graph& g) throw (NoSuchElementException)')
            text = text.replace('const V& getHeadVertex () const throw (NoSuchElementException)', 'const V& getHeadVertex (const Graph& g) const throw (NoSuchElementException)')
            text = text.replace('V& getTailVertex () throw (NoSuchElementException)', 'V& getTailVertex (Graph& g) throw (NoSuchElementException)')
            text = text.replace('const V& getTailVertex () const throw (NoSuchElementException)', 'const V& getTailVertex (const Graph& g) const throw (NoSuchElementException)')
            text = text.replace('if (vertices.size () <= head)', 'if (g.vertices.size () <= head)')
            text = text.replace('return vertices[head];', 'return g.vertices[head];')
            text = text.replace('if (vertices.size () <= tail)', 'if (g.vertices.size () <= tail)')
            text = text.replace('return vertices[tail];', 'return g.vertices[tail];')
        if 'protected:\n\n      /**\n       * Assigns the endvertices with the right\'s content.' in text:
            text = text.replace('protected:\n\n      /**\n       * Assigns the endvertices with the right\'s content.', 'public:\n\n      /**\n       * Assigns the endvertices with the right\'s content.')
        graph_h.write_text(text, encoding='utf-8')

    server_socket_cc = repo_dir / 'lib' / 'ServerSocket.cc'
    if server_socket_cc.exists():
        text = server_socket_cc.read_text(encoding='utf-8', errors='ignore')
        text = text.replace('if (bind (socket_id, (sockaddr*) &sin, (socklen_t) sizeof (sin)) < 0) {', 'if (::bind (socket_id, (sockaddr*) &sin, (socklen_t) sizeof (sin)) < 0) {')
        server_socket_cc.write_text(text, encoding='utf-8')


def _find_mccore_artifacts(install_prefix: Path) -> tuple[Path | None, Path | None]:
    include_root = None
    lib_path = None

    inc = install_prefix / 'include'
    if (inc / 'mccore').exists():
        include_root = inc
    elif any((inc).glob('mccore-*')):
        include_root = inc

    lib_candidates: list[Path] = []
    for lib_dir in [install_prefix / 'lib', install_prefix / 'lib64']:
        if lib_dir.exists():
            lib_candidates.extend(sorted(lib_dir.glob('libmccore.*')))
    if not lib_candidates:
        lib_candidates = sorted(install_prefix.rglob('libmccore.*'))
    if lib_candidates:
        lib_path = lib_candidates[0]

    return include_root, lib_path


def build_mccore_from_repo(runtime_dir: Path) -> tuple[bool, str, Path | None, Path | None, Path | None]:
    """Attempt to fetch/build/install MCCORE from official source."""
    git_exe = find_tool(['git'])
    cmake_exe = find_tool(['cmake'])
    if not git_exe:
        return False, 'git is required to fetch MCCORE source.', None, None, None
    if not cmake_exe:
        return False, 'cmake is required to build MCCORE source.', None, None, None

    ext_root = runtime_dir / 'src' / 'external'
    ext_root.mkdir(parents=True, exist_ok=True)
    repo_dir = ext_root / 'mccore'

    if not repo_dir.exists():
        ok, msg = run_cmd([git_exe, 'clone', '--depth', '1', 'https://github.com/major-lab/mccore.git', str(repo_dir)], cwd=ext_root)
        if not ok:
            return False, f'Failed to clone MCCORE: {msg}', None, None, None

    _apply_mccore_compat_patches(repo_dir)

    build_dir = ext_root / f'mccore-build-{detect_platform_dir()}'
    install_prefix = ext_root / f'mccore-install-{detect_platform_dir()}'
    build_dir.mkdir(parents=True, exist_ok=True)
    install_prefix.mkdir(parents=True, exist_ok=True)

    ok_cfg, cfg_msg = run_cmd([
        cmake_exe,
        '-S', str(repo_dir),
        '-B', str(build_dir),
        '-DCMAKE_BUILD_TYPE=Release',
        '-DCMAKE_POLICY_VERSION_MINIMUM=3.5',
        '-DWITH-RNAML=OFF',
        '-DWITH-MYSQL=OFF',
        f'-DCMAKE_INSTALL_PREFIX={install_prefix}',
    ], cwd=repo_dir)
    if not ok_cfg:
        return False, f'MCCORE configure failed: {cfg_msg}', None, None, None

    ok_build, build_msg = run_cmd([
        cmake_exe,
        '--build', str(build_dir),
        '--config', 'Release',
        '--target', 'mccore',
    ], cwd=repo_dir)
    if not ok_build:
        return False, f'MCCORE build failed: {build_msg}', None, None, None

    ok_install, install_msg = run_cmd([
        cmake_exe,
        '--install', str(build_dir),
        '--config', 'Release',
    ], cwd=repo_dir)
    if not ok_install:
        return False, f'MCCORE install failed: {install_msg}', None, None, None

    include_root, lib_path = _find_mccore_artifacts(install_prefix)
    if not include_root or not lib_path:
        return False, 'MCCORE built but include/library artifacts were not found after install.', None, None, None

    return True, 'MCCORE built from official source.', install_prefix, include_root, lib_path


def build_mc_annotate_from_repo(runtime_dir: Path, fetch_mccore: bool = False) -> tuple[bool, str, Path | None]:
    """Attempt to fetch and build MC-Annotate from official source."""
    git_exe = find_tool(['git'])
    cmake_exe = find_tool(['cmake'])
    if not git_exe:
        return False, 'git is required to fetch MC-Annotate source.', None
    if not cmake_exe:
        return False, 'cmake is required to build MC-Annotate source.', None

    ext_root = runtime_dir / 'src' / 'external'
    ext_root.mkdir(parents=True, exist_ok=True)
    repo_dir = ext_root / 'MC-Annotate'

    if not repo_dir.exists():
        ok, msg = run_cmd([git_exe, 'clone', '--depth', '1', 'https://github.com/major-lab/MC-Annotate.git', str(repo_dir)], cwd=ext_root)
        if not ok:
            return False, f'Failed to clone MC-Annotate: {msg}', None

    build_dir = ext_root / f'mc-annotate-build-{detect_platform_dir()}'
    build_dir.mkdir(parents=True, exist_ok=True)

    mccore_prefix = None
    mccore_include = None
    mccore_lib = None
    if fetch_mccore:
        ok_mccore, msg_mccore, mccore_prefix, mccore_include, mccore_lib = build_mccore_from_repo(runtime_dir)
        if not ok_mccore:
            return False, f'MC-Annotate prerequisite failed: {msg_mccore}', None

    cfg_cmd = [
        cmake_exe,
        '-S', str(repo_dir),
        '-B', str(build_dir),
        '-DCMAKE_BUILD_TYPE=Release',
        '-DCMAKE_POLICY_VERSION_MINIMUM=3.5',
    ]
    if mccore_prefix and mccore_include and mccore_lib:
        cfg_cmd.extend([
            f'-DCMAKE_PREFIX_PATH={mccore_prefix}',
            f'-DMCCORE_INCLUDE_DIR={mccore_include}',
            f'-DMCCORE_LIBRARY={mccore_lib}',
        ])

    ok_cfg, cfg_msg = run_cmd(cfg_cmd, cwd=repo_dir)
    if not ok_cfg:
        if 'MCCORE' in cfg_msg.upper() or 'MCCORE' in cfg_msg:
            return False, (
                'MC-Annotate configure failed: missing MCCORE dependency. '
                'Build/install MCCORE first, then retry setup. '
                f'Details: {cfg_msg}'
            ), None
        return False, f'MC-Annotate configure failed: {cfg_msg}', None

    ok_build, build_msg = run_cmd([
        cmake_exe,
        '--build', str(build_dir),
        '--config', 'Release',
        '--target', 'mcannotate',
    ], cwd=repo_dir)
    if not ok_build:
        return False, f'MC-Annotate build failed: {build_msg}', None

    exe_candidates = [
        build_dir / 'src' / 'mcannotate',
        build_dir / 'src' / 'Release' / 'mcannotate.exe',
        build_dir / 'Release' / 'mcannotate.exe',
        build_dir / 'mcannotate',
        build_dir / 'mcannotate.exe',
    ]
    exe = find_first(exe_candidates)
    if not exe:
        return False, 'MC-Annotate build succeeded but executable was not found.', None

    return True, 'MC-Annotate built from official source.', exe


def copy_if_exists(src_candidates: list[Path], dst: Path) -> bool:
    src = find_first(src_candidates)
    if not src:
        return False
    if not is_binary_compatible(src):
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    ensure_executable(dst)
    return True


def _patch_macos_mccore_runtime(mc_bin: Path, runtime_dir: Path, platform_dir: str) -> tuple[bool, str]:
    """Ensure MC-Annotate can load libmccore from the same bin folder on macOS."""
    ext_root = runtime_dir / 'src' / 'external'
    mccore_install = ext_root / f'mccore-install-{platform_dir}'
    lib_candidates = [
        mccore_install / 'lib64' / 'libmccore.2.0.0.dylib',
        mccore_install / 'lib' / 'libmccore.2.0.0.dylib',
    ]
    mccore_lib = find_first(lib_candidates)
    if not mccore_lib:
        return False, 'MC-Annotate built, but libmccore.2.0.0.dylib was not found for macOS runtime packaging.'

    dst = mc_bin.parent / 'libmccore.2.0.0.dylib'
    shutil.copy2(mccore_lib, dst)

    install_name_tool = find_tool(['install_name_tool'])
    if not install_name_tool:
        return False, 'install_name_tool not found; cannot patch MC-Annotate rpath on macOS.'

    # Add executable-relative rpath so @rpath/libmccore.2.0.0.dylib resolves.
    run_cmd([install_name_tool, '-add_rpath', '@executable_path', str(mc_bin)])
    return True, 'Bundled libmccore and patched MC-Annotate rpath.'


# RNAVIEW source files, in the order defined by its bundled Makefile.
RNAVIEW_SOURCE_FILES = [
    'rnaview.c', 'fpair.c', 'fpair_sub.c', 'pair_type.c', 'nrutil.c',
    'ps-xy.c', 'ps-xy-sub.c', 'vrml.c', 'rnaxml-new.c', 'analyze.c',
    'pattern.c', 'xml2ps.c', 'multiple.c', 'statistics.c',
]

# Compiler flags validated by the build-feasibility audit. RNAVIEW is late-1990s
# ANSI/K&R C: modern clang/gcc reject implicit declarations by default, and
# macOS' fortified libc aborts on a benign overlapping strcpy in rna(). These
# flags make it build and run on Linux and macOS without touching the source.
RNAVIEW_BUILD_FLAGS = [
    '-D_FORTIFY_SOURCE=0',
    '-Wno-implicit-function-declaration',
    '-Wno-implicit-int',
    '-Wno-return-type',
]


def find_rnaview_binary(runtime_dir: Path) -> Path | None:
    """Return a platform-compatible compiled rnaview binary in the runtime bin dir."""
    bin_root = runtime_dir / 'bin'
    if not bin_root.exists():
        return None
    names = ['rnaview.exe', 'rnaview'] if platform.system().lower() == 'windows' else ['rnaview', 'rnaview.exe']
    search_dirs = [bin_root / detect_platform_dir()]
    for p in sorted(bin_root.iterdir()):
        if p.is_dir() and p not in search_dirs:
            search_dirs.append(p)
    for bdir in search_dirs:
        for name in names:
            cand = bdir / name
            if cand.exists() and cand.is_file() and is_binary_compatible(cand):
                return cand
    return None


def build_rnaview_from_source(runtime_dir: Path, force: bool = False) -> tuple[bool, str]:
    """Compile the bundled ThirdParty/RNAVIEW source into the runtime bin dir.

    Supported by the build-feasibility audit on Linux and macOS, and attempted
    best-effort on Windows via MinGW/GCC. MSVC is intentionally not supported.
    Failure is non-fatal: callers fall back to MC-Annotate-only preprocessing.
    """
    sys_name = platform.system().lower()
    src_root = runtime_dir / 'src' / 'RNAMotifScanX_src' / 'ThirdParty' / 'RNAVIEW'
    src_dir = src_root / 'src'
    include_dir = src_root / 'include'
    if not src_dir.exists() or not include_dir.exists():
        return False, f'RNAVIEW source not found under {src_root}.'

    platform_dir = detect_platform_dir()
    bin_dir = runtime_dir / 'bin' / platform_dir
    out_name = 'rnaview.exe' if sys_name == 'windows' else 'rnaview'
    out_path = bin_dir / out_name

    # Detect an existing, platform-compatible build and skip unless forced.
    if not force:
        existing = find_rnaview_binary(runtime_dir)
        if existing is not None:
            return True, f'RNAVIEW already present ({existing}); skipping rebuild.'

    # Select a compiler. On Windows, only MinGW/GCC-style compilers are used;
    # MSVC (cl) is not supported and is deliberately not searched for.
    if sys_name == 'windows':
        cc = find_tool(['gcc', 'cc', 'x86_64-w64-mingw32-gcc', 'clang'])
        if not cc:
            return False, (
                'RNAVIEW build skipped: no MinGW/GCC compiler found on Windows. '
                'MSVC is not supported for RNAVIEW. RNAMotifScanX preprocessing will '
                'use the MC-Annotate-only fallback.'
            )
    else:
        cc = find_tool(['cc', 'clang', 'gcc'])
        if not cc:
            return False, 'RNAVIEW build skipped: no C compiler (cc/clang/gcc) found.'

    sources = [str(src_dir / name) for name in RNAVIEW_SOURCE_FILES]
    missing_sources = [s for s in sources if not Path(s).exists()]
    if missing_sources:
        return False, f'RNAVIEW build aborted: missing source files: {missing_sources}.'

    bin_dir.mkdir(parents=True, exist_ok=True)
    cmd = (
        [cc]
        + RNAVIEW_BUILD_FLAGS
        + ['-I', str(include_dir)]
        + sources
        + ['-o', str(out_path), '-lm']
    )
    ok, msg = run_cmd(cmd, cwd=src_root)
    if not ok:
        return False, f'RNAVIEW compilation failed ({Path(cc).name}): {msg[-1200:]}'

    if not out_path.exists():
        return False, 'RNAVIEW compilation reported success but the binary was not produced.'
    ensure_executable(out_path)
    if not is_binary_compatible(out_path):
        return False, f'RNAVIEW built but the binary is not compatible with this platform ({platform_dir}).'

    return True, f'RNAVIEW built from source into {out_path}.'


def validate_runtime(runtime_dir: Path, query_file: str = '') -> dict:
    platform_dir = detect_platform_dir()
    bin_root = runtime_dir / 'bin'
    bin_dir = bin_root / platform_dir
    queries_dir = runtime_dir / 'queries'
    mat_dir = runtime_dir / 'mat'
    query_path = Path(query_file).expanduser() if query_file else None

    candidate_bin_dirs = []
    if bin_dir.exists():
        candidate_bin_dirs.append(bin_dir)
    if bin_root.exists():
        for p in sorted(bin_root.iterdir()):
            if p.is_dir() and p not in candidate_bin_dirs:
                candidate_bin_dirs.append(p)

    rmsx = None
    mc = None
    resolved_bin_dir = bin_dir
    for bdir in candidate_bin_dirs:
        if not rmsx:
            rmsx = find_first([
                bdir / 'RNAMotifScanX', bdir / 'scan',
                bdir / 'RNAMotifScanX.exe', bdir / 'scan.exe',
            ])
            if rmsx:
                resolved_bin_dir = bdir
        if not mc:
            mc = find_first([
                bdir / 'MC-Annotate', bdir / 'mc-annotate',
                bdir / 'MC-Annotate.exe', bdir / 'mc-annotate.exe',
            ])
            if mc:
                resolved_bin_dir = bdir
        if rmsx and mc:
            break

    if query_path:
        missing_queries = [] if query_path.is_file() else [str(query_path)]
    else:
        expected_queries = [
            'k-turn_consensus.struct',
            'c-loop_consensus.struct',
            'sarcin-ricin_consensus.struct',
            'reverse-kturn_consensus.struct',
            'e-loop_consensus.struct',
        ]
        missing_queries = [q for q in expected_queries if not (queries_dir / q).exists()]

    missing = []
    if not rmsx:
        missing.append('RMSX executable')
    if not mc:
        missing.append('MC-Annotate executable')
    if missing_queries:
        missing.append(f'query files ({len(missing_queries)} missing)')
    if not mat_dir.exists() or not any(mat_dir.iterdir()):
        missing.append('matrix files')

    # RNAVIEW is optional: when absent, preprocessing falls back to MC-Annotate
    # only, so it never contributes to the 'missing'/'ok' gating.
    rnaview = find_rnaview_binary(runtime_dir)

    return {
        'ok': not missing,
        'platform_dir': platform_dir,
        'bin_dir': str(resolved_bin_dir),
        'rmsx_executable': str(rmsx) if rmsx else '',
        'mc_annotate_executable': str(mc) if mc else '',
        'rnaview_executable': str(rnaview) if rnaview else '',
        'missing': missing,
        'missing_queries': missing_queries,
    }


def build_from_source(runtime_dir: Path, fetch_mc_annotate: bool = False, fetch_mccore: bool = False) -> tuple[bool, str]:
    source_dir = runtime_dir / 'src' / 'RNAMotifScanX_src'
    platform_dir = detect_platform_dir()
    bin_dir = runtime_dir / 'bin' / platform_dir
    bin_dir.mkdir(parents=True, exist_ok=True)

    if not source_dir.exists():
        return False, f'Source tree not found: {source_dir}'

    sys_name = platform.system().lower()

    cmake_exe = find_tool(['cmake'])
    if not cmake_exe:
        return False, 'Build tool not found: CMake is required.'

    boost_inc, boost_lib = detect_boost_prefix()
    if not boost_inc or not boost_lib:
        mc_note = ''
        if fetch_mc_annotate:
            ok_mc, msg_mc, _ = build_mc_annotate_from_repo(runtime_dir, fetch_mccore=fetch_mccore)
            if ok_mc:
                mc_note = ' MC-Annotate source build is available.'
            else:
                mc_note = f' MC-Annotate check: {msg_mc}'
        if sys_name == 'darwin':
            return False, 'Boost headers/libs not found. Install Boost (for example: brew install boost) and retry.' + mc_note
        if sys_name == 'windows':
            return False, 'Boost headers/libs not found. Install Boost (for example via vcpkg/conan) and retry.' + mc_note
        return False, 'Boost headers/libs not found. Install Boost (for example: apt install libboost-all-dev) and retry.' + mc_note

    build_dir = runtime_dir / 'build' / detect_platform_dir()
    build_dir.mkdir(parents=True, exist_ok=True)

    configure_cmd = [
        cmake_exe,
        '-S', str(source_dir),
        '-B', str(build_dir),
        '-DCMAKE_BUILD_TYPE=Release',
        '-DBoost_NO_BOOST_CMAKE=ON',
        f'-DBoost_INCLUDE_DIR={boost_inc}',
        f'-DBoost_LIBRARY_DIRS={boost_lib}',
    ]

    boost_prefix = str(Path(boost_inc).parent)
    configure_cmd.extend([
        f'-DBOOST_ROOT={boost_prefix}',
        f'-DCMAKE_PREFIX_PATH={boost_prefix}',
    ])
    if sys_name == 'windows':
        configure_cmd.extend(['-A', 'x64'])

    try:
        subprocess.run(configure_cmd, cwd=str(source_dir), check=True, capture_output=True, text=True)
        subprocess.run(
            [cmake_exe, '--build', str(build_dir), '--config', 'Release', '--target', 'scan'],
            cwd=str(source_dir),
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        msg = (exc.stderr or exc.stdout or '').strip()
        return False, f'CMake build failed: {msg[-1500:]}'

    scan_dst = bin_dir / ('scan.exe' if sys_name == 'windows' else 'scan')
    copied_scan = copy_if_exists([
        build_dir / 'scan',
        build_dir / 'scan.exe',
        build_dir / 'Release' / 'scan.exe',
    ], scan_dst)

    mc_dst = bin_dir / ('MC-Annotate.exe' if sys_name == 'windows' else 'MC-Annotate')
    copied_mc = copy_if_exists([
        source_dir / 'StructureAnnotation' / 'MC-Annotate',
        source_dir / 'ThirdParty' / 'MC-Annotate',
    ], mc_dst)

    mc_note = ''
    if (not copied_mc) and fetch_mc_annotate:
        ok_mc, msg_mc, mc_exe = build_mc_annotate_from_repo(runtime_dir, fetch_mccore=fetch_mccore)
        mc_note = msg_mc
        if ok_mc and mc_exe:
            copied_mc = copy_if_exists(
                [mc_exe],
                mc_dst,
            )

    if copied_mc and sys_name == 'darwin':
        ok_patch, patch_msg = _patch_macos_mccore_runtime(mc_dst, runtime_dir, platform_dir)
        if not ok_patch:
            return False, f'MC-Annotate built but macOS runtime patch failed: {patch_msg}'

    if not copied_scan:
        return False, 'Build completed but scan executable was not found.'
    if not copied_mc:
        tail = f' {mc_note}' if mc_note else ''
        return False, (
            'scan built successfully, but no compatible MC-Annotate binary is available in the provided source bundle. '
            f'Current platform: {detect_platform_dir()}. '
            'The shipped MC-Annotate files are Linux-only ELF binaries.' + tail
        )

    return True, f'Installed runtime binaries into {bin_dir}'


def main() -> int:
    parser = argparse.ArgumentParser(description='RMSX integrated runtime setup')
    parser.add_argument('--runtime-dir', required=True, help='Path to rmsx_runtime directory')
    parser.add_argument('--build', action='store_true', help='Attempt build/install if runtime is incomplete')
    parser.add_argument('--fetch-mc-annotate', action='store_true', help='Allow cloning/building MC-Annotate from official repository')
    parser.add_argument('--fetch-mccore', action='store_true', help='Allow cloning/building MCCORE from official repository')
    parser.add_argument('--query-file', default='', help='Optional custom RNAMotifScanX query file to validate')
    parser.add_argument('--skip-rnaview', action='store_true', help='Do not attempt to build the bundled RNAVIEW during setup')
    parser.add_argument('--rebuild-rnaview', action='store_true', help='Force rebuilding RNAVIEW even if a compatible binary already exists')
    parser.add_argument('--json', action='store_true', help='Print JSON result')
    args = parser.parse_args()

    runtime_dir = Path(args.runtime_dir).resolve()
    result = validate_runtime(runtime_dir, query_file=args.query_file)
    setup_message = ''

    if args.build and not result['ok']:
        ok, setup_message = build_from_source(
            runtime_dir,
            fetch_mc_annotate=bool(args.fetch_mc_annotate),
            fetch_mccore=bool(args.fetch_mccore),
        )
        result = validate_runtime(runtime_dir, query_file=args.query_file)
        result['setup_attempted'] = True
        result['setup_ok'] = ok
        result['setup_message'] = setup_message
    else:
        result['setup_attempted'] = False

    # RNAVIEW is an optional preprocessing component. Build it during setup on
    # supported platforms (Linux/macOS, and Windows via MinGW). A failure here is
    # never fatal: preprocessing falls back to MC-Annotate only.
    if (args.build or args.rebuild_rnaview) and not args.skip_rnaview:
        rnaview_ok, rnaview_message = build_rnaview_from_source(
            runtime_dir, force=bool(args.rebuild_rnaview)
        )
        result['rnaview_build_attempted'] = True
        result['rnaview_build_ok'] = rnaview_ok
        result['rnaview_message'] = rnaview_message
        result['rnaview_executable'] = validate_runtime(
            runtime_dir, query_file=args.query_file
        ).get('rnaview_executable', '')
    else:
        result['rnaview_build_attempted'] = False

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"Platform dir: {result['platform_dir']}")
        print(f"Binary dir: {result['bin_dir']}")
        if result['ok']:
            print('RMSX runtime is ready')
        else:
            print('RMSX runtime is incomplete')
            for item in result.get('missing', []):
                print(f'  - missing: {item}')
            if setup_message:
                print(f'Setup: {setup_message}')
        if result.get('rnaview_executable'):
            print(f"RNAVIEW: {result['rnaview_executable']}")
        elif result.get('rnaview_build_attempted'):
            print('RNAVIEW: not available (MC-Annotate-only preprocessing will be used)')
        if result.get('rnaview_message'):
            print(f"RNAVIEW setup: {result['rnaview_message']}")

    return 0 if result['ok'] else 2


if __name__ == '__main__':
    sys.exit(main())
