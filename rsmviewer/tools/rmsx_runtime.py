"""Find, build, verify and invoke the RNAMotifScanX ``scan`` executable.

RNAMotifScanX is a C++ program. The copy shipped in this repository is a
Linux x86-64 binary, so every other platform needs another way to run it. This
module provides three "runtimes" and picks the first one that works:

  native  a ``scan`` that runs directly on this OS/CPU: the bundled ELF on
          Linux x86-64, or a copy built from the bundled source by
          :func:`setup` on macOS / other Linux (fast, no virtualisation);
  wsl     Windows only: the Linux binary run through WSL2;
  docker  any OS: the Linux binary inside an x86-64 container.

:func:`setup` prepares one (``rmv_setup RNAMotifScanX``), :func:`diagnose`
explains what is and is not available (``rmv_rmsx_doctor``), and
:func:`build_command` turns a scan request into the exact command line for the
chosen runtime, translating file paths where the runtime needs it.
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import os
import platform
import re
import shutil
import ssl
import subprocess
import tarfile
import tempfile
import threading
import time
import urllib.error
import urllib.request
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from .runtime_layout import get_runtime_platform_dir

DOCKER_IMAGE = "ubuntu:22.04"
DEFAULT_BUNDLE_URL = "https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz"
SCAN_SOURCES = [
    "structural_motif", "annotated_structure", "scan_structure",
    "motif_graph_matching", "simulation", "parameters", "map_PDB_index",
    "main_scan",
]
BOOST_LIBS = ["iostreams", "filesystem", "thread", "program_options"]
SCAN_PVALUE = "1.0"

Progress = Callable[[str], None]


# ── small helpers ───────────────────────────────────────────────────────────

def _decode(data: bytes) -> str:
    if not data:
        return ""
    if b"\x00" in data:  # wsl.exe prints UTF-16
        return data.decode("utf-16-le", errors="ignore")
    return data.decode("utf-8", errors="replace")


def run_capture(argv: List[str], timeout: float = 60, env: Optional[dict] = None,
                cwd: Optional[str] = None) -> Tuple[int, str, str]:
    """Run a command, returning (returncode, stdout, stderr); never raises."""
    try:
        proc = subprocess.run(argv, capture_output=True, timeout=timeout, env=env, cwd=cwd)
    except subprocess.TimeoutExpired:
        return 124, "", f"timed out after {timeout}s"
    except OSError as exc:
        return 127, "", f"{type(exc).__name__}: {exc}"
    return proc.returncode, _decode(proc.stdout), _decode(proc.stderr)


def _last_line(text: str, default: str = "") -> str:
    lines = [line.strip() for line in (text or "").strip().splitlines() if line.strip()]
    return lines[-1][:240] if lines else default


_HOST_MACHINE: Optional[str] = None


def host_machine() -> str:
    """CPU architecture of the *hardware*, lower-case (``arm64``, ``x86_64``...).

    ``platform.machine()`` reports the architecture of the running process, so an
    Intel build of PyMOL running under Rosetta on an Apple-silicon Mac says
    ``x86_64`` although the machine is ``arm64``. Homebrew's libraries and the
    scanner we build are arm64 there, so the hardware is what matters.
    """
    global _HOST_MACHINE
    if _HOST_MACHINE is None:
        machine = platform.machine().lower()
        if platform.system() == "Darwin":
            rc, out, _ = run_capture(["/usr/sbin/sysctl", "-n", "hw.optional.arm64"], timeout=10)
            if rc == 0 and out.strip() == "1":
                machine = "arm64"
        _HOST_MACHINE = "arm64" if machine == "aarch64" else machine
    return _HOST_MACHINE


def is_translated() -> bool:
    """True when this process is an Intel build running under Rosetta on an ARM Mac."""
    return platform.system() == "Darwin" and host_machine() == "arm64" \
        and platform.machine().lower() != "arm64"


def platform_dir() -> str:
    if platform.system() == "Darwin":
        return "macos-arm64" if host_machine() == "arm64" else "macos-x86_64"
    return get_runtime_platform_dir()


def exe_name() -> str:
    return "scan.exe" if platform.system() == "Windows" else "scan"


def open_url(url: str, timeout: float = 60, method: str = "GET"):
    """Open a URL. Returns ``(response, insecure_tls)``.

    Verified TLS first. Some Python installs (notably python.org builds on
    macOS) ship no CA bundle and reject valid certificates, so retry with
    certifi's bundle, and only as a last resort without verification.
    """
    request = urllib.request.Request(url, headers={"User-Agent": "RSMViewer"}, method=method)
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
            return urllib.request.urlopen(request, timeout=timeout, context=context), mode == "unverified"
        except urllib.error.URLError as exc:
            if isinstance(exc.reason, ssl.SSLCertVerificationError) and index < len(contexts) - 1:
                continue
            raise
    raise RuntimeError("unreachable")


# ── layout ──────────────────────────────────────────────────────────────────

def layout(runtime_dir: str) -> Dict[str, Path]:
    root = Path(os.path.expanduser(str(runtime_dir))).resolve()
    return {"root": root, "src": root / "RNAMotifScanX_src", "bin": root / "bin" / platform_dir()}


def find_src_root(runtime_dir: str) -> Optional[Path]:
    """The RNAMotifScanX source tree; it holds the ``mat/`` scoring matrices."""
    src = layout(runtime_dir)["src"]
    return src if (src / "mat").is_dir() else None


def query_dirs(config: dict, runtime_dir: str) -> List[Path]:
    """Directories searched for ``<family>_consensus.struct`` query models.

    ``Queries/reduced`` comes before ``Queries``: the reduced models are the
    ones that reproduce the published preannotated results (verified against
    the shipped 1S72 logs); the full set finds additional marginal hits.
    """
    dirs: List[Path] = []
    configured = os.path.expanduser(str(config.get("query_motifs_dir") or "").strip())
    if configured:
        dirs.append(Path(configured))
    src = layout(runtime_dir)["src"]
    dirs += [src / "Queries" / "reduced", src / "Queries"]
    seen, out = set(), []
    for directory in dirs:
        if directory.is_dir() and directory not in seen:
            seen.add(directory)
            out.append(directory)
    return out


def find_query(family: str, dirs: List[Path]) -> str:
    for directory in dirs:
        for name in (f"{family}_consensus.struct", f"{family}.struct",
                     f"{family}_consensus.txt", f"{family}.txt"):
            candidate = directory / name
            if candidate.is_file():
                return str(candidate)
    return ""


# ── executables ─────────────────────────────────────────────────────────────

def _magic(path: Path) -> bytes:
    try:
        with open(path, "rb") as fh:
            return fh.read(20)
    except OSError:
        return b""


def is_linux_x86_elf(path: Path) -> bool:
    head = _magic(path)
    return head[:4] == b"\x7fELF" and len(head) >= 20 and head[4] == 2 and head[18:20] == b"\x3e\x00"


def binary_matches_host(path: Path) -> bool:
    """True when ``path`` is an executable format this OS/CPU can run directly."""
    head = _magic(path)
    system, machine = platform.system(), platform.machine().lower()
    if head[:2] == b"#!":
        return system != "Windows"
    if head[:4] == b"\x7fELF":
        if system != "Linux":
            return False
        elf_machine = int.from_bytes(head[18:20], "little") if len(head) >= 20 else 0
        return {0x3E: machine in ("x86_64", "amd64"), 0xB7: machine in ("aarch64", "arm64")}.get(elf_machine, False)
    if head[:2] == b"MZ":
        return system == "Windows"
    if head[:4] == b"\xcf\xfa\xed\xfe":  # thin 64-bit Mach-O: check its CPU type
        if system != "Darwin":
            return False
        cputype = int.from_bytes(head[4:8], "little")
        allowed = {0x01000007}                      # x86_64 (Rosetta runs it on ARM too)
        if host_machine() == "arm64":
            allowed.add(0x0100000C)                 # arm64
        return cputype in allowed
    if head[:4] in (b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca"):
        return system == "Darwin"
    return False


def probe(argv_prefix: List[str], env: Optional[dict] = None, timeout: float = 90) -> Tuple[bool, str]:
    """Start the scanner with ``--help`` to prove it launches and loads its libraries."""
    rc, out, err = run_capture(list(argv_prefix) + ["--help"], timeout=timeout, env=env)
    text = out + err
    if "Allowed options" in text:
        return True, "starts and responds to --help"
    return False, f"{_last_line(text, 'no output')} (exit {rc})"


def linux_elf(config: dict, runtime_dir: str) -> Optional[Path]:
    """A Linux x86-64 ``scan`` binary, for the WSL and Docker runtimes."""
    lay = layout(runtime_dir)
    configured = os.path.expanduser(str(config.get("rmsx_executable") or "").strip())
    for candidate in ([Path(configured)] if configured else []) + [
        lay["src"] / "scan", lay["root"] / "bin" / "linux-x86_64" / "scan"
    ]:
        if candidate.is_file() and is_linux_x86_elf(candidate):
            return candidate
    return None


def native_candidates(config: dict, runtime_dir: str) -> List[Path]:
    lay = layout(runtime_dir)
    configured = os.path.expanduser(str(config.get("rmsx_executable") or "").strip())
    found = ([Path(configured)] if configured else []) + [lay["bin"] / exe_name(), lay["src"] / "scan"]
    out: List[Path] = []
    for path in found:
        if path.is_file() and path not in out:
            out.append(path)
    return out


# ── runtime bundle (downloaded on first use) ────────────────────────────────
# The scanner's source, scoring matrices, query models and Linux binary are not
# shipped in the plugin repository. They are downloaded once from the project's
# server into ``external/rmsx/`` the first time run_from_scratch needs them.

def bundle_ready(runtime_dir: str) -> bool:
    """True when the runtime folder already holds the scanner source, scoring
    matrices and query models (so nothing needs downloading)."""
    src = layout(runtime_dir)["src"]
    return ((src / "mat" / "iso.mat").is_file()
            and (src / "Queries" / "reduced").is_dir()
            and (src / "main_scan.cc").is_file())


def _read_sha256(text: str) -> str:
    for token in (text or "").split():
        if re.fullmatch(r"[0-9a-fA-F]{64}", token):
            return token.lower()
    return ""


def _download_bundle(config: dict, root: Path, tmp: Path, progress: Progress) -> Tuple[Path, str]:
    """Download ``rmsx.tar.gz``, verify its SHA-256, and unpack it into ``tmp/unpacked``.

    Members are checked before anything is written: absolute or ``..`` paths are
    rejected and symbolic links / device files are skipped. Returns
    ``(unpacked_dir, url)``; raises on any failure.
    """
    url = str(config.get("rmsx_bundle_url") or DEFAULT_BUNDLE_URL).strip()
    archive = tmp / "rmsx.tar.gz"
    progress(f"Downloading the RNAMotifScanX files from {url} (about 30 MB)...")
    response, insecure = open_url(url, timeout=120)
    digest = hashlib.sha256()
    with response, open(archive, "wb") as out:
        for chunk in iter(lambda: response.read(1 << 20), b""):
            out.write(chunk)
            digest.update(chunk)
    if insecure:
        progress("Warning: this Python has no CA certificates, so the server certificate "
                 "was not verified (the checksum below still is).")

    sha_response, _ = open_url(url + ".sha256", timeout=60)
    with sha_response:
        expected = _read_sha256(sha_response.read().decode("utf-8", errors="replace"))
    if not expected:
        raise ValueError(f"could not read a SHA-256 checksum from {url}.sha256")
    if digest.hexdigest() != expected:
        raise ValueError("checksum mismatch: the downloaded archive is corrupt or was altered "
                         f"(expected {expected[:12]}..., got {digest.hexdigest()[:12]}...)")
    progress("Checksum verified. Unpacking...")

    staging = tmp / "unpacked"
    staging.mkdir()
    staging_resolved = staging.resolve()
    with tarfile.open(archive, "r:*") as tar:
        members = tar.getmembers()
        names = [Path(m.name).parts for m in members if Path(m.name).parts]
        strip = 1 if names and all(parts[0] == "rmsx" for parts in names) else 0
        keep = []
        for member in members:
            parts = Path(member.name).parts
            if len(parts) <= strip:
                continue
            if Path(member.name).is_absolute() or ".." in parts:
                raise ValueError(f"unsafe path in archive: {member.name}")
            if member.issym() or member.islnk() or member.isdev() or member.isfifo():
                continue                     # e.g. a stray libstdc++.a symlink
            member.name = "/".join(parts[strip:])
            target = (staging / member.name).resolve()
            if staging_resolved not in target.parents and target != staging_resolved:
                raise ValueError(f"unsafe path in archive: {member.name}")
            keep.append(member)
        tar.extractall(staging, members=keep)
    if not (staging / "RNAMotifScanX_src" / "mat").is_dir():
        raise ValueError("archive does not contain RNAMotifScanX_src/mat")
    return staging, url


def ensure_runtime_bundle(config: dict, runtime_dir: str, progress: Progress) -> dict:
    """Download and unpack ``rmsx.tar.gz`` into ``runtime_dir`` unless it is already there.

    The archive is verified against ``<url>.sha256`` before anything is unpacked
    (see :func:`_download_bundle`) and its files are moved into ``runtime_dir``
    without overwriting any file that already exists, so a rerun never clobbers
    your own files. Returns ``{'ok', 'downloaded', 'url', 'error'}``; never raises.
    """
    result = {"ok": False, "downloaded": False, "url": "", "error": ""}
    if bundle_ready(runtime_dir):
        result["ok"] = True
        return result

    result["url"] = str(config.get("rmsx_bundle_url") or DEFAULT_BUNDLE_URL).strip()
    root = layout(runtime_dir)["root"]
    tmp = None
    try:
        root.mkdir(parents=True, exist_ok=True)
        tmp = Path(tempfile.mkdtemp(prefix=".rmsx_dl_", dir=str(root)))
        progress(f"RNAMotifScanX runtime not found under {root}.")
        staging, _ = _download_bundle(config, root, tmp, progress)
        moved = 0
        for path in sorted(staging.rglob("*")):
            destination = root / path.relative_to(staging)
            if path.is_dir():
                destination.mkdir(parents=True, exist_ok=True)
            elif not destination.exists():
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(path), str(destination))
                moved += 1
        progress(f"RNAMotifScanX runtime unpacked into {root} ({moved} files).")
        result.update(ok=True, downloaded=True)
    except urllib.error.HTTPError as exc:
        result["error"] = f"HTTP {exc.code} {exc.reason} for {exc.url or result['url']}"
    except Exception as exc:  # noqa: BLE001 - reported to the caller
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)
    return result


# ── files macOS moved to iCloud ─────────────────────────────────────────────
# With "Desktop & Documents in iCloud" and "Optimize Mac Storage", macOS can
# replace files under ~/Desktop or ~/Documents by empty "dataless" placeholders
# when the disk is nearly full. Reading one blocks until iCloud brings it back,
# which can take forever (or never finish), and the scanner then hangs with
# PyMOL waiting on it. These helpers find such files and repair them.

_SF_DATALESS = 0x40000000


def is_dataless(path) -> bool:
    try:
        return bool(getattr(os.lstat(path), "st_flags", 0) & _SF_DATALESS)
    except OSError:
        return False


def wait_for_local_files(paths, timeout: float = 45.0) -> List[str]:
    """Make offloaded files readable again. Returns the ones still unavailable
    after ``timeout`` seconds (an empty list when everything is on disk)."""
    pending = [str(p) for p in dict.fromkeys(str(p) for p in paths) if is_dataless(p)]
    if not pending:
        return []

    def touch(path: str) -> None:
        try:
            with open(path, "rb") as fh:
                fh.read(1)
        except OSError:
            pass

    threads = []
    for path in pending:
        thread = threading.Thread(target=touch, args=(path,), daemon=True)
        thread.start()
        threads.append((path, thread))
    deadline = time.time() + timeout
    for _path, thread in threads:
        thread.join(max(0.0, deadline - time.time()))
    return [path for path, thread in threads if thread.is_alive() or is_dataless(path)]


def restore_runtime_files(config: dict, runtime_dir: str, files, progress: Progress) -> dict:
    """Replace files under the runtime folder (e.g. offloaded ones) with fresh
    copies from the verified ``rmsx.tar.gz``. Returns ``{'ok', 'restored', 'error'}``."""
    result = {"ok": False, "restored": 0, "error": ""}
    root = layout(runtime_dir)["root"].resolve()
    tmp = None
    try:
        tmp = Path(tempfile.mkdtemp(prefix=".rmsx_fix_", dir=str(root)))
        staging, _ = _download_bundle(config, root, tmp, progress)
        for file in files:
            try:
                relative = Path(file).resolve().relative_to(root)
            except ValueError:
                continue
            source = staging / relative
            if not source.is_file():
                continue
            destination = root / relative
            try:
                destination.unlink()
            except OSError:
                pass
            shutil.copy2(source, destination)
            result["restored"] += 1
        result["ok"] = True
    except Exception as exc:  # noqa: BLE001
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)
    return result


# ── docker / wsl ────────────────────────────────────────────────────────────

def docker_state() -> Tuple[str, bool, str]:
    """(docker path or '', daemon reachable, detail)."""
    docker = shutil.which("docker")
    if not docker:
        return "", False, "docker is not installed"
    rc, out, err = run_capture([docker, "info", "--format", "{{.ServerVersion}}"], timeout=25)
    if rc == 0 and out.strip():
        return docker, True, f"daemon reachable (server {out.strip()})"
    return docker, False, _last_line(err, "daemon not reachable")


def wsl_state() -> Tuple[str, bool, str]:
    """(wsl.exe path or '', a distro is usable, detail); Windows only."""
    if platform.system() != "Windows":
        return "", False, "not Windows"
    wsl = shutil.which("wsl") or shutil.which("wsl.exe")
    if not wsl:
        return "", False, "wsl.exe not found (WSL is not installed)"
    rc, out, err = run_capture([wsl, "-e", "sh", "-c", "echo $HOME"], timeout=45)
    if rc == 0 and out.strip().startswith("/"):
        return wsl, True, f"WSL distro usable (home {out.strip()})"
    return wsl, False, _last_line(err or out, "no WSL distribution installed or it did not start")


def to_wsl_path(path) -> str:
    text = str(path)
    match = re.match(r"^([A-Za-z]):[\\/](.*)$", text)
    if not match:
        return text.replace("\\", "/")
    return f"/mnt/{match.group(1).lower()}/{match.group(2).replace(chr(92), '/')}"


class _Mounts:
    """Maps host paths to container paths, bind-mounting each directory once."""

    def __init__(self) -> None:
        self.mounts: List[Tuple[Path, str]] = []

    def map(self, path, is_dir: bool = False) -> str:
        p = Path(path).resolve()
        base = p if is_dir else p.parent
        for host, target in self.mounts:
            if base == host or host in base.parents:
                rel = p.relative_to(host).as_posix()
                return target if rel == "." else f"{target}/{rel}"
        target = f"/mnt/m{len(self.mounts)}"
        self.mounts.append((base, target))
        return target if is_dir else f"{target}/{p.name}"

    def args(self) -> List[str]:
        out: List[str] = []
        for host, target in self.mounts:
            out += ["--mount", f"type=bind,source={host},target={target},readonly"]
        return out


@dataclass
class Runtime:
    kind: str            # 'native' | 'wsl' | 'docker'
    exe: str             # native: the executable; wsl/docker: the Linux ELF on the host
    note: str = ""
    wsl_exe: str = ""    # wsl: path of the binary inside the distro (copied there by setup)


def _wsl_prefix(rt: Runtime) -> List[str]:
    wsl = shutil.which("wsl") or shutil.which("wsl.exe") or "wsl.exe"
    return [wsl, "-e", rt.wsl_exe or to_wsl_path(rt.exe)]


def build_command(rt: Runtime, src_root: str, query: str, in_file: str, nch_file: str,
                  threads: int) -> Tuple[List[str], dict, str]:
    """Return ``(argv, env, container_name)`` for one scan on the given runtime.

    ``container_name`` is set for Docker so a timed-out run can be removed.
    """
    def tail(q: str, i: str, n: str) -> List[str]:
        return [q, i, f"--map_pdb={n}", "--pvalue", SCAN_PVALUE, "--num_threads", str(threads),
                "--write_alignment"]

    env = os.environ.copy()
    if rt.kind == "native":
        env["RNAMOTIFSCANX_PATH"] = str(src_root)
        return [rt.exe] + tail(query, in_file, nch_file), env, ""
    if rt.kind == "wsl":
        argv = _wsl_prefix(rt)[:2] + ["env", f"RNAMOTIFSCANX_PATH={to_wsl_path(src_root)}"] + _wsl_prefix(rt)[2:]
        argv += tail(to_wsl_path(query), to_wsl_path(in_file), to_wsl_path(nch_file))
        return argv, env, ""
    docker = shutil.which("docker") or "docker"
    mounts = _Mounts()
    ctr_src = mounts.map(src_root, is_dir=True)   # first, so the exe and queries reuse it
    ctr_exe = mounts.map(rt.exe)
    ctr = tail(mounts.map(query), mounts.map(in_file), mounts.map(nch_file))
    name = f"rmsx_scan_{os.getpid()}_{uuid.uuid4().hex[:8]}"
    argv = [docker, "run", "--rm", "--platform", "linux/amd64", "--name", name,
            "-e", f"RNAMOTIFSCANX_PATH={ctr_src}", *mounts.args(), DOCKER_IMAGE, ctr_exe, *ctr]
    return argv, env, name


# ── choosing a runtime ──────────────────────────────────────────────────────

def _try_native(config: dict, runtime_dir: str) -> Tuple[Optional[Runtime], str]:
    candidates = native_candidates(config, runtime_dir)
    if not candidates:
        return None, f"no scan executable found for this platform ({platform_dir()})"
    reasons = []
    for path in candidates:
        if not binary_matches_host(path):
            reasons.append(f"{path} is not a {platform.system()}/{host_machine()} executable")
            continue
        if os.name != "nt" and not os.access(path, os.X_OK):
            try:
                os.chmod(path, os.stat(path).st_mode | 0o111)
            except OSError:
                pass
        ok, detail = probe([str(path)])
        if ok:
            return Runtime("native", str(path), f"native scan: {path}"), ""
        reasons.append(f"{path}: {detail}")
    return None, "; ".join(reasons)


def _try_wsl(config: dict, runtime_dir: str, verify: bool) -> Tuple[Optional[Runtime], str]:
    wsl, usable, detail = wsl_state()
    if not usable:
        return None, detail
    elf = linux_elf(config, runtime_dir)
    if elf is None:
        return None, "no Linux x86-64 scan binary to run inside WSL"
    home = detail.split("home ")[-1].rstrip(")")
    copied = f"{home}/.rsmviewer/scan"
    rc, _, _ = run_capture([wsl, "-e", "test", "-x", copied], timeout=30)
    rt = Runtime("wsl", str(elf), f"WSL2: {copied if rc == 0 else to_wsl_path(elf)}",
                 wsl_exe=copied if rc == 0 else "")
    if verify:
        ok, why = probe(_wsl_prefix(rt))
        if not ok:
            return None, why
    return rt, ""


def _try_docker(config: dict, runtime_dir: str, verify: bool) -> Tuple[Optional[Runtime], str]:
    docker, up, detail = docker_state()
    if not docker:
        return None, detail
    if not up:
        return None, f"docker installed but {detail}"
    elf = linux_elf(config, runtime_dir)
    if elf is None:
        return None, "no Linux x86-64 scan binary to run inside the container"
    rt = Runtime("docker", str(elf), f"Docker ({DOCKER_IMAGE}, linux/amd64)")
    if verify:
        mounts = _Mounts()
        ctr_exe = mounts.map(elf)
        argv = [docker, "run", "--rm", "--platform", "linux/amd64", *mounts.args(), DOCKER_IMAGE, ctr_exe]
        ok, why = probe(argv, timeout=900)                 # the first run pulls the image
        if not ok:
            return None, why
    return rt, ""


def resolve_runtime(config: dict, runtime_dir: str, verify: bool = False) -> Tuple[Optional[Runtime], List[str]]:
    """First working runtime, or ``(None, reasons)``. ``config['scan_runtime']``
    may force ``native``, ``wsl`` or ``docker``; the default ``auto`` tries them in that order."""
    mode = str(config.get("scan_runtime") or "auto").strip().lower()
    order = {"auto": ["native", "wsl", "docker"], "native": ["native"], "wsl": ["wsl"], "docker": ["docker"]}
    reasons: List[str] = []
    if mode not in order:
        reasons.append(f"unknown scan_runtime '{mode}' (use auto, native, wsl or docker); using auto")
        mode = "auto"
    for kind in order[mode]:
        if kind == "native":
            rt, why = _try_native(config, runtime_dir)
        elif kind == "wsl":
            if platform.system() != "Windows":
                rt, why = None, "WSL exists only on Windows"
            else:
                rt, why = _try_wsl(config, runtime_dir, verify)
        else:
            rt, why = _try_docker(config, runtime_dir, verify)
        if rt:
            return rt, reasons
        reasons.append(f"{kind}: {why}")
    return None, reasons


# ── building from source ────────────────────────────────────────────────────

def _find_compiler() -> str:
    for name in (os.environ.get("CXX", ""), "clang++", "g++", "c++"):
        if name and shutil.which(name):
            return shutil.which(name)
    return ""


def find_boost() -> Tuple[Optional[Path], Optional[Path], List[str]]:
    """(include dir, lib dir, missing libs) for the first usable Boost, else (None, None, [])."""
    roots: List[Path] = []
    if os.environ.get("BOOST_ROOT"):
        roots.append(Path(os.environ["BOOST_ROOT"]))
    if shutil.which("brew"):
        rc, out, _ = run_capture([shutil.which("brew"), "--prefix", "boost"], timeout=60)
        if rc == 0 and out.strip():
            roots.append(Path(out.strip()))
    roots += [Path(p) for p in ("/opt/homebrew/opt/boost", "/usr/local/opt/boost", "/opt/local", "/usr", "/usr/local")]
    best: Tuple[Optional[Path], Optional[Path], List[str]] = (None, None, [])
    for root in roots:
        include = root / "include"
        if not (include / "boost" / "program_options.hpp").is_file():
            continue
        for lib in (root / "lib", root / "lib64", root / "lib" / "x86_64-linux-gnu", root / "lib" / "aarch64-linux-gnu"):
            if not lib.is_dir():
                continue
            missing = [n for n in BOOST_LIBS if not any(lib.glob(f"libboost_{n}.*"))]
            if not missing:
                return include, lib, []
            best = (include, lib, missing)
    return best


def _macos_sdks() -> List[str]:
    """SDK roots to try when linking; the default first, then older ones (a
    Command Line Tools SDK newer than the OS can have unreadable .tbd stubs)."""
    sdks: List[str] = [""]
    base = Path("/Library/Developer/CommandLineTools/SDKs")
    if base.is_dir():
        def version(p: Path):
            m = re.search(r"MacOSX(\d+(?:\.\d+)?)", p.name)
            return float(m.group(1)) if m else 0.0
        real = sorted({p.resolve() for p in base.glob("MacOSX*.sdk")}, key=version, reverse=True)
        sdks += [str(p) for p in real]
    return sdks


def build_native_scan(runtime_dir: str, progress: Progress, install_deps: bool = True) -> dict:
    """Compile ``scan`` from the bundled source into ``bin/<platform>/scan``."""
    result = {"ok": False, "path": "", "error": "", "log": ""}
    system = platform.system()
    if system not in ("Darwin", "Linux"):
        result["error"] = "building from source is only supported on macOS and Linux"
        return result
    src = layout(runtime_dir)["src"]
    threadpool = src / "Libraries" / "threadpool"
    missing = [n for n in SCAN_SOURCES if not (src / f"{n}.cc").is_file()]
    if missing or not threadpool.is_dir():
        result["error"] = f"RNAMotifScanX source is incomplete under {src} (missing: {', '.join(missing) or 'threadpool'})"
        return result

    cxx = _find_compiler()
    if not cxx:
        result["error"] = ("no C++ compiler found. " + (
            "Install the Xcode command line tools: xcode-select --install" if system == "Darwin"
            else "Install one, e.g.: sudo apt-get install g++"))
        return result
    progress(f"C++ compiler: {cxx}")

    include, libdir, missing_libs = find_boost()
    if (include is None or missing_libs) and system == "Darwin" and install_deps and shutil.which("brew"):
        progress("Boost is not installed; running: brew install boost  (this can take several minutes)")
        brew_cmd = [shutil.which("brew"), "install", "boost"]
        if is_translated():   # Homebrew refuses to run as x86_64 in its arm64 prefix
            brew_cmd = ["/usr/bin/arch", "-arm64"] + brew_cmd
        rc, out, err = run_capture(brew_cmd, timeout=3600)
        if rc != 0:
            result["error"] = f"brew install boost failed: {_last_line(err or out)}"
            return result
        include, libdir, missing_libs = find_boost()
    if include is None or missing_libs:
        need = ("brew install boost" if system == "Darwin"
                else "sudo apt-get install libboost-all-dev zlib1g-dev")
        result["error"] = (f"Boost libraries not found"
                           + (f" (missing: {', '.join(missing_libs)})" if missing_libs else "")
                           + f". Install them with: {need}")
        return result
    progress(f"Boost: {include.parent} (libs in {libdir})")
    arch_flags = ["-arch", "arm64"] if system == "Darwin" and host_machine() == "arm64" else []

    bin_dir = layout(runtime_dir)["bin"]
    bin_dir.mkdir(parents=True, exist_ok=True)
    log_lines: List[str] = []
    build_dir = Path(tempfile.mkdtemp(prefix="rmsx_build_"))
    try:
        def compile_one(name: str) -> Tuple[str, int, str]:
            obj = build_dir / f"{name}.o"
            argv = [cxx] + arch_flags + ["-std=c++14", "-O2", "-w", f"-I{include}", f"-I{threadpool}",
                    "-c", str(src / f"{name}.cc"), "-o", str(obj)]
            rc, out, err = run_capture(argv, timeout=900)
            return name, rc, (err or out)

        progress(f"Compiling {len(SCAN_SOURCES)} source files...")
        workers = max(1, min(4, os.cpu_count() or 1))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            for name, rc, text in pool.map(compile_one, SCAN_SOURCES):
                log_lines.append(f"[compile {name}] exit {rc}\n{text}")
                if rc != 0:
                    result["error"] = f"compiling {name}.cc failed: {_last_line(text)}"
                    return result

        libs = [f"-lboost_{n}" for n in BOOST_LIBS]
        if any(libdir.glob("libboost_system.*")):
            libs.append("-lboost_system")            # header-only in Boost >= 1.69 on some builds
        libs += ["-lz", "-lm", "-lpthread"] + (["-lrt"] if system == "Linux" else [])
        objects = [str(build_dir / f"{n}.o") for n in SCAN_SOURCES]
        target = bin_dir / exe_name()
        linked = False
        for sdk in (_macos_sdks() if system == "Darwin" else [""]):
            argv = [cxx] + arch_flags + ["-std=c++14", "-O2", "-w"] + (["-isysroot", sdk] if sdk else []) + objects + [
                "-o", str(target), f"-L{libdir}", f"-Wl,-rpath,{libdir}"] + libs
            progress("Linking" + (f" (SDK {Path(sdk).name})" if sdk else "") + "...")
            rc, out, err = run_capture(argv, timeout=900)
            log_lines.append(f"[link sdk={sdk or 'default'}] exit {rc}\n{err or out}")
            if rc == 0:
                linked = True
                break
        if not linked:
            result["error"] = f"linking failed: {_last_line(log_lines[-1])}"
            return result
        os.chmod(target, 0o755)
        ok, detail = probe([str(target)])
        if not ok:
            result["error"] = f"built {target} but it does not start: {detail}"
            return result
        result.update(ok=True, path=str(target))
        return result
    finally:
        result["log"] = "\n".join(log_lines)
        try:
            (bin_dir / "build.log").write_text(result["log"], encoding="utf-8")
        except OSError:
            pass
        shutil.rmtree(build_dir, ignore_errors=True)


# ── setup (rmv_setup RNAMotifScanX) ─────────────────────────────────────────

def _start_colima(progress: Progress) -> bool:
    colima = shutil.which("colima")
    if not colima:
        return False
    progress("Docker is not running; starting Colima (this takes about a minute)...")
    rc, out, err = run_capture([colima, "start"], timeout=900)
    if rc != 0:
        progress(f"colima start failed: {_last_line(err or out)}")
        return False
    return docker_state()[1]


def setup(config: dict, runtime_dir: str, progress: Progress, install_deps: bool = True) -> dict:
    """Make sure a scanner runtime works on this machine.

    Order: an existing native scan; a native build from source (macOS/Linux);
    WSL2 (Windows); Docker. Returns ``{'ok', 'runtime', 'problems', 'next'}``.
    """
    system = platform.system()
    report = {"ok": False, "runtime": None, "problems": [], "next": []}
    progress(f"Platform: {system} {host_machine()} ({platform_dir()})"
             + ("; this PyMOL is an Intel build running under Rosetta" if is_translated() else ""))

    bundle = ensure_runtime_bundle(config, runtime_dir, progress)
    if not bundle["ok"]:
        report["problems"].append(f"could not download the RNAMotifScanX runtime: {bundle['error']}")
        report["next"].append(f"Download {bundle['url']} yourself and extract it so that "
                              f"{layout(runtime_dir)['src']} exists (see external/rmsx_setup.md)")
        return report

    rt, why = _try_native(config, runtime_dir)
    if rt:
        report.update(ok=True, runtime=rt)
        return report
    progress(f"No ready-to-run native scanner ({why}).")

    if system in ("Darwin", "Linux"):
        progress("Building the scanner from the bundled source...")
        built = build_native_scan(runtime_dir, progress, install_deps=install_deps)
        if built["ok"]:
            rt, _ = _try_native(config, runtime_dir)
            report.update(ok=rt is not None, runtime=rt)
            if rt:
                return report
        else:
            report["problems"].append(f"native build: {built['error']}")
            progress(f"Native build not possible: {built['error']}")

    if system == "Windows":
        wsl, usable, detail = wsl_state()
        elf = linux_elf(config, runtime_dir)
        if usable and elf:
            home = detail.split("home ")[-1].rstrip(")")
            dest = f"{home}/.rsmviewer/scan"
            progress(f"Copying the Linux scanner into WSL ({dest})...")
            run_capture([wsl, "-e", "sh", "-c",
                         f'mkdir -p "{home}/.rsmviewer" && cp "{to_wsl_path(elf)}" "{dest}" && chmod +x "{dest}"'],
                        timeout=120)
            rt, why = _try_wsl(config, runtime_dir, verify=True)
            if rt:
                report.update(ok=True, runtime=rt)
                return report
            report["problems"].append(f"wsl: {why}")
        else:
            report["problems"].append(f"wsl: {detail}")
            report["next"].append("Install WSL2 (PowerShell as Administrator): wsl --install   then reboot, and run rmv_setup RNAMotifScanX again")

    docker, up, detail = docker_state()
    if docker and not up and system in ("Darwin", "Linux"):
        up = _start_colima(progress)
        detail = "daemon reachable" if up else detail
    if docker and up:
        progress(f"Checking the Docker route (first run downloads the {DOCKER_IMAGE} image)...")
        rt, why = _try_docker(config, runtime_dir, verify=True)
        if rt:
            report.update(ok=True, runtime=rt)
            return report
        report["problems"].append(f"docker: {why}")
    else:
        report["problems"].append(f"docker: {detail}")
        if system == "Darwin":
            report["next"].append("Or install Docker: brew install docker colima   then: colima start")
        elif system == "Windows":
            report["next"].append("Or install Docker Desktop: https://www.docker.com/products/docker-desktop/")
        else:
            report["next"].append("Or install Docker: https://docs.docker.com/engine/install/")
    return report


# ── diagnostics (rmv_rmsx_doctor) ───────────────────────────────────────────

def diagnose(config: dict, runtime_dir: str, pdb_id: str = "") -> dict:
    """Structured health report: ``{'lines': [(status, label, detail)], 'ready': bool}``.

    status is one of ``ok``, ``warn``, ``fail``, ``info``.
    """
    lines: List[Tuple[str, str, str]] = []
    add = lambda status, label, detail="": lines.append((status, label, detail))
    system, machine = platform.system(), host_machine()
    lay = layout(runtime_dir)

    add("info", "Platform", f"{system} {machine} ({platform_dir()})")
    if is_translated():
        add("info", "Rosetta", "this PyMOL is an Intel build running under Rosetta; "
            "the scanner is built and run natively for arm64")
    add("info", "Data mode", str(config.get("data_mode", "preannotated")))

    # Scanner runtimes -----------------------------------------------------
    add("info", "--- Scanner (needed only for run_from_scratch) ---")
    if bundle_ready(runtime_dir):
        add("ok", "RNAMotifScanX runtime files", str(lay["root"]))
    else:
        add("warn", "RNAMotifScanX runtime files",
            "not downloaded yet; downloaded automatically from "
            f"{config.get('rmsx_bundle_url') or DEFAULT_BUNDLE_URL} the first time you run from scratch "
            "(or run: rmv_setup RNAMotifScanX)")
    native, why = _try_native(config, runtime_dir)
    add("ok" if native else "warn", "Native scanner", native.exe if native else why)
    if system == "Windows":
        wsl, usable, detail = wsl_state()
        add("ok" if usable else "warn", "WSL2", detail)
    docker, up, detail = docker_state()
    add("ok" if up else "warn", "Docker", detail)
    if docker and not up and shutil.which("colima"):
        add("info", "Colima", "installed but not running: start it with `colima start`")
    elf = linux_elf(config, runtime_dir)
    add("ok" if elf else "warn", "Linux scan binary (for WSL/Docker)", str(elf) if elf else f"not found under {lay['src']}")

    rt, reasons = resolve_runtime(config, runtime_dir)
    if rt:
        add("ok", "Selected runtime", f"{rt.kind}: {rt.note}")
    else:
        add("fail", "Selected runtime", "none available" + ("".join(f"\n      - {r}" for r in reasons)))
        add("info", "Fix", "run: rmv_setup RNAMotifScanX")

    # Build toolchain ------------------------------------------------------
    if system in ("Darwin", "Linux"):
        add("info", "--- Native build toolchain (used by rmv_setup RNAMotifScanX) ---")
        cxx = _find_compiler()
        add("ok" if cxx else "warn", "C++ compiler", cxx or "not found")
        include, libdir, missing = find_boost()
        if include and not missing:
            add("ok", "Boost", f"{include.parent}")
        else:
            add("warn", "Boost", "not found" + (f" (missing {', '.join(missing)})" if missing else ""))

    # Data -----------------------------------------------------------------
    add("info", "--- Data ---")
    src = find_src_root(runtime_dir)
    add("ok" if src else "fail", "Scoring matrices (mat/)", str(src / "mat") if src else f"not found under {lay['src']}")
    qdirs = query_dirs(config, runtime_dir)
    families = list(config.get("motif_families") or [])
    if src:
        offloaded = [q for q in (find_query(f, qdirs) for f in families) if q and is_dataless(q)]
        offloaded += [str(p) for p in (src / "mat").glob("*") if is_dataless(p)]
        if offloaded:
            add("warn", "Files offloaded to iCloud",
                f"{len(offloaded)} required file(s), e.g. {offloaded[0]}; macOS moved them off disk "
                "(disk nearly full). They are restored automatically on the next run, or free up "
                "disk space / choose 'Download Now' in Finder")
    if qdirs:
        gaps = [f for f in families if not find_query(f, qdirs)]
        add("ok" if not gaps else "warn", "Query models", f"{qdirs[0]}" + (f" (missing: {', '.join(gaps)})" if gaps else ""))
    else:
        add("fail", "Query models", "no query directory found")
    prebuild = Path(os.path.expanduser(str(config.get("pdb_prebuild_dir") or ""))) if config.get("pdb_prebuild_dir") else None
    add("ok" if prebuild and prebuild.is_dir() else "warn", "Prepared-input / results folder", str(prebuild or "not configured"))
    if pdb_id and prebuild:
        entry = prebuild / pdb_id.lower()
        chains = sorted(d.name for d in entry.iterdir() if d.is_dir() and not d.name.startswith("_")) if entry.is_dir() else []
        add("ok" if chains else "info", f"{pdb_id.upper()} local data", f"chains {', '.join(chains)}" if chains else "not downloaded yet")

    base = str(config.get("preannotated_base_url") or "").rstrip("/")
    if base:
        try:
            response, insecure = open_url(f"{base}/{(pdb_id or '1s72').lower()}.tar.gz", timeout=10, method="HEAD")
            response.close()
            add("warn" if insecure else "ok", "Results server", "reachable" + (" (TLS not verified: this Python has no CA certificates)" if insecure else ""))
        except urllib.error.HTTPError as exc:
            add("ok" if exc.code == 404 else "warn", "Results server", f"reachable (HTTP {exc.code} for this PDB)")
        except Exception as exc:  # noqa: BLE001
            add("warn", "Results server", f"unreachable: {type(exc).__name__}: {exc}")

    ready = rt is not None and src is not None and bool(qdirs)
    return {"lines": lines, "ready": ready, "runtime": rt}
