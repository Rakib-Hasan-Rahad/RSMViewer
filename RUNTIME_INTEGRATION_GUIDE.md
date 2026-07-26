# RSMViewer Runtime Integration Guide

How the two external structural-annotation engines — **RNAMotifScanX (RMSX, Source 7)**
and **FR3D (Source 5)** — were taken from their raw, platform-limited source form and
turned into fully integrated, one-command features inside the RSMViewer PyMOL plugin.

This document is written for both developers and curious users. It explains, in plain
language and with examples:

1. what each engine is,
2. what raw materials we started from,
3. every step we took to integrate them,
4. how a user actually uses them,
5. what a user needs to know, and
6. the biggest challenges we hit and how we solved them.

> **Reading order:** Part 1 covers RNAMotifScanX (the harder one, built from raw C/C++).
> Part 2 covers FR3D (a Python toolkit). If you only care about *using* the tools, jump
> to [How a User Uses RMSX](#16-how-a-user-uses-rmsx-step-by-step) and
> [How a User Uses FR3D](#25-how-a-user-uses-fr3d-step-by-step).

---

## Table of Contents

- [Part 1 — RNAMotifScanX (Source 7)](#part-1--rnamotifscanx-source-7)
  - [1.1 What RNAMotifScanX Is](#11-what-rnamotifscanx-is)
  - [1.2 The Starting Point: Raw Materials](#12-the-starting-point-raw-materials)
  - [1.3 The Goal](#13-the-goal)
  - [1.4 Architecture Overview](#14-architecture-overview)
  - [1.5 Step-by-Step: How We Integrated It](#15-step-by-step-how-we-integrated-it)
  - [1.6 How a User Uses RMSX (Step by Step)](#16-how-a-user-uses-rmsx-step-by-step)
  - [1.7 What a User Needs to Know](#17-what-a-user-needs-to-know)
  - [1.8 Biggest Challenges We Faced](#18-biggest-challenges-we-faced)
- [Part 2 — FR3D (Source 5)](#part-2--fr3d-source-5)
  - [2.1 What FR3D Is](#21-what-fr3d-is)
  - [2.2 The Starting Point: Raw Materials](#22-the-starting-point-raw-materials)
  - [2.3 Architecture Overview](#23-architecture-overview)
  - [2.4 Step-by-Step: How We Integrated It](#24-step-by-step-how-we-integrated-it)
  - [2.5 How a User Uses FR3D (Step by Step)](#25-how-a-user-uses-fr3d-step-by-step)
  - [2.6 What a User Needs to Know](#26-what-a-user-needs-to-know)
  - [2.7 Biggest Challenges We Faced](#27-biggest-challenges-we-faced)
- [Appendix A — Command Cheat Sheet](#appendix-a--command-cheat-sheet)
- [Appendix B — File & Directory Map](#appendix-b--file--directory-map)

---

# Part 1 — RNAMotifScanX (Source 7)

## 1.1 What RNAMotifScanX Is

**RNAMotifScanX (RMSX)** is a research tool that searches an RNA 3D structure for known
structural motif families (like k-turns, sarcin-ricin loops, c-loops, e-loops, and
reverse k-turns). Instead of matching sequence letters, it matches the *interaction
network* of a structure — which bases pair with which, and which bases stack on which.

RMSX does not read a PDB/mmCIF file directly. It reads a pre-digested description of the
structure's base pairs and stacking interactions. Producing that description is a
multi-tool preprocessing job, which is the core of what we had to integrate.

The pipeline, at the highest level, is:

```text
3D structure  ──►  annotate interactions  ──►  format for RMSX  ──►  RMSX search  ──►  motif hits
   (PDB)          (MC-Annotate + RNAVIEW)        (.rmsx.in)          (scan binary)      (result logs)
```

## 1.2 The Starting Point: Raw Materials

What we received was **not a ready-to-run program**. It was a bundle of source code and
data that only compiled and ran on Linux. Concretely:

| Raw material | Language | Original state | Role |
|---|---|---|---|
| **RNAMotifScanX source** (`main_scan.cc`, `structural_motif.*`, `scan_structure.*`, `motif_graph_matching.*`, …) | C++ | Linux-only, hand Makefile, needs Boost + zlib | The actual motif search engine → builds the `scan` binary |
| **MC-Annotate source** | C++ | Linux-only | Annotates base pairs + stacking from a structure |
| **MCCORE source** | C++ | Linux-only, dependency of MC-Annotate | Core library MC-Annotate links against |
| **RNAVIEW source** (14 `.c` files) | C (late-1990s ANSI/K&R) | Linux-only, plain `Makefile`, shipped with a **Linux 32-bit ELF binary** | A *second* base-pair annotator; its output is merged with MC-Annotate |
| **RNAVIEW `BASEPARS/`** | data (11 `Atomic_*.pdb` files) | reference geometry | RNAVIEW loads these at runtime to classify pairs |
| **Query motifs** (`*_consensus.struct`) | data | 5 families | The motif templates RMSX searches for |
| **Scoring matrices** (`iso.mat`, `nuc.mat`, `stk.mat`) | data | — | Scoring tables RMSX needs |
| **`pdb_seqres.na.fa`** | data | reference | Reference sequences used when formatting RMSX input |
| Upstream `PrepareInput.py`, `ParseStructureAnnotation.py` | Python **2** | reference scripts | Show *how* upstream merged MC-Annotate + RNAVIEW into `.rmsx.in` |

Key problems baked into these raw materials:

- **Linux-only.** The bundled `scan`, MC-Annotate, and RNAVIEW binaries were Linux ELF.
  On macOS/Windows they simply cannot run.
- **The RNAVIEW binary was 32-bit Linux ELF** — unusable on modern macOS (arm64/x86_64)
  and Windows.
- **The reference preprocessing scripts were Python 2** — they raise `TabError`/syntax
  errors under the Python 3 that PyMOL ships, so we could not just call them.
- **No unified configuration.** Running RMSX by hand meant writing JSON config, locating
  query folders, wiring output directories, and setting environment variables.

## 1.3 The Goal

Turn all of the above into a single user action:

```pymol
rmv_db 7
rmv_load_motif
```

…that works on Linux and macOS (and best-effort on Windows), builds whatever binaries it
can from source on first run, and **falls back gracefully** with clear messages when a
piece is missing — without the user ever writing a config file.

## 1.4 Architecture Overview

```text
                     PyMOL  (user types rmv_db 7 ; rmv_load_motif)
                        │
                        ▼
          rsmviewer/gui.py   ── Source-7 command routing
            ├─ ensure_rmsx_runtime_ready()  ─►  runs setup_runtime.py --json [--build]
            ├─ _build_internal_rmsx_config() ─►  builds config from bundled paths (no JSON needed)
            └─ run_rmsx_wrapper()            ─►  calls the pipeline runner
                        │
                        ▼
   rsmviewer/tools/rmsx_runner.py   ── the actual pipeline
            1. run MC-Annotate  ─────────────►  <PDB>_mc_annotate.out
            2. run RNAVIEW      ─────────────►  <PDB>.pdb.out
            3. merge (union, MC-Annotate wins) ─►  merged interactions
            4. write <PDB>_<chain>.rmsx.in
            5. run scan once per motif family ─►  result_0_100_withbs.log
                        │
                        ▼
          loaded back into RSMViewer as Source-7 motifs
```

Supporting build logic lives in
[rsmviewer/tools/rmsx_runtime/setup_runtime.py](rsmviewer/tools/rmsx_runtime/setup_runtime.py):

- `detect_platform_dir()` → `macos-arm64` / `macos-x86_64` / `linux-x86_64` / `windows-x86_64`
- `build_from_source()` → builds the `scan` binary via CMake
- `build_mccore_from_repo()` / `build_mc_annotate_from_repo()` → optional source builds
- `build_rnaview_from_source()` / `find_rnaview_binary()` → the RNAVIEW C build
- `validate_runtime()` → the "doctor" health check

## 1.5 Step-by-Step: How We Integrated It

### Step A — Vendor the source trees as normal files

We placed every upstream tree *inside* the plugin as ordinary tracked files (not Git
submodules, so there are no fragile network dependencies at clone time):

```text
rsmviewer/tools/rmsx_runtime/
├── src/
│   ├── RNAMotifScanX_src/          ← the RMSX engine source
│   │   ├── main_scan.cc, structural_motif.*, scan_structure.*, ...
│   │   ├── CMakeLists.txt          ← ADDED by us (upstream had only a Makefile)
│   │   ├── StructureAnnotation/    ← pdb_seqres.na.fa + annotation helpers
│   │   └── ThirdParty/RNAVIEW/     ← RNAVIEW C source (+ BASEPARS/, include/, src/)
│   └── external/
│       ├── MC-Annotate/            ← MC-Annotate source
│       └── mccore/                 ← MCCORE source (MC-Annotate depends on this)
├── queries/                        ← 5 consensus query motifs
├── mat/                            ← iso.mat, nuc.mat, stk.mat
└── bin/<platform>/                 ← where built binaries land
```

**Why:** a self-contained tree means the plugin can rebuild what it needs offline, and a
`git clone` never breaks because of a missing submodule.

### Step B — Add a CMake build for the `scan` engine

Upstream RMSX shipped a hand-written Makefile tuned for one Linux box. We added
`src/RNAMotifScanX_src/CMakeLists.txt` so the same source builds portably. `setup_runtime.py`'s
`build_from_source()` invokes CMake to produce the platform `scan` binary.

Build dependencies (documented so a user can install them):

- CMake, a C++ compiler
- **Boost** (`program_options`, `thread`, `system`, `filesystem`, `iostreams`)
- zlib

The setup helper auto-detects Boost in common locations (`BOOST_ROOT`, `/opt/homebrew`,
`/usr/local`, `/usr`, vcpkg paths).

Example (what setup runs under the hood, from the repo root):

```bash
python rsmviewer/tools/rmsx_runtime/setup_runtime.py \
    --runtime-dir rsmviewer/tools/rmsx_runtime \
    --build --json
```

### Step C — Build MC-Annotate and MCCORE (optional, release-time)

MC-Annotate is the primary interaction annotator, and it links against MCCORE. Because
these are heavy C++ builds, we made them **opt-in** for release engineering rather than
forcing them on every user:

```bash
python rsmviewer/tools/rmsx_runtime/setup_runtime.py \
    --runtime-dir rsmviewer/tools/rmsx_runtime \
    --build --fetch-mccore --fetch-mc-annotate --json
```

On macOS, MCCORE produces a dynamic library (`libmccore.2.0.0.dylib`). We colocate it in
the same `bin/<platform>/` folder as the executables and set the executable rpath to the
executable's own directory so it loads without extra environment setup.

### Step D — Build RNAVIEW automatically from its C source

This was the newest and trickiest piece. The shipped RNAVIEW binary was Linux 32-bit ELF,
so on macOS/Windows we had to compile RNAVIEW from its 14 `.c` files. We added
`build_rnaview_from_source()` to `setup_runtime.py`:

- **Source list** (compiled in Makefile order):
  `rnaview.c fpair.c fpair_sub.c pair_type.c nrutil.c ps-xy.c ps-xy-sub.c vrml.c
  rnaxml-new.c analyze.c pattern.c xml2ps.c multiple.c statistics.c`
- **Validated compiler flags** (found by a build-feasibility audit — see challenges):
  ```text
  -D_FORTIFY_SOURCE=0 -Wno-implicit-function-declaration -Wno-implicit-int -Wno-return-type
  ```
  RNAVIEW is late-1990s C, so modern clang/gcc reject its implicit declarations by
  default; these flags let it compile **without modifying a single line of RNAVIEW source**.
- **Compiler selection:** Linux/macOS use `cc`/`clang`/`gcc`. Windows is attempted **only**
  with a MinGW/GCC toolchain; **MSVC is deliberately not supported** (it is not practical
  for this codebase).
- **Skip-if-present:** `find_rnaview_binary()` checks for an existing platform-compatible
  `rnaview` and skips rebuilding unless forced.
- **Non-fatal:** if the build fails, setup continues and the pipeline falls back to
  MC-Annotate-only annotation.

The single compile command it runs looks like:

```bash
cc -D_FORTIFY_SOURCE=0 -Wno-implicit-function-declaration -Wno-implicit-int -Wno-return-type \
   -I <RNAVIEW>/include <14 .c files> -o bin/<platform>/rnaview -lm
```

CLI switches we added: `--rebuild-rnaview` (force) and `--skip-rnaview` (bypass).

### Step E — Reproduce the preprocessing pipeline (MC-Annotate + RNAVIEW → `.rmsx.in`)

The reference `PrepareInput.py`/`ParseStructureAnnotation.py` scripts were Python 2, so we
could not call them. Instead we **ported their exact logic** into
[rsmviewer/tools/rmsx_runner.py](rsmviewer/tools/rmsx_runner.py), and verified the port is
**byte-for-byte faithful** with a differential test. The pipeline:

1. **Run MC-Annotate** on the structure → `<PDB>_mc_annotate.out`; parse residues, base
   pairs, and stacking interactions.
2. **Run RNAVIEW** on the same structure → `<PDB>.pdb.out`; parse its base pairs
   (`_parse_rnaview_output`).
3. **Merge** the two annotations as a **union** with **MC-Annotate precedence** on shared
   residue keys (`_merge_interactions`) — exactly matching the upstream behavior. RNAVIEW
   residue IDs are formatted as `<chain><resnum>` so their keys line up with MC-Annotate's.
4. **Write** the merged result into one or more `<PDB>_<chain>.rmsx.in` files (sequence
   header, reference sequence, base pairs, stacking). The reference sequence comes from
   `pdb_seqres.na.fa` when available, otherwise the parsed structure sequence is used.
5. **Run `scan` once per family**, writing `result_0_100_withbs.log` for each.

If RNAVIEW is unavailable, the runner logs a clear diagnostic and proceeds with
MC-Annotate-only annotation (unless the config sets `rnaview_required: true`).

### Step F — Platform binary layout and discovery

All built binaries are discovered from a predictable place first, with a best-effort
fallback to other platform folders:

```text
bin/macos-arm64/   bin/macos-x86_64/   bin/linux-x86_64/   bin/windows-x86_64/
```

Recognized names — RMSX engine: `RNAMotifScanX`, `scan` (`.exe` on Windows). MC-Annotate:
`MC-Annotate`, `mc-annotate` (`.exe`). RNAVIEW: `rnaview` (`.exe`).

### Step G — Generate the runtime config internally (no JSON for users)

`gui.py::_build_internal_rmsx_config()` assembles the whole configuration from bundled
plugin paths at runtime, so a normal user never writes JSON. It fills in:

- `rmsx_executable`, `mc_annotate_executable`, `rnaview_executable` (discovered in `bin/`)
- `rnaview_dir` (the `ThirdParty/RNAVIEW` folder that holds `BASEPARS/`)
- `query_motifs_dir` → `queries/`, `output_dir` → `database/user_annotations/RNAMotifScanX/`
- `motif_families` (the default five), plus threading and strand limits
- `incorporate_rnaview` (default **on**), `rnaview_required` (default **off**)

If a user passes an external JSON config, Source 7 warns that it uses the integrated
config and ignores the external file for normal execution.

### Step H — Wire it to PyMOL commands (Source 7)

Finally, we exposed everything through PyMOL commands. `rmv_db 7` selects the source, and
`rmv_load_motif` runs the pipeline and loads results just like any other source. Diagnostic
and control commands: `rmv_rmsx_doctor`, `rmv_rmsx setup|status|test|run|run_current`.

## 1.6 How a User Uses RMSX (Step by Step)

```pymol
# 1. Check the runtime is healthy (first run may build binaries)
rmv_rmsx_doctor

# 2. If it says NOT READY, run setup, then re-check
rmv_rmsx setup
rmv_rmsx_doctor

# 3. Load a structure
rmv_fetch 1S72

# 4. Select RNAMotifScanX as the data source
rmv_db 7

# 5. Run the search and load the detected motifs
rmv_load_motif

# 6. Inspect and visualize the results
rmv_summary
rmv_show K-TURN
rmv_show SARCIN-RICIN
rmv_view C-LOOP
rmv_save K-TURN cif
```

The five families you can find are: **k-turn, c-loop, sarcin-ricin, reverse-kturn, e-loop.**

You can also run it explicitly without changing the active source:

```pymol
rmv_rmsx run 1S72
rmv_rmsx run_current
```

## 1.7 What a User Needs to Know

- **First run may compile things.** On a fresh machine, `rmv_rmsx setup` can build the
  `scan` engine (needs CMake + Boost + zlib) and RNAVIEW (needs a C compiler). Subsequent
  runs reuse the built binaries.
- **Platform support.** Linux and macOS are fully supported. Windows is best-effort and
  needs a MinGW/GCC toolchain (MSVC is not supported for RNAVIEW).
- **It runs fresh each time.** Source 7 re-runs the pipeline and clears old family logs, so
  results always reflect the current structure — no stale cache surprises.
- **RNAVIEW is optional but recommended.** With it, base-pair coverage matches the
  reference pipeline. Without it, you still get MC-Annotate-only results plus a clear note.
- **No config files needed.** The plugin builds the configuration for you from bundled
  assets.
- **Output location.** Generated files live under
  `rsmviewer/database/user_annotations/RNAMotifScanX/`
  (`<PDB>_mc_annotate.out`, `<PDB>_<chain>.rmsx.in`, `<family>_consensus/result_0_100_withbs.log`).

If something fails, run `rmv_rmsx_doctor` — it prints exactly which binary, query, or
matrix file is missing.

## 1.8 Biggest Challenges We Faced

1. **The RNAVIEW stack-overflow crash (the hardest bug).**
   After we could build RNAVIEW, the binary still aborted with exit code 134
   (`SIGABRT`) at real install paths. An `lldb` backtrace pointed to
   `__stack_chk_fail` inside `get_reference_pdb()`. The cause: RNAVIEW composes its
   resource path into a **fixed-size `char spdb[80]` buffer** via
   `sprintf(spdb, "%sAtomic_%c.pdb", BDIR, ...)`. Our deep install path (~130+ characters)
   overflowed that buffer every time. Early "successful" test builds had only worked
   because they lived under a short `/tmp/...` path — a classic red herring that made us
   first suspect compiler flags and the compilation model.
   **Fix (without touching RNAVIEW source):** `rmsx_runner.py::_stage_short_rnaview_dir()`
   detects when the path would overflow and transparently exposes the RNAVIEW directory
   through a short staging path (a symlink, or a copy of `BASEPARS/` where symlinks aren't
   available), sets the `RNAVIEW` environment variable to it, runs the binary, and cleans
   up afterward. Verified end-to-end: a 137-char directory auto-staged and produced valid
   base-pair output.

2. **Modern compilers rejecting 1990s C.**
   RNAVIEW's implicit function declarations and implicit `int` are hard errors on current
   clang/gcc. We ran a **build-feasibility audit** to find the minimal flag set
   (`-D_FORTIFY_SOURCE=0 -Wno-implicit-function-declaration -Wno-implicit-int
   -Wno-return-type`) that compiles and *runs* correctly — with the constraint of **not
   editing the RNAVIEW source at all**. `-D_FORTIFY_SOURCE=0` was specifically required on
   macOS to avoid a benign overlapping-`strcpy` fortify abort in `rna()`.

3. **Python 2 reference scripts.**
   The upstream preprocessing (`PrepareInput.py`, `ParseStructureAnnotation.py`) is Python
   2 and won't even parse under PyMOL's Python 3. We re-implemented the merge logic in
   Python 3 and proved it byte-for-byte identical to the originals with a differential test,
   rather than shelling out to unusable scripts.

4. **Linux-only binaries in a cross-platform plugin.**
   Every shipped binary (`scan`, MC-Annotate, RNAVIEW) was Linux ELF. We had to add a
   portable CMake build for `scan`, optional source builds for MC-Annotate/MCCORE, and a
   from-source RNAVIEW build — each with platform detection and graceful fallback.

5. **Committed stale build artifacts.**
   The RNAVIEW tree shipped with stale Linux `.o` files in `obj/`. A plain `make` links
   those wrong-platform objects and fails. Our builder compiles the `.c` sources directly
   into the platform `bin/` and ignores the committed objects entirely.

6. **Heavy, optional dependencies.**
   Building MC-Annotate/MCCORE needs Boost and a full C++ toolchain. Forcing that on every
   user would be hostile, so we made those builds opt-in flags and kept the everyday path
   lightweight, with the doctor command explaining any gaps.

---

# Part 2 — FR3D (Source 5)

## 2.1 What FR3D Is

**FR3D ("Find RNA 3D")** is a geometry-based RNA annotation toolkit written in Python. It
inspects a 3D structure and reports loops and motifs directly from geometry. RSMViewer uses
it to generate local loop annotations for the loaded structure and then visualize them with
the standard command set.

Unlike RMSX, FR3D is Python, so there is **no C/C++ compilation** — but it has its own
integration challenge: Python dependencies and environment isolation.

## 2.2 The Starting Point: Raw Materials

| Raw material | Language | Original state | Role |
|---|---|---|---|
| **fr3d-python source tree** | Python 3 | a source package, not installed | The FR3D engine (`fr3d/…`) |
| Python dependencies | — | not present | `numpy`, `scipy`, `mmcif-pdbx` (imported as `pdbx`) |

The core problem: FR3D needs a Python environment with specific scientific packages. We
cannot assume PyMOL's bundled Python has `scipy` and `pdbx`, and we must not pollute the
user's system Python.

## 2.3 Architecture Overview

```text
                     PyMOL  (user types rmv_db 5 ; rmv_load_motif)
                        │
                        ▼
          rsmviewer/gui.py   ── Source-5 command routing
            ├─ ensure_fr3d_runtime_ready()  ─►  runs setup_runtime.py (build venv)
            ├─ _build_internal_fr3d_config()
            └─ run_fr3d_wrapper() → _run_fr3d_local_pipeline()
                        │
                        ▼
   rsmviewer/tools/fr3d_loop_extractor.py   ── runs FR3D, converts output
                        │
                        ▼
        <PDB>_fr3d_loops.csv  ─►  loaded into RSMViewer as Source-5 motifs
```

## 2.4 Step-by-Step: How We Integrated It

### Step A — Vendor the FR3D source tree

The FR3D source lives inside the plugin as normal tracked files:

```text
rsmviewer/tools/fr3d_runtime/
├── setup_runtime.py        ← setup + doctor helper
├── python/
│   ├── fr3d-python/         ← vendored FR3D source (tracked)
│   └── venv/                ← created locally on first run (NOT tracked)
├── bin/  cache/  scripts/
└── README.md
```

### Step B — Create an isolated runtime environment on first run

`setup_runtime.py::_bootstrap_runtime()` builds a **local virtual environment** so FR3D and
its dependencies never touch the user's system Python:

1. Find a host Python (`FR3D_BOOTSTRAP_PYTHON` env, then `sys.executable`, then `python3`).
2. Create `python/venv/` once (`python -m venv`).
3. Upgrade `pip`, then install `numpy`, `scipy`, and `mmcif-pdbx`.
4. Install the vendored FR3D package into the venv in editable mode
   (`pip install -e python/fr3d-python`).

The venv is intentionally **git-ignored** because it is machine-specific and large.

> **Release policy — strict and offline:** setup **never clones FR3D from the network** and
> **never fetches external source**. It only prepares the Python environment around the
> already-vendored source. If the vendored tree is missing, it reports that clearly instead
> of downloading anything.

### Step C — Detect readiness with a doctor check

`_detect_fr3d_root()` confirms `fr3d/__init__.py` exists, and `_detect_python()` confirms a
Python interpreter that can `import scipy, pdbx`. The report is emitted as JSON for `gui.py`
and surfaced through `rmv_fr3d doctor` / `rmv_fr3d status`.

### Step D — Run the pipeline and convert output

`fr3d_loop_extractor.py` runs FR3D against the loaded structure and converts its output into
RSMViewer's loop CSV format (`<PDB>_fr3d_loops.csv`), which then flows through the normal
loading path.

### Step E — Strict Source-5 execution policy

To keep behavior reproducible for reviewers and users, Source 5:

- always runs the bundled local FR3D pipeline (no remote/cache substitution),
- removes prior output for the requested PDB before rerunning,
- never falls back to BGSU download behavior in the Source-5 flow,
- keeps generated annotations under `rsmviewer/database/user_annotations/fr3d/`.

### Step F — Wire it to PyMOL commands (Source 5)

`rmv_db 5` selects the source; `rmv_load_motif` runs FR3D and loads the loops. Control and
diagnostics: `rmv_fr3d doctor|setup|status|run|run_current|refresh|config`.

## 2.5 How a User Uses FR3D (Step by Step)

```pymol
# 1. Check the runtime (first run may create a venv and install packages)
rmv_fr3d doctor

# 2. If not ready, run setup, then re-check
rmv_fr3d setup
rmv_fr3d status

# 3. Load a structure
rmv_fetch 1S72

# 4. Select FR3D as the data source
rmv_db 5

# 5. Run FR3D and load detected loops
rmv_load_motif

# 6. Inspect and visualize
rmv_summary
rmv_summary HL
rmv_show HL
rmv_view HL
rmv_save HL cif
```

Run it explicitly for a specific or current structure:

```pymol
rmv_fr3d run 1S72
rmv_fr3d run_current
rmv_fr3d refresh 1S72     # force a completely fresh rerun
```

## 2.6 What a User Needs to Know

- **First run installs packages.** `rmv_fr3d setup` creates a local virtual environment and
  installs `numpy`, `scipy`, and `mmcif-pdbx`. This needs the ability to create a venv and
  reach PyPI once. After that, it runs offline.
- **Your system Python is untouched.** Everything installs into
  `rsmviewer/tools/fr3d_runtime/python/venv/`.
- **It always runs fresh.** Source 5 re-runs FR3D and clears prior output for that PDB.
- **Write access is required** to `rsmviewer/tools/fr3d_runtime/` and
  `rsmviewer/database/user_annotations/fr3d/`.
- **Advanced override exists but is rarely needed:** `rmv_fr3d config <FR3D_ROOT>
  [OUTPUT_DIR] [AUTO_ON_FETCH]` points at a custom FR3D tree; running `rmv_fr3d setup`
  returns to the bundled runtime.

## 2.7 Biggest Challenges We Faced

1. **Environment isolation without breaking the user's Python.**
   FR3D needs `scipy`/`pdbx`, which PyMOL's Python may not have. Installing globally would
   risk the user's setup. Solution: a dedicated, git-ignored virtual environment built on
   first run, with FR3D installed editable inside it.
2. **Dependency detection vs. installation.**
   We needed to distinguish "environment is ready" from "environment must be built."
   `_python_has_fr3d_deps()` probes `import scipy, pdbx` and only triggers a build when
   required, so warm runs are instant.
3. **Reproducible, strict behavior.**
   Reviewers needed Source 5 to mean *local FR3D execution*, not a mixed local/remote/cache
   path. We enforced a strict policy: fresh run each time, no remote fallback, output always
   in the same place.
4. **Keeping the vendored tree clean.**
   The FR3D source is tracked as plain files (no nested `.git`, no submodule), and the venv
   plus generated outputs are excluded, so clones stay reproducible and the tree stays lean.

---

# Appendix A — Command Cheat Sheet

| Task | RNAMotifScanX (Source 7) | FR3D (Source 5) |
|---|---|---|
| Health check | `rmv_rmsx_doctor` | `rmv_fr3d doctor` |
| First-time setup | `rmv_rmsx setup` | `rmv_fr3d setup` |
| Status | `rmv_rmsx status` | `rmv_fr3d status` |
| Select source | `rmv_db 7` | `rmv_db 5` |
| Run + load | `rmv_load_motif` | `rmv_load_motif` |
| Explicit run (PDB) | `rmv_rmsx run 1S72` | `rmv_fr3d run 1S72` |
| Explicit run (current) | `rmv_rmsx run_current` | `rmv_fr3d run_current` |
| Force fresh | (Source 7 always runs fresh) | `rmv_fr3d refresh 1S72` |

Common visualization commands (both sources): `rmv_summary`, `rmv_show <TYPE>`,
`rmv_view <TYPE>`, `rmv_super <TYPE>`, `rmv_save <TYPE> cif`.

# Appendix B — File & Directory Map

**RNAMotifScanX (Source 7)**

| Path | Purpose |
|---|---|
| `rsmviewer/tools/rmsx_runtime/setup_runtime.py` | Setup, doctor, and all build logic (`scan`, MC-Annotate, MCCORE, RNAVIEW) |
| `rsmviewer/tools/rmsx_runner.py` | The pipeline: MC-Annotate + RNAVIEW → merge → `.rmsx.in` → `scan`; short-path staging |
| `rsmviewer/tools/rmsx_runtime/src/RNAMotifScanX_src/` | RMSX engine source + our `CMakeLists.txt` |
| `rsmviewer/tools/rmsx_runtime/src/RNAMotifScanX_src/ThirdParty/RNAVIEW/` | RNAVIEW C source + `BASEPARS/` |
| `rsmviewer/tools/rmsx_runtime/src/external/{MC-Annotate,mccore}/` | MC-Annotate + MCCORE source |
| `rsmviewer/tools/rmsx_runtime/queries/` | Five consensus query motifs |
| `rsmviewer/tools/rmsx_runtime/mat/` | Scoring matrices (`iso.mat`, `nuc.mat`, `stk.mat`) |
| `rsmviewer/tools/rmsx_runtime/bin/<platform>/` | Built binaries (`scan`, `MC-Annotate`, `rnaview`, libs) |
| `rsmviewer/database/user_annotations/RNAMotifScanX/` | Generated outputs |

**FR3D (Source 5)**

| Path | Purpose |
|---|---|
| `rsmviewer/tools/fr3d_runtime/setup_runtime.py` | Setup + doctor; builds the local venv |
| `rsmviewer/tools/fr3d_runtime/python/fr3d-python/` | Vendored FR3D source |
| `rsmviewer/tools/fr3d_runtime/python/venv/` | Local virtual environment (not tracked) |
| `rsmviewer/tools/fr3d_loop_extractor.py` | Runs FR3D and converts output to loop CSV |
| `rsmviewer/database/user_annotations/fr3d/` | Generated outputs |

**PyMOL-facing integration (both):** `rsmviewer/gui.py`
(`ensure_rmsx_runtime_ready`, `_build_internal_rmsx_config`, `run_rmsx_wrapper`,
`ensure_fr3d_runtime_ready`, `_build_internal_fr3d_config`, `run_fr3d_wrapper`,
`_run_fr3d_local_pipeline`).
