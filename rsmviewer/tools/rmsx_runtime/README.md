# RNAMotifScanX Runtime Integration in RSMViewer

This directory contains the integrated RNAMotifScanX runtime support used by RSMViewer Source 7.

RNAMotifScanX (RMSX) is a structural motif search tool. RSMViewer uses it to scan an RNA structure against curated query motifs, convert the search output into RSMViewer-readable annotation files, and visualize the detected motifs in PyMOL through the normal RSMViewer command set.

Source 7 makes RMSX usable from inside RSMViewer without requiring the user to manually assemble configuration files, locate query folders, or wire output directories. The plugin discovers bundled runtime assets, runs the pipeline, and loads the generated motif annotations into the same viewer used for all other sources.

---

## Contents

1. [What Source 7 Does](#what-source-7-does)
2. [Quick Tutorial](#quick-tutorial)
3. [Runtime Folder Structure](#runtime-folder-structure)
4. [Important Files and Their Functions](#important-files-and-their-functions)
5. [How RSMViewer Integrates RMSX](#how-rsmviewer-integrates-rmsx)
6. [Pipeline Implementation](#pipeline-implementation)
7. [Runtime Setup and Doctor Commands](#runtime-setup-and-doctor-commands)
8. [Platform Binary Layout](#platform-binary-layout)
9. [Build System](#build-system)
10. [Modifications Added for RSMViewer](#modifications-added-for-rsmviewer)
11. [Generated Files](#generated-files)
12. [Troubleshooting](#troubleshooting)
13. [Developer Checklist](#developer-checklist)

---

## What Source 7 Does

RSMViewer exposes RNAMotifScanX as Source 7:

```pymol
rmv_db 7
```

When Source 7 is active, `rmv_load_motif` runs the integrated RMSX pipeline for the current structure and loads the result logs into RSMViewer. The supported motif families are:

- `k-turn`
- `c-loop`
- `sarcin-ricin`
- `reverse-kturn`
- `e-loop`

The output can then be inspected with normal RSMViewer commands:

```pymol
rmv_summary
rmv_show K-TURN
rmv_show SARCIN-RICIN
rmv_view K-TURN
rmv_save K-TURN cif
```

Source 7 is strict local-runtime mode:

- it uses bundled query and matrix files
- it discovers runtime executables from `bin/<platform>/`
- it runs fresh when invoked through Source 7
- it does not require an external JSON config for normal use
- it does not download a CIF during Source-7 execution
- it fails with explicit doctor/setup messages when runtime assets are missing

---

## Quick Tutorial

### 1. Diagnose the runtime

Start PyMOL, load RSMViewer, then run:

```pymol
rmv_rmsx_doctor
```

Expected ready output includes:

- runtime directory
- platform directory
- RMSX executable path
- MC-Annotate executable path
- `Status: OK`

If the status is `NOT READY`, run setup:

```pymol
rmv_rmsx setup
rmv_rmsx_doctor
```

### 2. Load a structure

```pymol
rmv_fetch 1S72
```

### 3. Select Source 7

```pymol
rmv_db 7
```

### 4. Run and load motifs

```pymol
rmv_load_motif
```

This calls the integrated Source-7 pipeline. RSMViewer prepares the runtime configuration internally, runs RMSX, then loads generated `result_0_100_withbs.log` files from the RNAMotifScanX output directory.

### 5. Inspect results

```pymol
rmv_summary
rmv_show K-TURN
rmv_show C-LOOP
rmv_show E-LOOP
```

### 6. Explicit runtime commands

You can also run Source 7 explicitly:

```pymol
rmv_rmsx run 1S72
rmv_rmsx run_current
```

`rmv_load_motif`, `rmv_rmsx run`, and `rmv_rmsx run_current` all use the integrated runtime path when Source 7 is available.

---

## Runtime Folder Structure

```text
rmsx_runtime/
|-- README.md
|-- setup_runtime.py
|-- bin/
|   |-- macos-arm64/
|   |   |-- scan
|   |   |-- MC-Annotate
|   |   `-- libmccore.2.0.0.dylib
|   |-- macos-x86_64/
|   |-- linux-x86_64/
|   `-- windows-x86_64/
|-- mat/
|   |-- iso.mat
|   |-- nuc.mat
|   `-- stk.mat
|-- queries/
|   |-- c-loop_consensus.struct
|   |-- e-loop_consensus.struct
|   |-- k-turn_consensus.struct
|   |-- reverse-kturn_consensus.struct
|   `-- sarcin-ricin_consensus.struct
|-- src/
|   |-- RNAMotifScanX_src/
|   `-- external/
|       |-- MC-Annotate/
|       `-- mccore/
`-- build/
	`-- <platform>/
```

Some local/generated folders may exist during development, such as `build/<platform>/` or `src/external/*-build-*`. Those are build artifacts and should not be treated as required source files.

---

## Important Files and Their Functions

### `setup_runtime.py`

First-run setup and doctor helper for Source 7. It performs four main jobs:

1. detects the current platform folder name
2. validates required runtime assets
3. optionally builds RMSX from source with CMake
4. reports readiness as JSON for `gui.py`

Platform names are normalized as:

- macOS Apple Silicon: `macos-arm64`
- macOS Intel: `macos-x86_64`
- Linux Intel/AMD 64-bit: `linux-x86_64`
- Windows Intel/AMD 64-bit: `windows-x86_64`

The doctor checks for:

- RMSX executable
- MC-Annotate executable
- required query files
- matrix files

### `bin/<platform>/`

Contains platform-specific executables. RSMViewer searches this folder first for the current platform, then falls back to other bin folders only as a best-effort discovery path.

Recognized RMSX executable names:

- `RNAMotifScanX`
- `scan`
- `RNAMotifScanX.exe`
- `scan.exe`

Recognized MC-Annotate executable names:

- `MC-Annotate`
- `mc-annotate`
- `MC-Annotate.exe`
- `mc-annotate.exe`

### `queries/`

Contains the motif query structures used by RMSX. Each file corresponds to one supported motif family:

- `k-turn_consensus.struct`
- `c-loop_consensus.struct`
- `sarcin-ricin_consensus.struct`
- `reverse-kturn_consensus.struct`
- `e-loop_consensus.struct`

`rmsx_runner.py` maps family names to these query files and runs RMSX once per family.

### `mat/`

Contains scoring matrix files used by RNAMotifScanX:

- `iso.mat`
- `nuc.mat`
- `stk.mat`

During RMSX execution, `rmsx_runner.py` sets `RNAMOTIFSCANX_PATH` to the runtime root when `mat/` is present so the scanner can locate these matrices.

### `src/RNAMotifScanX_src/`

Vendored RNAMotifScanX source tree. RSMViewer added a CMake build path here so setup can build the `scan` executable on supported platforms when dependencies are available.

Important source/build files include:

- `CMakeLists.txt`: CMake project added for RSMViewer packaging.
- `main_scan.cc`: RMSX scan entry point used to build `scan`.
- `structural_motif.*`, `scan_structure.*`, `motif_graph_matching.*`: core RMSX motif search logic.
- `Queries/`: upstream query files included with the source tree.
- `mat/`: upstream matrix files included with the source tree.
- `StructureAnnotation/`: upstream annotation-related files and reference sequence data.
- `StructureAnnotation/pdb_seqres.na.fa`: sequence reference used when converting MC-Annotate output into RMSX input format.

### `src/external/MC-Annotate/`

Vendored MC-Annotate source tree used for source provenance and optional source builds. It was added as normal tracked files, not a Git submodule or gitlink.

### `src/external/mccore/`

Vendored MCCORE source tree. MC-Annotate depends on MCCORE. This tree was also added as normal tracked files.

### `build/<platform>/`

Local CMake build output generated by setup. This is machine-specific and should not be committed.

---

## How RSMViewer Integrates RMSX

The PyMOL-facing integration is implemented in `rsmviewer/gui.py`.

Important methods:

- `_build_internal_rmsx_config(...)`
- `_run_rmsx_runtime_setup(...)`
- `ensure_rmsx_runtime_ready(...)`
- `rmsx_doctor(...)`
- `run_rmsx_wrapper(...)`

### Internal config generation

RSMViewer builds Source-7 config automatically from bundled plugin paths. Normal users do not need to provide `rmsx_pipeline_config.json`.

The internal config includes:

- `rmsx_executable`: detected from `bin/<platform>/`
- `mc_annotate_executable`: detected from `bin/<platform>/`
- `query_motifs_dir`: `rmsx_runtime/queries/`
- `output_dir`: `rsmviewer/database/user_annotations/RNAMotifScanX/`
- `auto_download_cif`: `False` in Source-7 flow
- `motif_families`: default five RMSX families
- `target_chains`: default `['0']`
- `max_strands`: `3`
- `num_threads`: `4`

### Source selection

The user selects Source 7 with:

```pymol
rmv_db 7
```

If a user passes an external JSON config to Source 7, RSMViewer warns that Source 7 uses the integrated runtime config and ignores the external JSON for normal Source-7 execution.

### Runtime readiness

Before execution, RSMViewer calls `ensure_rmsx_runtime_ready(auto_setup=True)`. This runs `setup_runtime.py --json` to validate assets. If something is missing, setup is attempted once with `--build`.

---

## Pipeline Implementation

The runtime pipeline is implemented in `rsmviewer/tools/rmsx_runner.py`.

High-level flow:

```text
rmv_db 7
	-> select Source 7

rmv_load_motif
	-> load_user_annotations_action('rnamotifscanx', pdb)
	-> ensure_rmsx_runtime_ready(auto_setup=True)
	-> build internal RMSX config
	-> rmsx_runner.run_pipeline(..., force_fresh=True)
	-> run MC-Annotate on the structure
	-> run RNAVIEW on the structure
	-> merge MC-Annotate and RNAVIEW interactions (union, MC-Annotate precedence)
	-> convert the merged annotation into .rmsx.in target files
	-> run RMSX once per motif family
	-> write result_0_100_withbs.log files
	-> load RMSX motifs into RSMViewer
```

### Step 1: Runtime validation

`setup_runtime.py` verifies that required binaries, query files, and matrix files exist.

### Step 2: Structure location

`rmsx_runner.py` receives a PDB ID and optional structure path. In Source-7 mode, RSMViewer sets `auto_download_cif=False`, so the structure must already be available through the RSMViewer workflow.

### Step 3: MC-Annotate preparation

MC-Annotate is run on a PDB-formatted structure file. Its output is saved as:

```text
rsmviewer/database/user_annotations/RNAMotifScanX/<PDB_ID>_mc_annotate.out
```

The runner parses residues, base pairs, and stacking interactions from this file.

### Step 3b: RNAVIEW annotation

RNAVIEW is run on the same structure to produce a second base-pair annotation. Its output is merged with the MC-Annotate result as a union (MC-Annotate takes precedence on shared residue keys), reproducing the reference `PrepareInput.py` behavior. If RNAVIEW is unavailable, the runner falls back to MC-Annotate-only annotation with a diagnostic message (unless `rnaview_required` is set).

### Step 4: RMSX input generation

The merged MC-Annotate and RNAVIEW output is converted into one or more `.rmsx.in` target files. RNAVIEW base pairs are unioned with the MC-Annotate base pairs (MC-Annotate takes precedence on shared residue keys), matching the reference `PrepareInput.py` behavior. These files encode:

- sequence header
- reference sequence
- base-pair interactions
- stacking interactions

The reference sequence is read from the upstream FASTA file when available. If no reference sequence is found, the parsed structure sequence is used as a fallback so the pipeline can continue.

### Step 5: Per-family RMSX search

The runner scans each configured family:

```text
k-turn        -> k-turn_consensus.struct
c-loop        -> c-loop_consensus.struct
sarcin-ricin  -> sarcin-ricin_consensus.struct
reverse-kturn -> reverse-kturn_consensus.struct
e-loop        -> e-loop_consensus.struct
```

For each family, output is written to:

```text
rsmviewer/database/user_annotations/RNAMotifScanX/<family>_consensus/result_0_100_withbs.log
```

### Step 6: RSMViewer loading

After logs are produced, RSMViewer loads them as Source-7 user annotations. From that point onward, they behave like other motif sources inside the viewer.

---

## Runtime Setup and Doctor Commands

### PyMOL commands

```pymol
rmv_rmsx_doctor
rmv_rmsx setup
rmv_rmsx status
rmv_rmsx test
rmv_rmsx run <PDB_ID>
rmv_rmsx run_current
```

### What `rmv_rmsx_doctor` reports

The doctor prints:

- runtime directory
- detected platform directory
- binary directory used
- RMSX executable path
- MC-Annotate executable path
- missing runtime assets
- missing query files
- setup/build messages

### CLI equivalent

From the repository root:

```bash
python rsmviewer/tools/rmsx_runtime/setup_runtime.py \
	--runtime-dir rsmviewer/tools/rmsx_runtime \
	--json
```

Attempt setup/build:

```bash
python rsmviewer/tools/rmsx_runtime/setup_runtime.py \
	--runtime-dir rsmviewer/tools/rmsx_runtime \
	--build \
	--json
```

### Linux CI build and packaging

For release packaging on a real Linux environment, use the repository workflow:

- `.github/workflows/build-linux-rmsx-runtime.yml`

This workflow runs on Ubuntu and:

1. installs build dependencies
2. runs `scripts/ci/build_linux_rmsx_runtime.sh`
3. builds and validates `scan`, `cut`, `align`, `rnaview`, and `MC-Annotate`
4. verifies runtime link dependencies with `ldd` and fails on unresolved libraries
5. uploads `dist/rmsx_runtime_linux_x86_64.tar.gz` as a downloadable artifact

The CI script also emits `dist/reports/ldd_summary.txt` so release validation records
exactly which shared libraries each binary expects on Ubuntu.

You can also run the same script manually on Linux:

```bash
chmod +x scripts/ci/build_linux_rmsx_runtime.sh
./scripts/ci/build_linux_rmsx_runtime.sh
```

Optional source build with MC-Annotate/MCCORE fetch/build allowed:

```bash
python rsmviewer/tools/rmsx_runtime/setup_runtime.py \
	--runtime-dir rsmviewer/tools/rmsx_runtime \
	--build \
	--fetch-mc-annotate \
	--fetch-mccore \
	--json
```

The optional fetch flags are setup-time developer/release tools. Normal Source-7 runtime execution remains strict.

---

## Platform Binary Layout

Runtime binaries are discovered under:

```text
bin/<platform>/
```

Expected platform folders:

- `bin/macos-arm64/`
- `bin/macos-x86_64/`
- `bin/linux-x86_64/`
- `bin/windows-x86_64/`

Current bundled executable set:

```text
bin/macos-arm64/scan
bin/macos-arm64/MC-Annotate
bin/macos-arm64/libmccore.2.0.0.dylib
```

For full cross-platform release coverage, matching native executables still need to be supplied for:

- `macos-x86_64`
- `linux-x86_64`
- `windows-x86_64`

Each platform folder should contain both:

- an RMSX scanner executable (`scan` or `scan.exe`)
- an MC-Annotate executable (`MC-Annotate`, `mc-annotate`, or `.exe` equivalent)

On macOS, any required dynamic library such as `libmccore.2.0.0.dylib` should be colocated in the same platform bin folder and the executable rpath should point to the executable directory.

---

## Build System

RSMViewer added CMake-based build support around the RNAMotifScanX source tree:

```text
src/RNAMotifScanX_src/CMakeLists.txt
```

`rmv_rmsx setup` calls `setup_runtime.py --build`, which uses this CMake project to build `scan` or `scan.exe` when possible.

Build dependencies:

- CMake
- a C++ compiler for the current platform
- Boost headers/libraries (`program_options`, `thread`, `system`, `filesystem`, `iostreams`)
- zlib

Typical platform dependency examples:

- macOS: Xcode command line tools, CMake, Boost, zlib
- Linux: GCC/Clang, CMake, `libboost-all-dev`, zlib development headers
- Windows: MSVC build tools, CMake, Boost from vcpkg/conan/manual install, zlib

The setup helper detects common Boost locations, including:

- `BOOST_ROOT`
- `BOOST_HOME`
- `/opt/homebrew`
- `/usr/local`
- `/usr`
- common Windows/vcpkg paths

### RNAVIEW build

Setup also builds RNAVIEW from the bundled C source at
`src/RNAMotifScanX_src/ThirdParty/RNAVIEW/` and installs the resulting binary into
`bin/<platform>/rnaview`. This build is independent of Boost and CMake:

- Linux and macOS are fully supported (compiled with `cc`/`clang`/`gcc`).
- Windows is attempted only when a MinGW/GCC toolchain is available; MSVC is not supported.
- The build applies audit-validated flags
  (`-D_FORTIFY_SOURCE=0 -Wno-implicit-function-declaration -Wno-implicit-int -Wno-return-type`)
  so the source compiles cleanly on modern toolchains without being modified.
- An existing platform-compatible `rnaview` binary is reused; rebuilding is skipped unless forced.
- If the RNAVIEW build fails, setup continues and the pipeline falls back gracefully.

RNAVIEW composes its BASEPARS resource path into a fixed-size internal buffer, so at
runtime `rmsx_runner.py` transparently stages a short path to the RNAVIEW directory
when the install path would be too long. RNAVIEW source is never modified.

---

## Modifications Added for RSMViewer

The RMSX integration is not just a copied binary folder. Several packaging and integration changes were added:

1. Integrated Source-7 command routing in `rsmviewer/gui.py`.
2. Internal runtime config generation so users do not need external JSON files.
3. Runtime doctor/setup commands through `rmv_rmsx_doctor` and `rmv_rmsx setup`.
4. Platform-aware executable discovery under `bin/<platform>/`.
5. Fresh-run behavior for Source 7 through `force_fresh=True`.
6. Strict Source-7 behavior with `auto_download_cif=False` during RSMViewer runs.
7. CMake build support for RNAMotifScanX source.
8. Optional MC-Annotate/MCCORE source build support for release engineering.
9. macOS MCCORE runtime library colocation/rpath support.
10. Conversion from the merged MC-Annotate and RNAVIEW output into `.rmsx.in` target files.
11. RMSX CLI invocation with a short-option fallback for compatibility across RMSX builds.
12. Automatic RNAVIEW build from bundled C source during setup, with a runtime short-path workaround for RNAVIEW's fixed BASEPARS buffer.
13. Git packaging cleanup so MC-Annotate and MCCORE source trees are tracked as normal files, not embedded gitlinks.

---

## Generated Files

During execution, Source 7 may generate files under:

```text
rsmviewer/database/user_annotations/RNAMotifScanX/
```

Common generated files:

- `<PDB_ID>_mc_annotate.out`
- `<PDB_ID>_<chain>.rmsx.in`
- `<family>_consensus/result_0_100_withbs.log`

During setup/build, local build artifacts may be generated under:

```text
rsmviewer/tools/rmsx_runtime/build/<platform>/
rsmviewer/tools/rmsx_runtime/src/external/*-build-*/
rsmviewer/tools/rmsx_runtime/src/external/*-install-*/
```

These build directories are machine-specific. They are not required for users who receive packaged binaries and should not be committed as release source documentation.

---

## Troubleshooting

### `RMSX runtime is incomplete`

Run:

```pymol
rmv_rmsx_doctor
rmv_rmsx setup
rmv_rmsx_doctor
```

Check the missing items printed by doctor. The usual causes are:

- missing `scan` executable
- missing `MC-Annotate` executable
- missing query files
- missing matrix files
- unsupported platform binary set

### Missing platform binaries

Doctor reports the expected platform directory. Put native executables into:

```text
rsmviewer/tools/rmsx_runtime/bin/<platform>/
```

Then rerun:

```pymol
rmv_rmsx_doctor
```

### `scan built successfully, but no compatible MC-Annotate binary is available`

This means RNAMotifScanX built, but setup could not find a native MC-Annotate executable for the current platform.

Fix options:

- provide a native MC-Annotate executable in `bin/<platform>/`
- build MC-Annotate and MCCORE from source during release packaging
- use the optional setup flags for development builds

### Boost not found

Install Boost for the current platform, then rerun setup.

Examples:

```bash
brew install boost
```

or on Debian/Ubuntu:

```bash
sudo apt install libboost-all-dev zlib1g-dev cmake build-essential
```

### No RMSX result files after run

Check:

- `rmv_fetch <PDB_ID>` was run first
- Source 7 was selected with `rmv_db 7`
- `rmv_rmsx_doctor` reports `Status: OK`
- generated output directory is writable
- MC-Annotate produced a non-empty `<PDB_ID>_mc_annotate.out`
- query files exist in `queries/`

Then retry:

```pymol
rmv_rmsx run_current
rmv_summary
```

### MC-Annotate produces empty output

MC-Annotate typically expects PDB-formatted structure input. If annotation fails:

- confirm the PDB structure file exists for the active PDB
- check the MC-Annotate executable runs on the current platform
- inspect the generated output directory for `<PDB_ID>_mc_annotate.out`

### Old result files appear reused

Source 7 passes `force_fresh=True` during the integrated flow. The runner removes old family logs before rerunning. If you see old output, check whether you are using a standalone CLI path without `--fresh`.

Standalone fresh run:

```bash
python rsmviewer/tools/rmsx_runner.py --config <config.json> --pdb 1S72 --fresh
```

### External JSON config ignored

This is expected for normal Source 7. RSMViewer builds the Source-7 runtime config internally from bundled assets. External JSON config files are optional developer paths, not the default user workflow.

---

## Developer Checklist

Before calling Source 7 release-ready, verify:

1. `rmv_rmsx_doctor` reports `Status: OK` on the target platform.
2. `bin/<platform>/` contains native `scan` and MC-Annotate executables.
3. `queries/` contains the five expected consensus query files.
4. `mat/` contains matrix files.
5. `rmv_fetch 1S72`, `rmv_db 7`, `rmv_load_motif`, and `rmv_summary` run without manual path edits.
6. Generated logs are written under `database/user_annotations/RNAMotifScanX/`.
7. Build directories and local output files are not committed accidentally.
8. `src/external/MC-Annotate/` and `src/external/mccore/` remain normal tracked source directories, not nested Git repositories.

Useful Git packaging checks:

```bash
git ls-tree -r --full-tree HEAD rsmviewer/tools/rmsx_runtime/src | grep '^160000'
git check-ignore -v rsmviewer/tools/rmsx_runtime/bin/macos-arm64/scan || true
git status --short
```

The first command should print nothing. If it prints a `160000` entry, a nested Git repository was committed as a gitlink instead of as real source files.
