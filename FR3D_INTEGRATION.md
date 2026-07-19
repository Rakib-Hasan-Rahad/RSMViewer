# RSMViewer: FR3D Integration

Comprehensive reference for FR3D runtime integration in RSMViewer v1.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Runtime Architecture](#2-runtime-architecture)
3. [Installation Guide](#3-installation-guide)
4. [User Workflow](#4-user-workflow)
5. [Internal Pipeline](#5-internal-pipeline)
6. [Runtime Commands](#6-runtime-commands)
7. [Troubleshooting](#7-troubleshooting)
8. [Developer Notes](#8-developer-notes)

---

## 1. Overview

FR3D (Find RNA 3D) is a geometry-based RNA structural annotation toolkit.

RSMViewer integrates FR3D as **Source 5** to provide structure-driven loop annotations directly inside PyMOL. Source 5 loads FR3D-generated motif data into the same visualization, summary, superimposition, and export workflow used by all other sources.

Core goals of the integration:

- Keep Source 5 usage simple: `rmv_db 5`, `rmv_fetch <PDB_ID>`, `rmv_load_motif`
- Run a local FR3D pipeline from within RSMViewer
- Keep runtime setup reproducible in a bundled runtime directory
- Provide clear diagnostics via doctor/status commands

---

## 2. Runtime Architecture

### 2.1 Runtime directory

FR3D runtime assets are managed under:

- `rsmviewer/tools/fr3d_runtime/`

Key runtime paths:

- `rsmviewer/tools/fr3d_runtime/setup_runtime.py` — doctor/setup helper
- `rsmviewer/tools/fr3d_runtime/python/venv/` — runtime virtual environment (created by setup)
- `rsmviewer/tools/fr3d_runtime/python/fr3d-python/` — FR3D checkout (created by setup)
- `rsmviewer/database/user_annotations/fr3d/` — Source 5 output and data directory

### 2.2 Automatic runtime discovery

RSMViewer discovers FR3D runtime components by checking:

- runtime-local FR3D checkout (`python/fr3d-python`)
- environment variables when provided (`FR3D_ROOT`, `FR3D_PYTHON`)
- Python executables that can import required dependencies (`scipy`, `pdbx`)

### 2.3 Automatic setup

If runtime requirements are missing, Source 5 can bootstrap automatically through setup/doctor paths:

- create runtime venv
- install Python dependencies
- clone FR3D source
- install FR3D package into runtime environment

Commands:

- `rmv_fr3d setup`
- `rmv_fr3d doctor`

### 2.4 Cache behavior and fallback behavior

Source 5 is strict local-pipeline mode:

- local FR3D pipeline is required
- no BGSU fallback execution path
- fresh pipeline execution for each run
- existing `*_fr3d_loops.csv` for the same PDB is removed before a new run

This behavior applies to:

- `rmv_load_motif` when Source 5 is selected
- `rmv_fr3d run`
- `rmv_fr3d run_current`

---

## 3. Installation Guide

### 3.1 Prerequisites

- PyMOL with RSMViewer installed
- Internet access for first-time runtime bootstrap (dependency install and FR3D clone)
- A usable Python interpreter available on the host machine

### 3.2 First-time setup (recommended)

Run in PyMOL:

```pymol
rmv_fr3d setup
rmv_fr3d doctor
```

Expected doctor result:

- `Status: LOCAL PIPELINE READY`
- runtime mode reported as `local_pipeline`

### 3.3 Verify setup in the normal workflow

```pymol
rmv_db 5
rmv_fetch 1S72
rmv_load_motif
rmv_summary
```

Expected behavior:

- Source 5 starts a fresh local FR3D run from bundled runtime assets
- output is written to `rsmviewer/database/user_annotations/fr3d/`
- motifs become available to `rmv_summary`, `rmv_show`, `rmv_view`, `rmv_super`, `rmv_save`

### 3.4 Optional advanced configuration

Advanced users can override Source 5 runtime paths when using a custom FR3D tree:

```pymol
rmv_fr3d config /absolute/path/to/fr3d-python /absolute/path/to/output on
```

Syntax:

- `rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]`

---

## 4. User Workflow

Standard Source 5 workflow:

```pymol
rmv_db 5
rmv_fetch 1S72
rmv_load_motif
rmv_summary
rmv_show HAIRPIN LOOP
```

Step-by-step:

1. `rmv_db 5` selects FR3D integration source.
2. `rmv_fetch 1S72` loads the structure.
3. `rmv_load_motif` runs the local FR3D pipeline fresh for the active PDB and loads Source 5 motifs.
4. `rmv_summary` prints motif type counts.
5. `rmv_show HAIRPIN LOOP` renders one motif family on structure.

Useful follow-up commands:

```pymol
rmv_summary HL
rmv_show HL 1
rmv_view HL
rmv_save HL cif
```

---

## 5. Internal Pipeline

Source 5 pipeline flow inside RSMViewer:

1. Validate Source 5 runtime readiness.
2. Build internal runtime config (`_build_internal_fr3d_config`).
3. Run `run_pipeline(...)` from `rsmviewer/tools/fr3d_loop_extractor.py`.
4. Generate loop CSV for the target PDB in FR3D annotation directory.
5. Load generated motifs through user-annotation loading path.
6. Expose motifs through standard visualization manager APIs.

Data flow summary:

```text
rmv_load_motif (source=5)
  -> ensure_fr3d_runtime_ready()
  -> run_fr3d_wrapper(..., force_fresh)
  -> fr3d_loop_extractor.run_pipeline(...)
  -> <PDB>_fr3d_loops.csv
  -> load_user_annotations_action('fr3d', <PDB>)
  -> motif_loader.loaded_motifs
```

Pipeline outputs:

- `<PDB_ID>_fr3d_loops.csv`
- FR3D pairwise/basepair intermediate files produced in the FR3D output directory

---

## 6. Runtime Commands

### 6.1 `rmv_fr3d status`

Purpose:

- show runtime mode, discovered paths, and Source 5 usage hints

Example:

```pymol
rmv_fr3d status
```

Expected output includes:

- runtime dir
- runtime mode
- FR3D tool path
- Python path
- output path

### 6.2 `rmv_fr3d doctor`

Purpose:

- inspect runtime readiness and report missing requirements

Example:

```pymol
rmv_fr3d doctor
```

Expected output includes:

- status line (`LOCAL PIPELINE READY` or `NOT READY`)
- missing items list when not ready

### 6.3 `rmv_fr3d setup`

Purpose:

- run first-time runtime setup/bootstrap

Example:

```pymol
rmv_fr3d setup
```

Expected output:

- setup logs from bootstrap steps
- runtime becomes ready for local pipeline runs

### 6.4 `rmv_fr3d run <PDB_ID>`

Purpose:

- run Source 5 local pipeline for a specific structure

Example:

```pymol
rmv_fr3d run 1S72
```

Expected behavior:

- removes prior `1S72_fr3d_loops.csv` if present
- executes local pipeline fresh
- loads generated motifs into RSMViewer

### 6.5 `rmv_fr3d run_current`

Purpose:

- run Source 5 for the currently loaded structure

Example:

```pymol
rmv_fetch 1S72
rmv_fr3d run_current
```

### 6.6 `rmv_fr3d refresh [PDB_ID]`

Purpose:

- force rerun for specified or current PDB

Example:

```pymol
rmv_fr3d refresh 1S72
```

### 6.7 `rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]`

Purpose:

- advanced path/config override for Source 5 settings

Example:

```pymol
rmv_fr3d config /absolute/path/to/fr3d-python /absolute/path/to/output on
```

---

## 7. Troubleshooting

### Runtime not detected

Symptom:

- `NOT READY` in `rmv_fr3d doctor`

Fix:

1. Run `rmv_fr3d setup`
2. Re-run `rmv_fr3d doctor`
3. If still not ready, verify Python and FR3D paths in doctor output

### Python dependency issues

Symptom:

- doctor reports missing Python/dependency readiness

Fix:

- run setup again
- ensure host Python can create venv and install packages
- verify environment has network access for package install

### FR3D root not found (advanced override mode)

Symptom:

- advanced config points to an invalid FR3D root

Fix:

- run `rmv_fr3d setup` to return to bundled runtime defaults
- if using custom override, set a valid path with `rmv_fr3d config <FR3D_ROOT> ...`

### Permission errors

Symptom:

- setup cannot create runtime directories or write outputs

Fix:

- ensure write permission for:
  - `rsmviewer/tools/fr3d_runtime/`
  - `rsmviewer/database/user_annotations/fr3d/`

### Pipeline execution failure

Symptom:

- `Local FR3D pipeline failed` in command output

Fix:

1. Run `rmv_fr3d doctor`
2. Confirm runtime is ready
3. Retry `rmv_fr3d run <PDB_ID>`
4. Check output directory write permissions

### Cache/output confusion

Symptom:

- expected old FR3D output reused

Explanation:

- Source 5 runs fresh each time and removes previous loop CSV for that PDB before run.

### Download/bootstrap failures

Symptom:

- setup fails while installing dependencies or cloning FR3D

Fix:

- ensure internet connectivity
- ensure `git` is available in PATH
- re-run `rmv_fr3d setup`

---

## 8. Developer Notes

### 8.1 Core implementation points

Primary methods in `rsmviewer/gui.py`:

- `ensure_fr3d_runtime_ready(...)`
- `fr3d_doctor(...)`
- `_build_internal_fr3d_config(...)`
- `run_fr3d_wrapper(...)`
- `_run_fr3d_local_pipeline(...)`

Runtime setup helper:

- `rsmviewer/tools/fr3d_runtime/setup_runtime.py`

Pipeline executor:

- `rsmviewer/tools/fr3d_loop_extractor.py`

### 8.2 Execution flow for contributors

```text
rmv_fr3d / rmv_load_motif (source 5)
  -> runtime doctor/setup check
  -> internal config build
  -> local FR3D pipeline run
  -> CSV output write
  -> user annotation load
  -> motif visualization pipeline
```

### 8.3 Extension points

Common extension areas:

- setup/runtime discovery logic in `setup_runtime.py`
- loop extraction behavior in `fr3d_loop_extractor.py`
- command UX/help output in `gui.py` wrapper command handlers
- Source 5 post-processing before visualization load

### 8.4 Contributor validation checklist

1. `rmv_fr3d doctor`
2. `rmv_fr3d setup`
3. `rmv_db 5; rmv_fetch 1S72; rmv_load_motif`
4. `rmv_summary` and one render command (`rmv_show HL`)
5. verify Source 5 still runs fresh with no fallback path

---

For repository-wide command usage and general workflows, see:

- `README.md`
- `TUTORIAL.md`
- `DEVELOPER.md`
