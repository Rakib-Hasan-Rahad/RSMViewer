# RSMViewer: FR3D and RNAMotifScanX Integration

Step-by-step implementation and user guide for integrating FR3D (Source 5) and RNAMotifScanX/RMSX (Source 7) in RSMViewer v1.

---

## Table of Contents

1. Overview
2. What Was Integrated in RSMViewer
3. Integration Architecture (Code-Level)
4. First-Time User Setup
5. Step-by-Step User Workflows
6. Runtime Commands (Complete)
7. Troubleshooting
8. Quick Verification Checklist

---

## 1. Overview

RSMViewer integrates two annotation tools:

- FR3D as Source 5 (`rmv_db 5`)
- RNAMotifScanX (RMSX) as Source 7 (`rmv_db 7`)

Both tools are integrated into the same RSMViewer command pipeline:

- load structure
- select source
- run/load motifs
- summarize
- visualize
- export

Core user flow for both:

```pymol
rmv_fetch <PDB_ID>
rmv_db <SOURCE_ID>
rmv_load_motif
rmv_summary
rmv_show <MOTIF_TYPE>
```

---

## 2. What Was Integrated in RSMViewer

### FR3D (Source 5)

- Integrated local runtime under `rsmviewer/tools/fr3d_runtime/`
- First-run setup and diagnostics via `rmv_fr3d setup` and `rmv_fr3d doctor`
- Fresh local pipeline execution for each Source-5 run
- No fallback execution path for Source 5
- Output loaded directly into RSMViewer from `database/user_annotations/fr3d/`

### RMSX (Source 7)

- Integrated runtime under `rsmviewer/tools/rmsx_runtime/`
- Platform-aware binary discovery in `bin/<platform>/`
- First-run setup and diagnostics via `rmv_rmsx setup` / `rmv_rmsx_doctor` / `rmv_rmsx doctor`
- Fresh pipeline execution in Source-7 flow (`force_fresh=True`)
- Source-7 run path enforces local-only behavior (`auto_download_cif=False`)
- Output loaded from `database/user_annotations/RNAMotifScanX/`

---

## 3. Integration Architecture (Code-Level)

### Main integration points

- `rsmviewer/gui.py`
  - Source command routing (`rmv_db`, `rmv_load_motif`)
  - FR3D wrapper commands (`rmv_fr3d ...`)
  - RMSX wrapper commands (`rmv_rmsx ...`, `rmv_rmsx_doctor`)
  - Runtime doctor/setup hooks

- `rsmviewer/tools/fr3d_runtime/setup_runtime.py`
  - FR3D runtime detection and bootstrap

- `rsmviewer/tools/fr3d_loop_extractor.py`
  - FR3D loop extraction pipeline runner

- `rsmviewer/tools/rmsx_runtime/setup_runtime.py`
  - RMSX runtime detection and setup/build helper

### Source 5 execution flow

```text
rmv_db 5
  -> select Source 5

rmv_load_motif
  -> load_user_annotations_action('fr3d', pdb)
  -> run_fr3d_wrapper(pdb, force_refresh=True)
  -> ensure_fr3d_runtime_ready(auto_setup=True)
  -> _run_fr3d_local_pipeline(...)
  -> write <PDB>_fr3d_loops.csv
  -> load FR3D motifs into motif_loader
```

### Source 7 execution flow

```text
rmv_db 7
  -> select Source 7

rmv_load_motif
  -> load_user_annotations_action('rnamotifscanx', pdb)
  -> ensure_rmsx_runtime_ready(auto_setup=True)
  -> rmsx_runner.run_pipeline(..., force_fresh=True)
     with auto_download_cif=False
  -> load RMSX motifs into motif_loader
```

---

## 4. First-Time User Setup

## 4.1 FR3D first-time setup

In PyMOL:

```pymol
rmv_fr3d doctor
rmv_fr3d setup
rmv_fr3d status
```

Expected result:

- doctor/status report runtime ready
- Source 5 can run with `rmv_db 5` + `rmv_load_motif`

## 4.2 RMSX first-time setup

In PyMOL:

```pymol
rmv_rmsx_doctor
rmv_rmsx setup
rmv_rmsx status
```

Expected result:

- runtime report shows discovered executable paths and readiness
- Source 7 can run through `rmv_load_motif` or `rmv_rmsx run_current`

---

## 5. Step-by-Step User Workflows

## 5.1 FR3D workflow (Source 5)

Step 1: load structure

```pymol
rmv_fetch 1S72
```

Step 2: select FR3D source

```pymol
rmv_db 5
```

Step 3: load motifs (runs local FR3D pipeline)

```pymol
rmv_load_motif
```

Step 4: inspect and render

```pymol
rmv_summary
rmv_summary HL
rmv_show HL
rmv_show HAIRPIN LOOP
```

Step 5 (optional): force rerun

```pymol
rmv_fr3d refresh
```

## 5.2 RMSX workflow (Source 7)

Step 1: load structure

```pymol
rmv_fetch 1S72
```

Step 2: select RMSX source

```pymol
rmv_db 7
```

Step 3: load motifs (runs integrated Source-7 runtime)

```pymol
rmv_load_motif
```

Step 4: inspect and render

```pymol
rmv_summary
rmv_show K-TURN
rmv_show SARCIN-RICIN
```

Step 5 (optional): explicit runtime run command

```pymol
rmv_rmsx run_current
```

## 5.3 Combined FR3D + RMSX workflow

```pymol
rmv_fetch 1S72
rmv_db 5
rmv_load_motif
rmv_db 7
rmv_db 5 7
rmv_load_motif
rmv_summary K-TURN
rmv_show K-TURN rmsx
rmv_show K-TURN shared
```

---

## 6. Runtime Commands (Complete)

### FR3D commands

- `rmv_fr3d status`
- `rmv_fr3d doctor`
- `rmv_fr3d setup`
- `rmv_fr3d refresh [PDB_ID]`
- `rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]` (advanced override)
- `rmv_fr3d run <PDB_ID>`
- `rmv_fr3d run_current`

Notes:

- `run`/`run_current` execute local Source-5 pipeline and load results
- `refresh` forces rerun for the active/specified PDB

### RMSX commands

- `rmv_rmsx status`
- `rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH]`
- `rmv_rmsx args <ARG_TEMPLATE>`
- `rmv_rmsx doctor`
- `rmv_rmsx setup`
- `rmv_rmsx test`
- `rmv_rmsx run <PDB_ID> [EXTRA_ARGS]`
- `rmv_rmsx run_current [EXTRA_ARGS]`
- `rmv_rmsx_doctor`

Notes:

- `rmv_rmsx_doctor` prints integrated runtime diagnostics directly
- `rmv_rmsx run*` uses integrated pipeline config path when available

---

## 7. Troubleshooting

### FR3D issues

`FR3D runtime is incomplete`:

```pymol
rmv_fr3d doctor
rmv_fr3d setup
rmv_fr3d status
```

`No motifs after Source 5 run`:

- ensure a structure is loaded (`rmv_fetch`)
- ensure Source 5 selected (`rmv_db 5`)
- rerun `rmv_load_motif`
- force rerun with `rmv_fr3d refresh`

### RMSX issues

`RMSX runtime is incomplete`:

```pymol
rmv_rmsx_doctor
rmv_rmsx setup
rmv_rmsx status
```

`No RMSX results loaded`:

- ensure structure is loaded (`rmv_fetch`)
- ensure Source 7 selected (`rmv_db 7`)
- run `rmv_load_motif` or `rmv_rmsx run_current`
- if macOS/Windows and binaries are unavailable, provide compatible runtime assets as reported by doctor/setup

### General checks

```pymol
rmv_sources
rmv_source info
rmv_summary
rmv_chains
```

---

## 8. Quick Verification Checklist

FR3D verification:

```pymol
rmv_fr3d setup
rmv_db 5
rmv_fetch 1S72
rmv_load_motif
rmv_summary
```

RMSX verification:

```pymol
rmv_rmsx setup
rmv_db 7
rmv_fetch 1S72
rmv_load_motif
rmv_summary
```

Combined verification:

```pymol
rmv_db 5 7
rmv_load_motif
rmv_summary K-TURN
rmv_show K-TURN shared
```

---

For FR3D runtime internals and deep troubleshooting, see `FR3D_INTEGRATION.md`.
