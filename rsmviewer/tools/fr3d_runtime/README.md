# FR3D Runtime Integration in RSMViewer

This directory contains the integrated FR3D runtime support used by RSMViewer Source 5.

FR3D (Find RNA 3D) is a geometry-based RNA structural annotation toolkit. RSMViewer uses it to generate local loop and motif annotations directly from the currently loaded structure, then loads those annotations into the standard RSMViewer visualization workflow.

---

## Runtime Role

RSMViewer exposes FR3D as Source 5:

```pymol
rmv_db 5
```

When Source 5 is active, `rmv_load_motif` runs the bundled local FR3D pipeline and loads the generated loop annotations into RSMViewer. The user can then inspect, render, superimpose, and export motifs using the same commands used for other RSMViewer sources.

Typical workflow:

```pymol
rmv_fetch 1S72
rmv_db 5
rmv_load_motif
rmv_summary
rmv_show HL
```

---

## Directory Layout

Important files and folders:

- `setup_runtime.py`: FR3D runtime setup and diagnostic helper.
- `python/fr3d-python/`: vendored FR3D source tree used by Source 5.
- `python/venv/`: runtime virtual environment created locally by setup; not tracked in Git.
- `../fr3d_loop_extractor.py`: RSMViewer pipeline runner that calls FR3D and converts output into RSMViewer loop CSV format.
- `../../database/user_annotations/fr3d/`: generated Source-5 annotation output directory.

The vendored FR3D source tree is tracked as normal files in this repository. It is not a Git submodule and it does not contain nested `.git` metadata.

---

## Runtime Setup

Run these commands inside PyMOL after installing RSMViewer:

```pymol
rmv_fr3d doctor
rmv_fr3d setup
rmv_fr3d status
```

Expected result:

- doctor/status reports the runtime as ready
- `python/fr3d-python/` is detected
- a runtime Python environment can import the required FR3D dependencies
- Source 5 can run with `rmv_db 5` and `rmv_load_motif`

`rmv_fr3d setup` may create a local virtual environment and install Python dependencies. The virtual environment is intentionally ignored by Git because it is machine-specific.

---

## Source-5 Execution Policy

Source 5 uses strict local execution:

- runs the bundled local FR3D pipeline
- does not use cached FR3D loop CSV files as a replacement for execution
- removes prior output for the requested PDB before rerunning
- does not fall back to BGSU download behavior in Source-5 flow
- keeps generated annotations under `rsmviewer/database/user_annotations/fr3d/`

This policy makes reviewer and user behavior reproducible: Source 5 means local FR3D execution, not a mixed local/remote/cache path.

---

## Pipeline Flow

High-level flow:

```text
rmv_db 5
  -> select Source 5

rmv_load_motif
  -> load_user_annotations_action('fr3d', pdb)
  -> run_fr3d_wrapper(pdb, force_refresh=True)
  -> ensure_fr3d_runtime_ready(auto_setup=True)
  -> build internal FR3D config
  -> run fr3d_loop_extractor.run_pipeline(...)
  -> write <PDB_ID>_fr3d_loops.csv
  -> load generated motifs into RSMViewer
```

The generated CSV is then available to standard commands:

```pymol
rmv_summary
rmv_summary HL
rmv_show HL
rmv_view HL
rmv_save HL cif
```

---

## Commands

### `rmv_fr3d doctor`

Diagnoses runtime readiness and reports missing requirements.

```pymol
rmv_fr3d doctor
```

Use this first when Source 5 does not run.

### `rmv_fr3d setup`

Performs first-time runtime setup.

```pymol
rmv_fr3d setup
```

This prepares the local runtime environment around the bundled FR3D source tree.

### `rmv_fr3d status`

Prints runtime status, output path, auto-run state, and runtime discovery details.

```pymol
rmv_fr3d status
```

### `rmv_fr3d run <PDB_ID>`

Runs the Source-5 local FR3D pipeline for a specific structure.

```pymol
rmv_fr3d run 1S72
```

### `rmv_fr3d run_current`

Runs Source 5 for the currently loaded structure.

```pymol
rmv_fetch 1S72
rmv_fr3d run_current
```

### `rmv_fr3d refresh [PDB_ID]`

Forces a fresh rerun for the specified or active PDB.

```pymol
rmv_fr3d refresh 1S72
```

### `rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]`

Advanced override for custom FR3D runtime paths.

```pymol
rmv_fr3d config /absolute/path/to/fr3d-python /absolute/path/to/output on
```

Normal users should not need this command because the bundled runtime is used by default.

---

## Troubleshooting

### Runtime is incomplete

Run:

```pymol
rmv_fr3d doctor
rmv_fr3d setup
rmv_fr3d status
```

If the runtime still reports incomplete, check that Python can create a virtual environment and install packages.

### No motifs after Source 5 run

Check:

- a structure is loaded with `rmv_fetch <PDB_ID>`
- Source 5 is selected with `rmv_db 5`
- `rmv_fr3d doctor` reports runtime readiness
- rerun with `rmv_fr3d refresh <PDB_ID>`

### Permission errors

The plugin must be able to write to:

- `rsmviewer/tools/fr3d_runtime/`
- `rsmviewer/database/user_annotations/fr3d/`

### Advanced override path is invalid

If a custom FR3D path was configured incorrectly, run setup again to return to bundled runtime behavior:

```pymol
rmv_fr3d setup
rmv_fr3d status
```

Or provide a valid override path:

```pymol
rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]
```

---

## Developer Notes

Primary implementation points:

- `rsmviewer/gui.py`
  - `ensure_fr3d_runtime_ready(...)`
  - `fr3d_doctor(...)`
  - `_build_internal_fr3d_config(...)`
  - `run_fr3d_wrapper(...)`
  - `_run_fr3d_local_pipeline(...)`
- `rsmviewer/tools/fr3d_runtime/setup_runtime.py`
- `rsmviewer/tools/fr3d_loop_extractor.py`

Release packaging notes:

- Keep `python/fr3d-python/` tracked as regular files.
- Do not commit `python/venv/`.
- Do not commit generated CIF, basepair, motif, or loop output files.
- Do not reintroduce nested `.git` metadata under the vendored FR3D tree.
