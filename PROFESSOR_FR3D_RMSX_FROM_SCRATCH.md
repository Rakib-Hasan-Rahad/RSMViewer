# Professor Guide: How FR3D and RNAMotifScanX Were Integrated into RSMViewer (From Scratch)

This document is a beginner-friendly, full technical explanation of how we integrated FR3D and RNAMotifScanX (RMSX) into RSMViewer v1.

It explains:

- what problem we solved,
- which files we implemented,
- how commands trigger each pipeline,
- how data moves from structure input to motif visualization,
- how we handled cross-platform behavior,
- and how to troubleshoot each stage.

---

## 1. Big Picture (Beginner Overview)

RSMViewer is a PyMOL-based RNA motif analysis plugin.

We integrated two external motif pipelines as user-controlled sources:

- Source 5: FR3D (loop extraction)
- Source 7: RNAMotifScanX / RMSX (motif-family search)

Goal of the integration:

1. User loads a structure (`rmv_fetch 1S72`)
2. User selects source (`rmv_db 5` or `rmv_db 7`)
3. User runs loading (`rmv_load_motif`)
4. RSMViewer runs the correct pipeline automatically
5. Results are converted to motif instances and rendered with `rmv_summary`, `rmv_show`, `rmv_view`

We designed both Source 5 and Source 7 around strict fresh execution:

- fresh pipeline run for current PDB
- no stale cached outputs reused in source run path
- source behavior stays predictable for teaching/demo/research

---

## 2. Repository Map (Where Integration Lives)

## 2.1 Main orchestration

- `rsmviewer/gui.py`

This is the main integration controller. It contains:

- source selection and command routing,
- runtime setup/doctor calls,
- FR3D and RMSX wrappers,
- loading results into RSMViewer motif loader.

## 2.2 FR3D runtime and pipeline files

- `rsmviewer/tools/fr3d_runtime/setup_runtime.py`
- `rsmviewer/tools/fr3d_loop_extractor.py`

Responsibilities:

- detect/bootstrap FR3D runtime (python + deps + FR3D source)
- run FR3D pairwise annotation
- extract HL/IL loops with flankSS logic
- write `<PDB>_fr3d_loops.csv`

## 2.3 RMSX runtime and pipeline files

- `rsmviewer/tools/rmsx_runtime/setup_runtime.py`
- `rsmviewer/tools/rmsx_runner.py`

Responsibilities:

- detect platform and binaries
- validate queries and matrix assets
- optional source build path for runtime tools
- run MC-Annotate + RNAMotifScanX family search
- write `result_0_100_withbs.log` per family

## 2.4 Output directories used by RSMViewer

- FR3D output: `rsmviewer/database/user_annotations/fr3d/`
- RMSX output: `rsmviewer/database/user_annotations/RNAMotifScanX/`

These outputs are consumed by the user annotation provider and loaded into motifs.

## 2.5 Raw files and resources we used

This section is the exact "inputs/resources" view for both integrations.

### FR3D raw resources

Bundled/runtime resources:

- `rsmviewer/tools/fr3d_runtime/setup_runtime.py`
- `rsmviewer/tools/fr3d_runtime/python/fr3d-python/` (FR3D source checkout)
- `rsmviewer/tools/fr3d_runtime/python/venv/` (runtime virtual environment)

Raw structure and interaction artifacts used/generated during runs:

- `rsmviewer/database/user_annotations/fr3d/1s72.cif.gz`
- `rsmviewer/database/user_annotations/fr3d/4v9f.cif.gz`
- `rsmviewer/database/user_annotations/fr3d/1S72_basepair.txt`
- `rsmviewer/database/user_annotations/fr3d/4V9F_basepair.txt`
- `rsmviewer/database/user_annotations/fr3d/5BTP_basepair.txt`
- `rsmviewer/database/user_annotations/fr3d/8GLP_basepair.txt`
- `rsmviewer/database/user_annotations/fr3d/1S72_fr3d_loops.csv`
- `rsmviewer/database/user_annotations/fr3d/4V9F_fr3d_loops.csv`
- `rsmviewer/database/user_annotations/fr3d/5BTP_fr3d_loops.csv`
- `rsmviewer/database/user_annotations/fr3d/8GLP_fr3d_loops.csv`

Persistent FR3D wrapper config:

- `rsmviewer/database/user_annotations/fr3d/.fr3d_wrapper.json`

### RMSX raw resources

Bundled/runtime resources:

- `rsmviewer/tools/rmsx_runtime/setup_runtime.py`
- `rsmviewer/tools/rmsx_runtime/queries/k-turn_consensus.struct`
- `rsmviewer/tools/rmsx_runtime/queries/c-loop_consensus.struct`
- `rsmviewer/tools/rmsx_runtime/queries/sarcin-ricin_consensus.struct`
- `rsmviewer/tools/rmsx_runtime/queries/reverse-kturn_consensus.struct`
- `rsmviewer/tools/rmsx_runtime/queries/e-loop_consensus.struct`
- `rsmviewer/tools/rmsx_runtime/mat/nuc.mat`
- `rsmviewer/tools/rmsx_runtime/mat/iso.mat`
- `rsmviewer/tools/rmsx_runtime/mat/stk.mat`
- `rsmviewer/tools/rmsx_runtime/src/RNAMotifScanX_src/` (RMSX source tree)
- `rsmviewer/tools/rmsx_runtime/src/external/` (external build workspace)

Raw structure and intermediate artifacts used/generated during runs:

- `rsmviewer/database/user_annotations/RNAMotifScanX/1S72.pdb`
- `rsmviewer/database/user_annotations/RNAMotifScanX/1S72_mc_annotate.out`
- `rsmviewer/database/user_annotations/RNAMotifScanX/1S72_0.rmsx.in`
- `rsmviewer/database/user_annotations/RNAMotifScanX/k-turn_consensus/result_0_100_withbs.log`
- `rsmviewer/database/user_annotations/RNAMotifScanX/c-loop_consensus/result_0_100_withbs.log`
- `rsmviewer/database/user_annotations/RNAMotifScanX/sarcin-ricin_consensus/result_0_100_withbs.log`
- `rsmviewer/database/user_annotations/RNAMotifScanX/reverse-kturn_consensus/result_0_100_withbs.log`
- `rsmviewer/database/user_annotations/RNAMotifScanX/e-loop_consensus/result_0_100_withbs.log`

Data-flow meaning (important for explanation):

1. Structure file (CIF/PDB) is the physical input.
2. Tool-specific intermediate files are generated:
   - FR3D: `_basepair.txt`
   - RMSX: `_mc_annotate.out` and `.rmsx.in`
3. Final annotation file format consumed by RSMViewer is produced:
   - FR3D: `_fr3d_loops.csv`
   - RMSX: `result_0_100_withbs.log`
4. `UserAnnotationProvider` parses these final files into motif instances.

---

## 3. Step-by-Step: Integration from Scratch

## Step 1: Add runtime state into GUI controller

In `rsmviewer/gui.py` initialization:

- FR3D config/runtime state is created:
  - `fr3d_config_file`, `fr3d_output_path`, `fr3d_runtime_dir`, etc.
- RMSX config/runtime state is created:
  - `rmsx_config_file`, `rmsx_output_path`, `rmsx_runtime_dir`, etc.

Why this matters:

- all pipeline settings are centralized in one controller object
- wrappers can persist paths and startup behavior

## Step 2: Implement runtime doctor/setup hooks

For FR3D:

- runtime script: `fr3d_runtime/setup_runtime.py`
- GUI calls it through runtime helper methods
- commands: `rmv_fr3d doctor`, `rmv_fr3d setup`, `rmv_fr3d status`

For RMSX:

- runtime script: `rmsx_runtime/setup_runtime.py`
- GUI method `_run_rmsx_runtime_setup(...)` invokes script with JSON output
- commands: `rmv_rmsx doctor`, `rmv_rmsx setup`, `rmv_rmsx status`, `rmv_rmsx_doctor`

Why this matters:

- users get one-command diagnostics
- first run can bootstrap runtime assets automatically

## Step 3: Connect source IDs to strict pipeline behavior

In source loading flow (`load_user_annotations_action(...)` in `gui.py`):

- Source 5 path (`tool == 'fr3d'`):
  - runs FR3D pipeline fresh using `run_fr3d_wrapper(..., force_refresh=True)`
  - returns immediately after wrapper load

- Source 7 path (`tool in ['rmsx', 'rnamotifscanx']`):
  - ensures runtime is ready
  - builds or uses integrated runtime config
  - forces fresh run using `rmsx_runner.run_pipeline(..., force_fresh=True)`
  - forces `auto_download_cif=False` in source-7 run path

Why this matters:

- `rmv_load_motif` becomes the single user action for both tools
- source behavior is deterministic and consistent

## Step 4: Implement explicit FR3D wrapper command API

In `gui.py`, FR3D command supports:

- `rmv_fr3d status`
- `rmv_fr3d doctor`
- `rmv_fr3d setup`
- `rmv_fr3d refresh [PDB_ID]`
- `rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]`
- `rmv_fr3d run <PDB_ID>`
- `rmv_fr3d run_current`

Internally, `run_fr3d_wrapper(...)` does:

1. validates PDB and output dir
2. removes previous `<PDB>_fr3d_loops.csv` for fresh execution
3. ensures runtime ready
4. builds local pipeline config (`force_fresh=True`)
5. calls `_run_fr3d_local_pipeline(...)`
6. points source 5 to output dir
7. loads motifs via user annotations provider

Why this matters:

- allows both auto and manual control
- guarantees Source 5 uses local pipeline path in current implementation

## Step 5: Implement explicit RMSX wrapper command API

In `gui.py`, RMSX command supports:

- `rmv_rmsx status`
- `rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH]`
- `rmv_rmsx args <ARG_TEMPLATE>`
- `rmv_rmsx doctor`
- `rmv_rmsx setup`
- `rmv_rmsx test`
- `rmv_rmsx run <PDB_ID> [EXTRA_ARGS]`
- `rmv_rmsx run_current [EXTRA_ARGS]`
- `rmv_rmsx_doctor`

Internally, `run_rmsx_wrapper(...)` does:

1. validates PDB
2. ensures runtime ready
3. builds integrated config from runtime assets
4. enforces `auto_download_cif=False` in strict source path
5. runs `rmsx_runner.run_pipeline(..., force_fresh=True)`
6. sets source 7 output path
7. loads motifs into RSMViewer and reports family counts

Special note from implementation:

- when integrated pipeline mode is active, extra CLI args are ignored (explicitly logged)

## Step 6: Register command UX and source flow

User commands used in teaching/demo sequence:

```pymol
rmv_fetch 1S72
rmv_db 5
rmv_load_motif
rmv_summary
rmv_show HL
```

and for RMSX:

```pymol
rmv_fetch 1S72
rmv_db 7
rmv_load_motif
rmv_summary
rmv_show K-TURN
```

## Step 7: What we modified to integrate (exact engineering changes)

This is the specific change list you can present as implementation work.

### A. Main controller changes in `rsmviewer/gui.py`

We added/updated integration behavior in these areas:

1. Runtime state and persistent config fields
   - FR3D: output/runtime/config paths
   - RMSX: executable/output/runtime/config/args fields

2. Runtime readiness methods
   - FR3D readiness + setup attempt path
   - RMSX readiness via `_run_rmsx_runtime_setup(...)`

3. Source 5 strict execution path
   - `run_fr3d_wrapper(...)` removes previous loop CSV
   - forces fresh local pipeline execution
   - disables fallback behavior in source path

4. Source 7 strict execution path
   - integrated config build from runtime assets (`_build_internal_rmsx_config`)
   - `rmsx_runner.run_pipeline(..., force_fresh=True)`
   - enforced `auto_download_cif=False` in source pipeline path

5. `rmv_load_motif` source-trigger integration (`load_user_annotations_action(...)`)
   - source 5 auto-runs FR3D fresh pipeline
   - source 7 auto-runs RMSX fresh pipeline before loading

6. Command surface exposed to users
   - FR3D wrapper command family (`rmv_fr3d ...`)
   - RMSX wrapper command family (`rmv_rmsx ...`, `rmv_rmsx_doctor`)

### B. FR3D pipeline adaptation in `rsmviewer/tools/fr3d_loop_extractor.py`

What we implemented/adapted for integration:

1. FR3D pairwise output parsing into internal structures.
2. Canonical cWW filtering logic for stem-defining basepairs.
3. Crossing-number handling to exclude pseudoknot crossings for loop-boundary definition.
4. Secondary-structure graph and flankSS-style loop boundary detection.
5. HL/IL extraction and export in RSMViewer-consumable CSV format.
6. A callable `run_pipeline(...)` API so GUI can run it programmatically.

### C. FR3D runtime bootstrap in `rsmviewer/tools/fr3d_runtime/setup_runtime.py`

What we implemented/adapted:

1. FR3D root discovery using environment variables and candidate paths.
2. Python interpreter discovery with dependency validation (`scipy`, `pdbx`).
3. Optional bootstrap flow:
   - create venv,
   - install runtime deps,
   - clone `fr3d-python`,
   - install FR3D package editable.
4. JSON doctor/setup report consumed by GUI.

### D. RMSX runtime adaptation in `rsmviewer/tools/rmsx_runtime/setup_runtime.py`

What we implemented/adapted:

1. Platform detection (macOS arm64/x86_64, Linux x86_64, Windows x86_64).
2. Runtime binary validation and compatibility checks.
3. Query/matrix asset validation (`queries/*.struct`, `mat/*.mat`).
4. Optional source build helpers for MCCORE and MC-Annotate.
5. Build compatibility patches and packaging support logic.
6. JSON doctor/setup report consumed by GUI and CLI commands.

### E. RMSX pipeline adaptation in `rsmviewer/tools/rmsx_runner.py`

What we implemented/adapted:

1. High-level `run_pipeline(...)` callable API for GUI integration.
2. Result existence checks per motif family.
3. MC-Annotate preprocessing and annotation file generation.
4. RMSX target preparation (`.rmsx.in`) from annotation output.
5. Per-family query execution and output log generation.
6. Modern and legacy RMSX CLI fallback invocation paths.
7. Runtime environment setup for matrix path resolution.

### F. Annotation ingestion integration

In the user annotation loading path, we aligned generated files with existing motif parsing so that:

- FR3D generated CSV is directly loaded as Source 5 motifs.
- RMSX family logs are directly loaded as Source 7 motifs.

This let us integrate new runtimes without replacing the entire visualization stack.

---

## 4. How FR3D Pipeline Works Internally (Detailed)

FR3D implementation file: `rsmviewer/tools/fr3d_loop_extractor.py`

Pipeline stages:

1. Input acquisition
   - pipeline receives PDB ID and config
   - CIF/PDB is located (and may be downloaded if enabled by standalone config path)

2. Pairwise annotation
   - executes FR3D pairwise interaction annotator (`NA_pairwise_interactions.py`)

3. Basepair filtering
   - parser reads `_basepairs.txt`
   - keeps canonical cWW stem-defining interactions
   - excludes near interactions and non-nested crossing cases

4. Secondary-structure graph assembly
   - chain-aware mapping of sequence positions and partner pairs

5. Loop boundary extraction (flankSS logic)
   - Hairpin loops (HL)
   - Internal loops (IL)

6. Output formatting
   - writes `<PDB>_fr3d_loops.csv` in motif-loader compatible format

7. GUI loading
   - source 5 path points to FR3D output dir
   - motif loader imports and indexes motif types/instances

FR3D input/output chain you can present:

1. Raw structure: `*.cif` or `*.cif.gz`
2. Raw FR3D interaction output: `<PDB>_basepair.txt`
3. Processed motif output: `<PDB>_fr3d_loops.csv`
4. In-memory motifs in RSMViewer loader

Practical example:

```pymol
rmv_fr3d setup
rmv_fetch 1S72
rmv_db 5
rmv_load_motif
rmv_summary HL
rmv_show HAIRPIN LOOP
```

---

## 5. How RMSX Pipeline Works Internally (Detailed)

RMSX pipeline file: `rsmviewer/tools/rmsx_runner.py`

Pipeline stages:

1. Resolve runtime config
   - executable paths
   - query motifs directory
   - output directory
   - family list (`k-turn`, `c-loop`, `sarcin-ricin`, `reverse-kturn`, `e-loop`)

2. Fresh-run check
   - if `force_fresh=True`, old result logs are ignored/removed as needed

3. Structure input stage
   - locate CIF/PDB from configured paths
   - in Source-7 strict path, auto external CIF download is disabled

4. Annotation stage
   - run MC-Annotate on PDB-formatted input
   - write `<PDB>_mc_annotate.out`

5. RMSX input preparation
   - parse annotation output
   - generate chain target file(s) (`<PDB>_<chain>.rmsx.in`)
   - align sequence references when available

6. Family scan stage
   - for each motif family:
     - find consensus query file
     - run RNAMotifScanX
     - save `result_0_100_withbs.log` in family subfolder

7. GUI loading stage
   - source 7 path is updated to runtime output dir
   - motifs are loaded and rendered by existing show/view commands

RMSX input/output chain you can present:

1. Raw structure: `<PDB>.pdb` (and CIF discovery path)
2. Raw annotation: `<PDB>_mc_annotate.out`
3. Prepared RMSX target: `<PDB>_<chain>.rmsx.in`
4. Processed family outputs: `<family>_consensus/result_0_100_withbs.log`
5. In-memory motifs in RSMViewer loader

Practical example:

```pymol
rmv_rmsx setup
rmv_fetch 1S72
rmv_db 7
rmv_load_motif
rmv_summary
rmv_show K-TURN
```

---

## 6. Cross-Platform Strategy (How We Made It Portable)

## 6.1 Platform detection

In `rsmviewer/tools/rmsx_runtime/setup_runtime.py`:

- detects OS + CPU and maps to runtime platform folder
  - `macos-arm64`
  - `macos-x86_64`
  - `windows-x86_64`
  - `linux-x86_64`

In `gui.py` (`_build_internal_rmsx_config`):

- chooses matching `bin/<platform>/` first
- if missing, scans other runtime bin folders as fallback discovery

## 6.2 Binary compatibility checks

RMSX runtime setup validates binary format (ELF/PE/Mach-O) and expected tools:

- RNAMotifScanX executable
- MC-Annotate executable

If missing or incompatible:

- doctor/setup reports clear missing items
- user can provide runtime binaries in platform bin folder

## 6.3 Source-build path for difficult platforms

`setup_runtime.py` includes helper builders:

- MCCORE build logic (with compatibility patches)
- MC-Annotate build logic using CMake
- optional packaging patching for macOS runtime libs

Why this matters:

- integration can still proceed where prebuilt binaries are not ready
- runtime script centralizes dependency diagnostics

## 6.4 macOS-specific practical policy for RMSX execution

In `run_rmsx_wrapper(...)` (`gui.py`):

- when no results are produced on macOS, message explains Linux x86-64 execution limitation
- recommended workflow:
  1. generate RMSX logs on Linux
  2. copy `result_0_100_withbs.log` files into output directory
  3. load in RSMViewer on macOS

## 6.5 FR3D portability

FR3D runtime setup uses Python + pip dependency checks:

- detects Python interpreters
- checks `scipy` and `pdbx` availability
- can bootstrap local virtual environment and clone/install FR3D source

This design is naturally cross-platform because it relies on Python runtime management.

---

## 7. End-to-End Demonstration Scripts

## 7.1 FR3D end-to-end demo

```pymol
rmv_fr3d doctor
rmv_fr3d setup
rmv_fetch 1S72
rmv_db 5
rmv_load_motif
rmv_summary
rmv_summary HL
rmv_show HL
rmv_view HL 1
```

## 7.2 RMSX end-to-end demo

```pymol
rmv_rmsx_doctor
rmv_rmsx setup
rmv_fetch 1S72
rmv_db 7
rmv_load_motif
rmv_summary
rmv_show K-TURN
rmv_view K-TURN 1
```

## 7.3 Combined-source demo

```pymol
rmv_fetch 1S72
rmv_db 5 7
rmv_load_motif
rmv_summary K-TURN
rmv_show K-TURN shared
rmv_show K-TURN rmsx
```

---

## 8. Troubleshooting Playbook (Stage-by-Stage)

## Stage A: runtime not ready

FR3D:

```pymol
rmv_fr3d doctor
rmv_fr3d setup
rmv_fr3d status
```

RMSX:

```pymol
rmv_rmsx_doctor
rmv_rmsx setup
rmv_rmsx status
```

What to check:

- missing binary/dependency names in doctor output
- runtime path printed by doctor
- platform folder under runtime `bin/`

## Stage B: pipeline ran but no motifs loaded

Checklist:

1. structure is loaded (`rmv_fetch <PDB>`)
2. source is selected (`rmv_db 5` or `rmv_db 7`)
3. rerun `rmv_load_motif`
4. inspect summary with `rmv_summary`

For FR3D:

- verify `<PDB>_fr3d_loops.csv` exists in FR3D output dir

For RMSX:

- verify family logs exist:
  - `<output>/<family>_consensus/result_0_100_withbs.log`

## Stage C: RMSX executable problems

Use:

```pymol
rmv_rmsx test
```

If fail:

- fix executable path with `rmv_rmsx config ...`
- ensure executable permissions
- validate working directory exists

## Stage D: source mismatch or stale expectation

Use quick status commands:

```pymol
rmv_sources
rmv_source info
rmv_summary
rmv_chains
```

---

## 9. Why This Integration Design Is Strong

1. One command flow for users (`rmv_db` + `rmv_load_motif`)
2. Dedicated wrappers for explicit control (`rmv_fr3d`, `rmv_rmsx`)
3. Runtime doctor/setup built into plugin
4. Clear source-specific behavior
5. Cross-platform detection and diagnostics
6. Fresh run policies for reproducibility in pipeline-centric workflows

---

## 10. Quick Viva/Professor Explanation Script

If your professor asks "how exactly did you integrate these tools?", explain it in this order:

1. We added FR3D and RMSX as dedicated user sources (5 and 7).
2. We added runtime setup/doctor scripts for each tool.
3. We connected source selection to automatic pipeline execution in `rmv_load_motif`.
4. We added wrapper command groups for manual execution and diagnostics.
5. We routed pipeline outputs into the existing user-annotation motif loader.
6. We added platform detection and binary/runtime validation for portability.
7. We ensured strict fresh-run behavior in source pipeline paths for consistency.

---

## 11. Related Documentation

- `FR3D_INTEGRATION.md` (FR3D deep runtime details)
- `FR3D_RMSX_INTEGRATION.md` (concise combined integration guide)
- `DEVELOPER.md` (plugin architecture and command map)
- `README.md` and `TUTORIAL.md` (user-facing workflows)
