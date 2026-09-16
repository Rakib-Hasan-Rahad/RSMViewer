# RSMViewer configuration (`config/`)

This folder holds the two JSON configuration files that control the external
motif-search pipelines:

| File | Controls | Used by |
| --- | --- | --- |
| `rmsx_config.json` | RNAMotifScanX (RMSX) source | `rmv_db RNAMotifScanX` / source 7 |
| `fr3d_config.json` | FR3D source | `rmv_db FR3D` / source 5 |

The RNA 3D Motif Atlas and Rfam sources use their local cache files when
available and only contact their public APIs when a cache entry is missing or
refresh is explicitly requested.

### How paths are resolved

All relative paths in these files are resolved **relative to this `config/`
folder**. For example, `../external/rmsx/bin/scan` points to
`<project>/external/rmsx/bin/scan`. You may also use absolute paths.

### How to load a config

The default files are picked up automatically. FR3D can be registered from a
custom config at runtime:

```text
rmv_fr3d register /absolute/path/to/fr3d_config.json     # FR3D
```

Edit a value, save the file, then re-run the `rmv_db` command (or `rmv_refresh`)
to apply it.

---

## `rmsx_config.json` — RNAMotifScanX

RMSX has two ways to provide annotations, selected by `data_mode`.

### Fields

| Field | Type | Meaning / possible values |
| --- | --- | --- |
| `data_mode` | string | `"preannotated"` (default) reads prebuilt results; `"run_from_scratch"` executes the RMSX binaries. |
| `rmsx_executable` | path | RNAMotifScanX `scan` binary. Only needed for `run_from_scratch`. |
| `mc_annotate_executable` | path | `MC-Annotate` binary (base-pair annotation). `run_from_scratch` only. |
| `rnaview_executable` | path | Optional `rnaview` binary; leave as-is to skip the RNAVIEW merge. |
| `rnaview_dir` | path | Directory with RNAVIEW `BASEPARS` resources (only if `rnaview_executable` is set). |
| `pdb_prebuild_archive` | path | `.tar.gz` of preannotated results (default input for `preannotated`). |
| `pdb_prebuild_dir` | path | An already-extracted preannotated directory; used in preference to the archive when it exists. |
| `output_dir` | path | Where per-family result logs are written/read (`../output/rmsx_results`). |
| `motif_families` | list | Families to load: `k-turn`, `c-loop`, `sarcin-ricin`, `reverse-kturn`, `e-loop`. |
| `pvalue_thresholds` | object | Per-family P-value cutoff. Lower = stricter (fewer, higher-confidence hits). Omit a family to keep its paper default. |

Note: the residue-set redundancy threshold (Jaccard 0.60) is **not** a config
value; it is defined in code (`rsmviewer/database/consolidated_table.py`,
`DEFAULT_JACCARD_THRESHOLD`).

### Preannotated mode (default, no binaries needed)

1. Obtain the preannotated bundle (`rmsx_preannotated_input_output.tar.gz`) from Figshare:
   **[https://doi.org/10.6084/m9.figshare.33826795](https://doi.org/10.6084/m9.figshare.33826795)**.
2. Place it at the path in `pdb_prebuild_archive`, i.e.
   `external/rmsx_preannotated/rmsx_preannotated_input_output.tar.gz`.
   Alternatively, extract it and point `pdb_prebuild_dir` at the resulting
   `rmsx_work_default/` folder.
3. Run:

   ```text
   rmv_fetch 1S72
   rmv_db RNA3DMotifAtlas,RNAMotifScanX
   rmv_select SR, 1S72, RNA3DMotifAtlas and RNAMotifScanX, as group_TP
   ```

RSMViewer copies the matching family logs for the requested PDB into
`output_dir` and loads them. The Jaccard/containment consolidation then aligns
them with the other sources.

### Preannotated data layout

Inside the archive (or `pdb_prebuild_dir`), data is keyed by PDB id and chain:

```text
rmsx_work_default/
└── <pdb_id_lowercase>/                 e.g. 1s72/
    ├── _prep_main/                      inputs used to build the targets
    │   ├── <PDB>.pdb                    coordinate file
    │   ├── <PDB>.pdb.mca                MC-Annotate output
    │   └── <PDB>_<chain>.rmsx.in/.nch   RMSX target input files
    └── <chain>/                         one folder per scanned chain, e.g. 0/
        ├── <pdb>_<chain>.rmsx.in/.nch   RMSX inputs for this chain
        ├── sarcin-ricin_consensus.log   alignment OUTPUT (one per family)
        ├── k-turn_consensus.log
        ├── c-loop_consensus.log
        ├── e-loop_consensus.log
        └── reverse-kturn_consensus.log
```

Each `*_consensus.log` holds one alignment block per hit; RSMViewer reads every
block (across all chains of the PDB), applies the family P-value threshold, then
consolidates. To add a new PDB, drop a `rmsx_work_default/<pdb_id>/` folder in
the same structure and it becomes available to `rmv_db RNAMotifScanX`.

### Run-from-scratch mode (optional, needs binaries)

Set `data_mode` to `"run_from_scratch"` and provide working `rmsx_executable`
and `mc_annotate_executable` paths (place the compiled binaries under
`external/rmsx/bin/`). RSMViewer then runs MC-Annotate + RNAMotifScanX for the
loaded PDB. The RNAMotifScanX software and its binaries are **not** distributed
with RSMViewer.

---

## `fr3d_config.json` — FR3D

FR3D has two supported modes:

| `data_mode` | Behavior |
| --- | --- |
| `cache` | Load previously generated FR3D results for the structure from `output/fr3d_runs/`; no FR3D search or network call is made. |
| `run_from_scratch` | Run the official `fr3d-python` pipeline using the configured checkout and queries. |

The repository default is `cache`. Cache mode serves results only from
`output/fr3d_runs/`; when none exist there for the structure it reports that and
asks you to run from scratch. There is no bundled/external cache.
To run FR3D itself, change only this field to:

```json
"data_mode": "run_from_scratch"
```

### Fields

| Field | Type | Meaning / possible values |
| --- | --- | --- |
| `data_mode` | string | `"cache"` serves prior results from `output/fr3d_runs/`; `"run_from_scratch"` runs the official FR3D pipeline. |
| `fr3d_python_path` | path | Root of the fr3d-python checkout (must contain `fr3d/__init__.py` and `fr3d/search/FR3D.py`). |
| `query_path` | path | A queries directory (filtered by `query_selection`) or a single query `.json` file. Defaults to the checkout's own `fr3d/search/queries`. |
| `python_path` | path | Optional. Interpreter used to run FR3D; must have `numpy`, `scipy`, `mmcif-pdbx`. Omit to auto-detect. |
| `interactions_path` | path | Optional local interaction data directory. |
| `run_output_path` | path | Where FR3D run outputs (CSV/provenance) are written. |
| `allow_network` | bool | `false` keeps FR3D fully offline. Set `true` only to let FR3D download reference **coordinate** structures (`.cif`) that some geometric queries use as their search template. 
| `query_selection` | string | `"selected"` (default) runs only the files listed in `query_families`; `"all"` runs every query file in `query_path`. |
| `query_families` | list | Exact query names (with or without `.json`) to run when `query_selection` is `"selected"`. |
| `query_timeout_seconds` | int | Per-query timeout (≥ 10). |

### Setup and dependencies

1. Paste the official fr3d-python software into
   `external/fr3d/fr3d-python-latest/` (the path in `fr3d_python_path`). It must
   contain `fr3d/__init__.py` and `fr3d/search/FR3D.py`.
2. Run the one-shot setup:

   ```text
   rmv_setup FR3D
   ```

   This finds a suitable Python (cross-platform), installs `numpy`, `scipy`,
   and `mmcif-pdbx`, and registers FR3D. To use a specific interpreter, pass
   its path: `rmv_setup FR3D /absolute/path/to/python`. Alternatively set
   `python_path` to an interpreter that already has the dependencies.
3. Check status:

   ```text
   rmv_fr3d status
   ```

### Running FR3D

For the repository's prepared cache, use:

```text
rmv_fetch 1S72
rmv_db FR3D
```

To run the FR3D Python pipeline instead, set `"data_mode":
"run_from_scratch"` in `config/fr3d_config.json` first, then run the same
commands.

```text
rmv_fetch 1S72
rmv_db FR3D
```

RSMViewer runs the FR3D queries under `query_path` against the loaded structure
and loads the resulting motif candidates. If `rmv_db FR3D` reports that FR3D is
not ready, it prints the precise reason and tells you to run `rmv_setup FR3D`.

Note on queries: FR3D's bundled **geometric** queries define their template from
a reference PDB and therefore need `allow_network: true` to download that
coordinate file. Otherwise a query that cannot reach its reference is skipped
with a message and the remaining queries continue. This only affects reference
structures, not annotations.

---

## Quick reference

```text
rmv_db                      List sources and usage (no argument)
rmv_db RNA3DMotifAtlas      Load Atlas (API, cached)
rmv_db Rfam                 Load Rfam (API, cached)
rmv_db RNAMotifScanX        Load RMSX (preannotated or from-scratch)
rmv_db FR3D                 Load the FR3D cache, or run FR3D when data_mode is run_from_scratch
rmv_refresh                 Bypass caches and re-fetch
```
