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
folder**. For example, `../external/rmsx_preannotated/rmsx_work_default` points to
`<project>/external/rmsx_preannotated/rmsx_work_default`. You may also use absolute paths.

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

RMSX has two modes, selected by `data_mode`:

- **`preannotated`** (default) — reads precomputed consensus logs. To provide the
  most up-to-date RNAMotifScanX annotations, RSMViewer collects them live from
  our server, downloading each requested PDB's results the first time; no binary
  is executed.
- **`run_from_scratch`** — runs RNAMotifScanX on the PDB. The first time you use
  it, all required files and third-party software are downloaded automatically
  from `rmsx_bundle_url` into `external/rmsx/`.

### Fields

| Field | Type | Meaning / possible values |
| --- | --- | --- |
| `data_mode` | string | `"preannotated"` (default) loads precomputed annotations; `"run_from_scratch"` runs RNAMotifScanX. |
| `scan_runtime` | string | How `run_from_scratch` runs the scanner: `"auto"` (default: native, then WSL2 on Windows, then Docker) or force `"native"`, `"wsl"`, `"docker"`. Inspect with `rmv_rmsx_doctor`. |
| `rmsx_bundle_url` | URL | Where the RNAMotifScanX files are downloaded from on first `run_from_scratch` use. Default: `https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz`. |
| `rmsx_executable` | path | Optional `scan` program of your own. Leave empty to use the one RSMViewer prepares. |
| `scan_chains` | list | Chains to run in `run_from_scratch`, e.g. `["9"]`. `[]` (default) = every chain. |
| `query_motifs_dir` | path | Optional directory of your own `<family>_consensus.struct` query models. Empty = the standard set. |
| `num_threads` | int | Threads passed to `scan` (default 4). |
| `scan_timeout_seconds` | int | Time limit per (chain, family) scan (default 21600 = 6 hours). |
| `preannotated_base_url` | URL | Server folder holding one `<pdb_lowercase>.tar.gz` per PDB. `preannotated` mode downloads `<base>/<pdb>.tar.gz` when the PDB is not already in `pdb_prebuild_dir`. Default: `https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default`. |
| `pdb_prebuild_dir` | path | Local per-PDB data folder (`external/rmsx_preannotated/rmsx_work_default`); downloaded results are extracted here. |
| `output_dir` | path | Where per-family result logs are written/read (`../output/rmsx_results`). |
| `motif_families` | list | Families to load: `k-turn`, `c-loop`, `sarcin-ricin`, `reverse-kturn`, `e-loop`. |
| `pvalue_thresholds` | object | Per-family P-value cutoff; a hit is kept when its P-value is `<=` the cutoff. Lower = stricter (fewer, higher-confidence hits). Omit a family to keep its paper default (0.05 if the family is unknown). Names are matched regardless of spelling (`KINK-TURN`/`K-TURN`, `REVERSE-KINK-TURN`/`REVERSE-K-TURN`, ...). Applied to preannotated and from-scratch results, and to single- and multi-source `rmv_db`; the file is re-read on every `rmv_db`. |

Note: the residue-set redundancy threshold (Jaccard 0.60) is **not** a config
value; it is defined in code (`rsmviewer/database/consolidated_table.py`,
`DEFAULT_JACCARD_THRESHOLD`).

### Preannotated mode (default, no binaries needed)

Run:

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas,RNAMotifScanX
rmv_select SR, 1S72, RNA3DMotifAtlas and RNAMotifScanX, as group_TP
```

For the requested PDB, RSMViewer uses the first of these that has data:

1. `pdb_prebuild_dir/<pdb>/` (already downloaded or placed by you);
2. a download of `<preannotated_base_url>/<pdb>.tar.gz` (lowercase PDB ID),
   extracted into `pdb_prebuild_dir` (needs an internet connection; later loads
   are local). It is attempted only for a 4-character PDB ID with no local
   results, and accepted only if the archive extracts safely and holds at least
   one `*_consensus.log`.

It then copies the matching family logs into `output_dir` and loads them. The
Jaccard/containment consolidation aligns them with the other sources. If the
server has no results for a PDB (the dataset may not cover every structure yet),
RSMViewer reports that; use `run_from_scratch` to scan it yourself.

### Preannotated data layout

Inside each downloaded archive (and in `pdb_prebuild_dir`), data is keyed by PDB id and chain:

```text
rmsx_work_default/
└── <pdb_id_lowercase>/                 e.g. 1s72/
    └── <chain>/                         one folder per chain, e.g. 0/
        ├── sarcin-ricin_consensus.log   alignment output (one per family)
        ├── k-turn_consensus.log
        ├── c-loop_consensus.log
        ├── e-loop_consensus.log
        └── reverse-kturn_consensus.log
```

Each `*_consensus.log` holds one alignment block per hit; RSMViewer reads every
block (across all chains of the PDB), applies the family P-value threshold, then
consolidates. To add a new PDB, drop a `rmsx_work_default/<pdb_id>/` folder in
the same structure and it becomes available to `rmv_db RNAMotifScanX`.

### Run-from-scratch mode (optional)

Set `data_mode` to `"run_from_scratch"`, then:

```text
rmv_fetch 1S72
rmv_db RNAMotifScanX
```

The first time you run it, RSMViewer downloads all required files and third-party
software by itself from `rmsx_bundle_url` into `external/rmsx/`, prepares them for
your computer, and runs. See [`external/rmsx_setup.md`](../external/rmsx_setup.md)
for the requirements on macOS, Windows and Linux and for troubleshooting.

---

## `fr3d_config.json` — FR3D

FR3D has two supported modes:

| `data_mode` | Behavior |
| --- | --- |
| `cache` | Load previously generated FR3D results for the structure from `output/fr3d_runs/`; no FR3D search or network call is made. |
| `run_from_scratch` | Run the official `fr3d-python` pipeline using the configured checkout and queries. |

The repository default is `run_from_scratch`, which needs your FR3D checkout
(see the setup below). Cache mode serves results only from `output/fr3d_runs/`;
when none exist there for the structure it reports that and asks you to run from
scratch. There is no bundled/external cache. To use cache mode, set:

```json
"data_mode": "cache"
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

To load FR3D results from an earlier run of the structure (cache mode), use:

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
