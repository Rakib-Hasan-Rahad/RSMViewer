# RSMViewer Developer Guide

This document describes the current implementation for maintainers. End-user
installation and command usage live in [README.md](README.md) and
[docs/TUTORIAL.md](docs/TUTORIAL.md). Configuration fields are in
[config/README.md](config/README.md). The detailed collection and table
semantics are in [docs/SUPPLEMENTARY_COLLECTION_AND_TABLE.md](docs/SUPPLEMENTARY_COLLECTION_AND_TABLE.md).

- [Architecture](#architecture)
- [Runtime entry point](#runtime-entry-point)
- [Public vs internal identity](#public-vs-internal-identity)
- [Core data model](#core-data-model)
- [Residue-set consolidation](#residue-set-consolidation)
- [Motif name normalization](#motif-name-normalization)
- [Command flow](#command-flow)
- [External pipelines](#external-pipelines)
- [Caching](#caching)
- [Module map](#module-map)
- [Testing](#testing)
- [Conventions](#conventions)

---

## Architecture

```text
PyMOL command
  -> rsmviewer/gui.py                     command handlers + session state
       -> source registry / provider selection
       -> per-structure ConsolidatedAnnotationTable
       -> PyMOL view / object / super operations
       -> export and alignment services

Provider data
  -> MotifInstance / ResidueSpec
  -> ResidueMerger (residue-set redundancy filtering)
  -> ConsolidatedAnnotationTable (stable IDs + source hierarchy)
  -> rmv_select query-group snapshot
```

---

## Runtime entry point

`rsmviewer/plugin.py` is loaded by PyMOL. `__init_plugin__` (called with `None`
in headless tests):

1. Configures UTF-8-safe console output.
2. Prints the startup banner (version 2.0.0, the four source names).
3. Initializes the database registry.
4. Registers all commands via `initialize_gui()`.

The GUI singleton is created by `get_gui()` and exposed as the module-level
`rsmviewer.gui.gui`.

---

## Public vs internal identity

Public source identity is deliberately separate from adapter identity:

| Public source | Internal adapter |
| --- | --- |
| `RNA3DMotifAtlas` | BGSU RNA 3D Hub API adapter |
| `Rfam` | Rfam API adapter |
| `FR3D` | external FR3D adapter/runner |
| `RNAMotifScanX` | external RMSX adapter/runner |

Internal adapters are routed through `SourceRegistry` and `SOURCE_ID_MAP` (keys
`3` Atlas, `4` Rfam, `5` FR3D, `7` RNAMotifScanX). No public command exposes
adapter names, numeric IDs, or numeric source suffixes.

---

## Core data model

### `ResidueSpec`

`rsmviewer/database/base_provider.py`. Stores chain, residue number, nucleotide,
insertion code, and model.

### `MotifInstance`

A provider result: instance ID, source motif label, structure ID, residue list,
annotation text, and source metadata.

### `ConsolidatedAnnotationTable`

`rsmviewer/database/consolidated_table.py`. One table per structure. Each row
(`AnnotationRow`) holds:

- deterministic `motif_id` (`<PDBID>_00001`)
- `structure_id`
- normalized residue set (chain, number, insertion code, model)
- per-source labels (`source_annotations`)
- per-source hierarchy (`source_hierarchy`)
- retrieval provenance

Rows are materialized deterministically: `_reindex_structure` sorts by
`(structure_id, residue_set)` and assigns sequential five-digit IDs, so IDs do
not depend on provider load order or cache sequence.

---

## Residue-set consolidation

Two residue sets are merged into one row when **either** condition holds:

- Jaccard similarity ≥ `DEFAULT_JACCARD_THRESHOLD` (0.60), or
- overlap coefficient `|A∩B| / min(|A|,|B|)` ≥ `DEFAULT_CONTAINMENT_THRESHOLD`
  (0.80).

The containment rule captures nested annotations — a tight motif core from one
source inside an extended region from another — which Jaccard alone would score
too low. Both constants are defined in `consolidated_table.py`;
`_find_matching_row` prefers the strongest match and breaks ties by `motif_id`.

`rsmviewer/database/residue_merger.py` (`ResidueMerger`) performs the
name-agnostic, residue-set-based merge used by the combined-source loader:
overlap is decided only by residues/chains (exact match, subset/superset, or
Jaccard ≥ threshold), never by motif label. Each source's own label is preserved
in the hierarchy cache. Sources are processed right-to-left so the leftmost
(highest-precedence) source wins ties.

> Terminology: this project uses **residue-based merging**. The former
> "cascade merge" name and module no longer exist.

---

## Motif name normalization

`rsmviewer/database/motif_aliases.py` is the single source of truth for motif
name handling:

- `canonical_motif(value)` — resolves abbreviations (`SR`, `KT`, `CL`, `EL`),
  separators, and family-index suffixes (`sarcin-ricin-1`) to a canonical family
  name; unknown families fall back to a normalized form.
- `labels_match_motif(query, labels)` — exact canonical match, with a guarded
  substring fallback for compound/novel labels. Distinct families such as
  `K-TURN` and `REVERSE-K-TURN` never cross-match.
- `is_known_family(value)` — whether a value resolves to a known family.

`query_parser.canonical_motif_name` and `gui.select_annotation_query` both
delegate here, so adding an alias in one place updates the whole system.

---

## Command flow

### `rmv_fetch`

`fetch_raw_pdb` accepts one PDB ID, arbitrary comma-separated IDs, or local
`.pdb`/`.cif` paths. Each structure is stored in `gui.loaded_structures` with one
table in `gui.annotation_tables`.

### `rmv_db`

`select_database` parses canonical names via `SourceRegistry`, records
`current_source_names`, selects adapters, and loads every fetched structure. Each
source's raw records are stored separately in the structure's
`ConsolidatedAnnotationTable` (`merge_enabled=False`): `rmv_db` performs no
containment/Jaccard merging, so overlapping and nested annotations are preserved
for later `rmv_select` / `rmv_combine_groups`. The load prints one compact table
per source per structure (`SELECTABLE NAME`, `ANNOTATION NAME`, `COUNT`) via
`_family_source_breakdown`. When RNAMotifScanX (source 7) is among the sources,
`_ensure_rmsx_results` first makes its results available per `data_mode`
(preannotated: local folder, else a per-PDB download from the results server;
run_from_scratch: a real scan; see below).

### `rmv_select`

`rsmviewer/database/query_parser.py` parses four clauses into a
`QueryExpression` (motif, structures, sources predicate, group, text).
`gui.select_annotation_query` gathers the raw rows any source labels as the
family, consolidates them per structure with `ResidueMerger` (containment +
Jaccard), evaluates the Boolean source predicate against each merged row's
contributing sources, and saves the resulting motif IDs under the group name.
Each merged row keeps every contributing database's original label.

### `rmv_list` / `rmv_view` / `rmv_hide`

Read `ConsolidatedAnnotationTable` rows. `rmv_list` resolves its argument as a
saved group, a stable motif ID, or a motif-family name (e.g. `SARCIN-RICIN`);
for a family it lists every row any source labels as that family. `rmv_view`
sets the whole target structure(s) to gray80 (all atoms, including protein and
ligands) and then highlights the parent-structure residues (optionally with
`color=` and `padding=`) without copying objects; for several targets it grays
the involved structures once so each highlight is preserved. `rmv_hide` recolors
the structure back to neutral gray.

### `rmv_create_object` / `rmv_super` / `rmv_combine_groups`

`create_annotation_objects` builds `motif_<id>` objects and colors each one by
its group. `rmv_combine_groups` unions the motif IDs of query groups into a new
group and preserves each database's original source labels as separate columns
(never the input group names); overlapping merged rows keep every contributing
database's label, including differing family assignments. It records each
member's origin group in `member_colors`, so per-source colors set with
`rmv_set_color` survive the combine: `_group_member_color_key` colors each object
by its origin group unless an explicit color is set on the combined group itself.
`rmv_super` copies objects temporarily, computes pairwise RMSD, chooses the
minimum-average-RMSD medoid, and transforms the selected objects. Every
superimposed instance inherits one color (the group's session color if set,
otherwise the motif-family color); the view auto-orients onto the superimposed
section and only the superimposed objects remain enabled.

### `rmv_set_color`

`rmv_set_color <group>, <color>` (no comma after the command) stores a color
preference and recolors any existing member objects; it also feeds
`rmv_create_object`. It is a thin alias over `rmv_set color, ...`.

---

## External pipelines

### FR3D

Runner: `rsmviewer/tools/fr3d_search_runner.py`. It runs the user-provided
official fr3d-python checkout **without modifying it**, using an in-memory
autofix for compile-blocking empty blocks and a `cif_local` target so FR3D
annotates the exact loaded structure. Configuration: `config/fr3d_config.json`.

`gui.register_fr3d_source` validates the checkout and an interpreter that can
import `numpy`, `scipy`, `mmcif-pdbx`, and `fr3d`. `run_fr3d_search` discovers
queries (`query_path`), stages the target CIF, runs each query, and ingests the
resulting CSVs. `FR3DConverter` converts CSV/TXT to `MotifInstanceSimple`.
Bundled geometric queries reference an external PDB template and require
`allow_network: true`.

### RNAMotifScanX

Modules: `rsmviewer/tools/rmsx_runner.py` (preannotated data and running scans) and `rsmviewer/tools/rmsx_runtime.py` (downloading the runtime bundle, finding,
building and invoking the `scan` executable per platform, plus setup and
diagnostics).
`config/rmsx_config.json` controls `data_mode`, `scan_runtime`, paths, motif
families, output directory, and `pvalue_thresholds`. Commands: `rmv_db
RNAMotifScanX` (loads results per `data_mode`), `rmv_setup RNAMotifScanX`
(`MotifVisualizerGUI.setup_rmsx_runtime`), and `rmv_rmsx_doctor`
(`MotifVisualizerGUI.rmsx_doctor`). Both single- and multi-source `rmv_db` go
through `MotifVisualizerGUI._ensure_rmsx_results`, which re-reads the config.

- **Preannotated mode:** `MotifVisualizerGUI._copy_preannotated_rmsx` resolves a
  PDB's results in this order and stops at the first source with data:
  1. the local `pdb_prebuild_dir` (`external/rmsx_preannotated/rmsx_work_default/<pdb>/`);
  2. `download_preannotated_pdb`, which fetches
     `<preannotated_base_url>/<pdb_lowercase>.tar.gz` (default
     `https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default`)
     and extracts it into `pdb_prebuild_dir`.

  The preannotated data is collected live from the project's server so that the
  most up-to-date RNAMotifScanX annotations are used. The download is staged in
  a temporary directory and moved into place only when complete; unsafe archive
  paths, links, and archives with no `*_consensus.log` are rejected; an existing
  `<pdb>/` folder is merged into, never replaced. TLS is verified (system store,
  then `certifi`); only if both fail is an unverified connection used, with a
  warning. An existing local folder is
  never re-downloaded, so to refresh a PDB delete `<pdb>/` under
  `pdb_prebuild_dir`. A PDB the server does not have (HTTP 404) produces an
  explicit message rather than a silent empty result.
  `copy_preannotated_results` then reads the matching `*_consensus.log`
  files straight from that folder (no extraction cache) and concatenates **all chains** of a family into one result file (so
  no chain overwrites another).
- **From-scratch mode:** `_run_rmsx_from_scratch` calls
  `rmsx_runner.run_scan_prepared`, synchronously, and never falls back to preannotated data. Output goes to
  `output/rmsx_results/run_from_scratch/<PDB>_<stamp>/`; a same-session repeat
  reuses it (`_rmsx_scan_runs`, cleared by `rmv_refresh` and `rmv_reset session`).
  Query models are searched in `Queries/reduced` before `Queries` because the
  reduced set reproduces the published results.
- **Runtime bundle (first use):** the scanner source, scoring matrices, query
  models, Linux binary and third-party tools are not in the repository
  (`external/rmsx/` is git-ignored except `README.md` and `bin/.gitkeep`).
  `rmsx_runtime.ensure_runtime_bundle` downloads `rmsx_bundle_url` (default
  `https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz`) when
  `bundle_ready` is false, verifies it against `<url>.sha256`, unpacks it in a
  temporary folder (unsafe paths rejected, symlinks/devices skipped), and moves
  files into `external/rmsx/` without overwriting existing ones.
  On Apple Silicon the archive includes a prebuilt `bin/macos-arm64/scan`;
  `_try_native` probes it and uses it if it starts (it needs Homebrew's arm64
  Boost), otherwise `setup` builds a native scanner.
  `run_scan_prepared` calls it, then `resolve_runtime`, and if no scanner runs yet
  calls `setup` automatically, so the first `rmv_db RNAMotifScanX` prepares
  everything. `rmv_setup RNAMotifScanX` does the same ahead of time.
- **Scanner runtime** (`rmsx_runtime.resolve_runtime`, `scan_runtime` config):
  *native* (a `scan` that runs on this OS/CPU: the downloaded ELF on Linux x86-64,
  or one built by `build_native_scan` into `bin/<platform>/`), then *WSL2*
  (Windows), then *Docker* (`ubuntu:22.04`, `linux/amd64`). `build_command`
  builds the exact command for each runtime and translates paths (`to_wsl_path`;
  Docker bind mounts, one per directory). `build_native_scan` compiles the 8
  `.cc` files with the system C++ compiler and Boost (Homebrew is installed on
  macOS when needed), links only the Boost libraries that exist (`boost_system`
  is header-only in current Boost), and retries the link against older macOS
  SDKs because a Command Line Tools SDK newer than the OS can have unreadable
  `.tbd` stubs. `setup` tries native, then a native build, then WSL2, then
  Docker (starting Colima if it is installed but stopped), verifying each by
  starting the scanner. The Windows routes are implemented but unverified on a
  real Windows machine.
- **Scanner crashes:** `scan` can segfault part-way through a large structure
  (also seen on the cluster that produced the preannotated logs, which end with
  `# scan failed rc=-11`). `_run_one_scan` keeps the complete alignments printed
  before the crash (`_complete_alignments`), marks the log the same way, and the
  run is reported in `partial_runs`; a crash with no complete alignment fails.
- **P-values** are random estimates (`scan` seeds a simulation with the clock),
  so borderline hits vary between runs; alignments and scores are deterministic.

The converter accepts both tabular RMSX rows and alignment-report logs
(`Aligning`, `Alignment score`, `P-value` blocks), applying per-family P-value
thresholds (a record is kept when `p_value <= threshold`).

P-value thresholds come from `pvalue_thresholds` in `config/rmsx_config.json`,
falling back to the paper defaults in `RMSX_PVALUE_THRESHOLDS`
(`user_annotations/converters.py`) and then to 0.05 for an unknown family.
Family names are matched by `rmsx_pvalue_threshold` regardless of spelling
(`KINK-TURN`/`K-TURN`, `REVERSE-KINK-TURN`/`REVERSE-K-TURN`,
`reverse-kturn_consensus`, ...). The config is re-read on every `rmv_db`, single-
or multi-source, and replaces any earlier values, so an edit takes effect on the
next `rmv_db` and a family deleted from the file returns to its paper default.

---

## Caching

- `rsmviewer/database/cache_manager.py` caches provider API responses outside the
  install directory (structure ID, adapter, retrieval time, expiry, provenance).
- `rsmviewer/database/motif_hierarchy_cache.py` is a text-keyed SQLite display
  cache using canonical `source_key` values. It auto-migrates a legacy
  integer-keyed schema, and exposes `close()` / `close_hierarchy_cache()` so
  `rmv_reset` can drop and reopen a fresh connection.
- RMSX preannotated data is **not** a cache: it is downloaded per PDB into
  `external/rmsx_preannotated/rmsx_work_default/<pdb>/` and read from there, so
  `rmv_reset cache` leaves it in place. (It only removes the
  `output/rmsx_results/.preannotated_cache/` folder that older versions created.)
- A snapshot of the SQLite hierarchy cache is shipped in the repository so a
  fresh clone works immediately; it is regenerated on demand.

`rmv_reset` requires an explicit subcommand to actually reset anything. With no
argument it only prints details about the two subcommands below and performs
no reset:

- `rmv_reset cache` — clears the caches: the SQLite hierarchy cache (data and
  file, including `-wal`/`-shm`), the on-disk API response cache, and each
  provider's in-process memory cache (RMSX results downloaded into
  `rmsx_work_default/` are data, not cache, and are kept). Loaded objects, query groups, and other session state are left
  untouched.
- `rmv_reset session` — deletes all PyMOL objects and resets session state
  (loaded structures, query groups, source selections, motif loader, custom
  colors) to defaults. Caches on disk are left untouched.

---

## Module map

```text
rsmviewer/
├── plugin.py                     PyMOL entry point + banner
├── gui.py                        command handlers, session state, pipelines
├── alignment.py                  rmv_super / rmv_align + motif alias helpers
├── colors.py                     color assignment and PyMOL coloring
├── structure_exporter.py         minimal mmCIF export
├── database/
│   ├── base_provider.py          MotifInstance / ResidueSpec
│   ├── consolidated_table.py     ConsolidatedAnnotationTable + thresholds
│   ├── residue_merger.py         ResidueMerger (residue-based merge)
│   ├── motif_aliases.py          canonical motif-name normalization
│   ├── query_parser.py           rmv_select grammar
│   ├── motif_hierarchy_cache.py  text-keyed display cache
│   ├── cache_manager.py          API response cache
│   ├── bgsu_api_provider.py      RNA3DMotifAtlas (BGSU API)
│   ├── rfam_api_provider.py      Rfam (API)
│   └── user_annotations/         FR3D / RMSX converters + provider
└── tools/
    ├── fr3d_search_runner.py     official FR3D runner (unmodified checkout)
    ├── rmsx_runner.py            RMSX preannotated data, scan execution
    └── rmsx_runtime.py           RMSX scanner runtime: find/build/verify (native, WSL2, Docker), setup, doctor
```

---

## Testing

From the project root:

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q rsmviewer tests
RSMVIEWER_ROOT="$PWD" pymol -cq tests/pymol_applications_e2e.py
RSMVIEWER_ROOT="$PWD" pymol -cq tests/pymol_smoke.py
```

Key tests:

- `tests/pymol_applications_e2e.py` — all six applications, PASS/FAIL report.
- `tests/test_consolidated_table.py` — Jaccard + containment merge, stable IDs.
- `tests/test_motif_aliases.py` — abbreviation/variant resolution, K-TURN vs
  REVERSE-K-TURN separation.
- `tests/test_query_parser.py` — selection grammar and precedence.
- `tests/pymol_fetch_multi_smoke.py` — comma-separated multi-fetch.

---

## Conventions

- Public commands stay limited to the registered named-source, motif-ID, and
  group interface. Do not expose adapter names or numeric IDs.
- Residue-overlap thresholds are code-defined (`consolidated_table.py`), not
  config options.
- P-value cutoffs live only in `config/rmsx_config.json` under
  `pvalue_thresholds`.
- User-facing text refers to Atlas/Rfam retrieval as **API** access.
- The merge concept is **residue-based merging**; the "cascade" name is retired.
