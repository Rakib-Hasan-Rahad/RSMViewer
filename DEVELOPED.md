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
`current_source_names`, selects adapters, and loads every fetched structure. In
combined mode, `_load_combined_motifs` fetches each source, runs within-source
dedup, populates the `ConsolidatedAnnotationTable`, and merges with
`ResidueMerger`. When RNAMotifScanX (source 7) is among the sources,
`_ensure_rmsx_preannotated` ingests its preannotated results first.

### `rmv_select`

`rsmviewer/database/query_parser.py` parses four clauses into a
`QueryExpression` (motif, structures, sources predicate, group, text).
`gui.select_annotation_query` matches rows by the source predicate and
`labels_match_motif`, then saves the motif IDs under the group name.

### `rmv_list` / `rmv_view` / `rmv_hide`

Read `ConsolidatedAnnotationTable` rows. `rmv_list` resolves its argument as a
saved group, a stable motif ID, or a motif-family name (e.g. `SARCIN-RICIN`);
for a family it lists every row any source labels as that family. `rmv_view`
highlights parent-structure residues (optionally with `color=` and `padding=`)
without copying objects; `rmv_hide` recolors to neutral gray.

### `rmv_create_object` / `rmv_super` / `rmv_combine_groups`

`create_annotation_objects` builds `motif_<id>` objects and colors each one by
its group. `rmv_combine_groups` unions the motif IDs of query groups into a new group
and records each member's origin group in `member_colors`, so per-source colors
set with `rmv_set_color` survive the combine: `_group_member_color_key` colors
each object by its origin group unless an explicit color is set on the combined
group itself. `rmv_super` copies objects temporarily, computes pairwise RMSD,
chooses the minimum-average-RMSD medoid, and transforms the selected objects.

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

Runner: `rsmviewer/tools/rmsx_runner.py`. `config/rmsx_config.json` controls
`data_mode`, executable paths, the preannotated archive/directory, motif
families, output directory, and `pvalue_thresholds`.

- **Preannotated mode:** `copy_preannotated_results` extracts the matching
  `*_consensus.log` outputs for the PDB. It concatenates **all chains** of a
  family into one result file (so no chain overwrites another).
- **From-scratch mode (`run_pipeline`):** prefers prebuilt `.rmsx.in` inputs
  from `pdb_prebuild_dir` then `pdb_prebuild_archive`
  (`_copy_prebuilt_targets_from_directory` / `_extract_prebuilt_targets_from_archive`);
  only when none are found does it run MC-Annotate to generate inputs. It then
  runs the RMSX `scan` step.

The converter accepts both tabular RMSX rows and alignment-report logs
(`Aligning`, `Alignment score`, `P-value` blocks), applying per-family P-value
thresholds.

---

## Caching

- `rsmviewer/database/cache_manager.py` caches provider API responses outside the
  install directory (structure ID, adapter, retrieval time, expiry, provenance).
- `rsmviewer/database/motif_hierarchy_cache.py` is a text-keyed SQLite display
  cache using canonical `source_key` values. It auto-migrates a legacy
  integer-keyed schema, and exposes `close()` / `close_hierarchy_cache()` so
  `rmv_reset` can drop and reopen a fresh connection.
- `rsmviewer/tools/rmsx_runner.py` maintains a per-PDB **preannotated extraction
  cache** at `output/rmsx_results/.preannotated_cache/<pdb_id>/`. Enumerating
  members of the large gzip archive requires decompressing the whole stream, so
  each PDB's small `*_consensus.log` files are extracted once and reused on
  later loads and across PyMOL sessions. The cache is stamped with the source
  archive/directory identity (path, mtime, size); a changed source invalidates
  it automatically.
- A snapshot of the SQLite hierarchy cache and the preannotated extraction cache
  is shipped in the repository so a fresh clone works immediately; both are
  regenerated on demand.

`rmv_reset` clears **all** caches and session state in one call: it deletes every
PyMOL object, resets session variables, clears the SQLite hierarchy cache (data
and file, including `-wal`/`-shm`), the on-disk API response cache, each
provider's in-process memory cache, the preannotated RMSX extraction cache, the
motif loader, and custom color assignments.

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
    └── rmsx_runner.py            RMSX preannotated + from-scratch pipeline
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
