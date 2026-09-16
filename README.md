# RSMViewer

**A PyMOL plugin for retrieving, integrating, visualizing, and comparing RNA structural motif annotations.**

![Version 3.0.0](https://img.shields.io/badge/version-2.0.0-blue)
![PyMOL 3.x+](https://img.shields.io/badge/PyMOL-2.x%2B-brightgreen)
![Platforms](https://img.shields.io/badge/platforms-macOS%20%7C%20Linux%20%7C%20Windows-lightgrey)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.6084%2Fm9.figshare.33826777-orange)](https://doi.org/10.6084/m9.figshare.33826777)

RSMViewer is a motif-centric PyMOL plugin for retrieving, consolidating,
querying, and comparing RNA structural motif (RSM) annotations. It integrates
annotations from multiple independent sources without forcing them into a single
universal classification, so agreement and disagreement between annotation
methods stay visible.

> **Who it's for:** RNA structural biologists and bioinformaticians who want to
> explore, compare, and validate RNA 3D motif annotations from different
> databases and model-based search tools directly in PyMOL — using short
> commands, without writing code.

A motif is represented internally as a **set of residues**, with no constraint
on the number of strands or chains. This lets motifs defined by sequence,
secondary structure, base-interaction patterns, or 3D geometry coexist in one
framework. Redundant residue sets are merged; distinct ones are preserved.

```text
Fetch structures  ->  Load named sources  ->  Query motifs  ->  View · Compare · Export
```

- **Version:** 2.0.0
- **Compatible with:** PyMOL 3.x and later (tested on PyMOL 3.x, macOS/Linux/Windows)
- **License:** see [LICENSE](LICENSE)

---

## Table of contents

- [Features](#features)
- [Installation in PyMOL](#installation-in-pymol)
- [Quick start](#quick-start)
- [Core concepts](#core-concepts)
- [Command reference](#command-reference)
- [Example applications](#example-applications)
- [Annotation sources](#annotation-sources)
- [FR3D setup](#fr3d-setup)
- [RNAMotifScanX setup](#rnamotifscanx-setup)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [Testing](#testing)
- [Documentation](#documentation)
- [Support](#support)
- [Citation](#citation)
- [License](#license)

---

## Features

- Load structures by PDB ID or from local `.pdb` / `.cif` files, one or many at a time.
- Four public annotation sources: `RNA3DMotifAtlas`, `Rfam`, `FR3D`, `RNAMotifScanX`.
- Atlas and Rfam retrieved from their public **APIs** with local response caching.
- FR3D executed through a user-provided checkout; RNAMotifScanX loaded from
  preannotated data or executed from scratch.
- Residue-set consolidation that preserves every source's own label and hierarchy.
- Deterministic, stable motif IDs (e.g. `1S72_00001`).
- Boolean selection across motifs, structures, and sources (`and`, `or`, `not`).
- Saved query groups, selectable PyMOL objects, medoid-based superimposition,
  and minimal coordinates-only mmCIF export.

External analysis software and large offline databases are **not** distributed
with this project; they are provided by the user under `external/`.

---

## Installation in PyMOL

1. Open **PyMOL**, then select **Plugin → Plugin Manager**.
2. Select **Install New Plugin**.
3. Under **Install from local file**, choose **Choose File** and navigate to the cloned or extracted `RSMViewer-main` directory.
4. Open the `rsmviewer` folder and select the `__init__.py` file.
5. Select **Open**, then confirm the installation when PyMOL prompts you.

If the plugin is loaded successfully, you should see the startup banner in your PyMOL console:

```text
================================================================================
RSMViewer
RNA Structural Motif Visualization and Comparative Analysis for PyMOL
Version 2.0.0 | Updated: 13 September 2026 | Compatible with PyMOL 2.x+
================================================================================
```

If the plugin does not load automatically, follow these steps:

1. Return to **Plugin → Plugin Manager → Settings**.
2. Under **Plugin override search path**, locate **RSMViewer** and move it to the top of the list.
3. Restart PyMOL if prompted. RSMViewer should then be available from the **Plugin** menu.

Verify the installation in PyMOL:

```text
rmv_help
rmv_db
```

The startup banner lists the four canonical source names.

### Terminal encoding (UTF-8)

RSMViewer prints a plain-ASCII banner and reconfigures standard output to UTF-8
on startup, so no special terminal setup is required. If your environment forces
a legacy code page (for example an Anaconda Prompt on Windows or a shell started
with `LANG=C`) and you see garbled characters in any output, enable UTF-8:

- **macOS/Linux:** start PyMOL from a UTF-8 locale, e.g.
  `LANG=en_US.UTF-8 pymol`, or run `export PYTHONUTF8=1` first.
- **Windows:** run `chcp 65001` in the prompt before launching PyMOL, or set the
  environment variable `PYTHONUTF8=1`.

---

## Quick start

Paste these seven commands into the PyMOL command line:

```text
rmv_fetch 1S72                                   # load a structure
rmv_db RNA3DMotifAtlas                            # load annotations from a source
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR
rmv_list group_SR                                 # inspect motif IDs and residues
rmv_view group_SR                                 # highlight the family
rmv_create_object group_SR                        # make selectable objects
rmv_super group_SR                                # medoid superimposition
```

**What you'll see:** `rmv_db` prints a motif-family table for 1S72; `rmv_select`
saves the Sarcin-Ricin family as `group_SR`; `rmv_list` shows each motif's stable
ID and residues; `rmv_view` colors them on the structure; and `rmv_super` reports
the medoid and per-instance RMSDs. New to the commands? Follow the
[tutorial](docs/TUTORIAL.md).

---

## Core concepts

**Structures.** `rmv_fetch` loads coordinates only; it does not retrieve
annotations. Multiple structures can be loaded in one session.

**Sources.** `rmv_db` loads one or more named sources for every fetched
structure. Source names are case-insensitive and always shown canonically. Each
source's raw records are kept separate at load time; `rmv_db` prints one compact
table per source per structure (`SELECTABLE NAME`, `ANNOTATION NAME`, `COUNT`)
and performs no overlap merging.

**Consolidation.** Residue-set consolidation is applied within the requested
scope during `rmv_select` (one family per structure) and, across families, only
when you run `rmv_combine_groups`. Two residue sets are treated as the **same**
fragment when their Jaccard index is at least `0.60`, **or** when one set is
largely nested inside the other (overlap coefficient at least `0.80`). This
containment rule lets a tight motif core from one source and an extended
annotation of the same motif from another source share a single row, while each
source keeps its own original label on that row.

**Stable motif IDs.** Each consolidated row has a deterministic ID of the form
`<PDBID>_00001`, independent of load order.

**Query groups.** `rmv_select ... as <group>` saves a snapshot of the matching
motif IDs and the original query text. Groups drive view, object creation,
superimposition, combination, coloring, and export.

---

## Command reference

| Command | Purpose |
| --- | --- |
| `rmv_fetch <ID\|path> [, <ID> ...]` | Fetch one or more structures, or load a local `.pdb`/`.cif`. |
| `rmv_db <source>[,<source>...]` | Select named sources and load annotations immediately. |
| `rmv_select <motif>, <structures>, <sources>, as <group>` | Save a Boolean query as stable motif IDs. |
| `rmv_list [ID\|group\|family]` | List all rows, one motif ID, a saved group, or every row in a family (e.g. `SARCIN-RICIN`). |
| `rmv_view <ID\|group> [, color=<name>] [, padding=<n>]` | Highlight motif residues on parent structures. |
| `rmv_hide <ID\|group\|all>` | Remove highlight (recolor the structure to neutral gray). |
| `rmv_create_object <ID\|group>` | Create selectable PyMOL objects (camera unchanged). |
| `rmv_super <ID\|group>` | Medoid-based superimposition (sequence-independent). |
| `rmv_align <ID\|group>` | Medoid-based superimposition (sequence-dependent). |
| `rmv_combine_groups <group>, <group>[, ...], as <group>` | Combine saved groups into a new group. |
| `rmv_set_color <group>, <color>` | Set a group color; preserved per source through `rmv_combine_groups`. |
| `rmv_color <motif>, <color>` | Set a motif-family color preference. |
| `rmv_bg_color <color>` | Change the background (non-motif) color. |
| `rmv_save <group\|ID\|ALL> cif` | Export minimal coordinates-only mmCIF files (prints the output directory). |
| `rmv_save current [file.png]` | Save the current PyMOL view as a high-resolution PNG (prints the path). |
| `rmv_save ALL [representation]` | Save an image of every motif (cartoon by default). |
| `rmv_colors` | List supported color names. |
| `rmv_db` (no args) | Show the four public sources and usage. |
| `rmv_source info [<N>]` | Show the active source configuration. |
| `rmv_fr3d status\|setup\|register\|run [PDB]` | Inspect / install / register / run the FR3D pipeline. |
| `rmv_rmsx status\|config\|doctor\|setup\|test\|run\|run_current` | Inspect or run the RNAMotifScanX integration. |
| `rmv_rmsx_doctor` | Diagnose the RMSX runtime and dependencies. |
| `rmv_pair <selection>` / `rmv_pair_batch <selection>` | Inspect base-pair interactions. |
| `rmv_chains` / `rmv_loaded` | Show chain diagnostics / loaded structure and source tags. |
| `rmv_refresh [PDB]` | Bypass caches and re-fetch. With no argument, refreshes every active structure; with a PDB ID, refreshes only that one. |
| `rmv_debug ON\|OFF` | Toggle verbose diagnostics (off by default). |
| `rmv_help` | Show the in-PyMOL command reference. |
| `rmv_reset` | Delete all objects and clear all caches and session state. |

The selection grammar has four comma-separated clauses:

```text
rmv_select <motif>, <structures>, <sources>, as <group>
```

- **motif** — a family name or abbreviation (`SR`, `SARCIN-RICIN`, `KT`,
  `K-TURN`, `CL`, `C-LOOP`, `EL`, `E-LOOP`, `GNRA`, …). Abbreviations,
  separators, and family-index suffixes (`sarcin-ricin-1`) all resolve to the
  same family.
- **structures** — one ID, an `and`/`or` list, or `all`.
- **sources** — a Boolean expression over source names with `not` > `and` > `or`
  precedence. A source `and` requires both sources on the **same** consolidated
  row.
- **group** — a name for the saved result (`as group_SR`).

---

## Example applications

These mirror the applications in the RSMViewer paper.

### 1 — Visualize an individual motif instance

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR
rmv_list group_SR
rmv_view 1S72_00001            # a concrete motif ID from rmv_list
```

### 2 — Visualize all instances within a family

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR
rmv_view group_SR
```

### 3 — Visualize multiple families across multiple structures

```text
rmv_fetch 1S72, 1FFK
rmv_db RNA3DMotifAtlas, Rfam
rmv_select SARCIN-RICIN, 1S72 and 1FFK, RNA3DMotifAtlas and Rfam, as group_SR
rmv_select KT, all, RNA3DMotifAtlas and Rfam, as group_KT
rmv_select CL, all, RNA3DMotifAtlas and Rfam, as group_CL
rmv_select EL, all, RNA3DMotifAtlas and Rfam, as group_EL
rmv_view group_SR, group_KT, group_CL, group_EL
```

### 4 — Superimpose a family to study structural variation

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_select SR, 1S72, RNA3DMotifAtlas, as group_SR
rmv_create_object group_SR
rmv_super group_SR
```

### 5 — Benchmark a search tool against a reference database

Uses RNA 3D Motif Atlas as ground truth for RNAMotifScanX predictions.

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas, RNAMotifScanX
rmv_select SR, 1S72, RNA3DMotifAtlas and RNAMotifScanX, as group_TP      # true positives
rmv_select SR, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP  # false positives
rmv_select SR, 1S72, RNA3DMotifAtlas and not RNAMotifScanX, as group_FN  # false negatives
```

### 6 — Inspect false positives for novel motifs

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas, RNAMotifScanX
rmv_select SR, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP
rmv_select SR, 1S72, RNA3DMotifAtlas, as group_known
rmv_set_color group_FP, red
rmv_set_color group_known, blue
rmv_combine_groups group_FP, group_known, as group_combined
rmv_create_object group_combined
rmv_super group_combined
```

> Note the `rmv_set_color` syntax: the group name follows the command with **no
> comma** after the command word — `rmv_set_color group_FP, red`.

Colors set on each source group are preserved per member when the groups are
combined: after `rmv_create_object group_combined`, the `group_FP` objects stay
red and the `group_known` objects stay blue. (An explicit color set on the
combined group itself overrides this and colors every member uniformly.)

---

## Annotation sources

| Source | Retrieval | Notes |
| --- | --- | --- |
| `RNA3DMotifAtlas` | BGSU RNA 3D Hub **API** | Online, cached locally. |
| `Rfam` | Rfam **API** | Online, cached locally. |
| `FR3D` | Local cache or user-provided checkout | Loads `external/fr3d/fr3d_cache/` by default, or runs official fr3d-python in `run_from_scratch` mode. |
| `RNAMotifScanX` | Preannotated data or from-scratch run | See setup below. |

Atlas and Rfam require network access on first retrieval; responses are cached
so repeated queries are offline. Use `rmv_refresh` to bypass the cache.

---

## FR3D setup

FR3D has two modes, selected by `data_mode` in `config/fr3d_config.json`:

- **`cache` (default):** `rmv_db FR3D` loads previously generated results — it
  first looks for a prior run under `output/fr3d_runs/<PDB>/`, and if none is
  found there it falls back to the bundled cache at `external/fr3d/fr3d_cache/`.
  It never runs FR3D or makes a network call.
- **`run_from_scratch`:** `rmv_db FR3D` always executes the official BGSU
  fr3d-python pipeline on the loaded structure and loads the freshly generated
  results from a new timestamped run directory under `output/fr3d_runs/`.

`rmv_reset` clears the FR3D run cache at `output/fr3d_runs/`.

To run the pipeline from scratch:

1. Paste the official fr3d-python checkout into `external/fr3d/fr3d-python-latest/`
   (download from <https://github.com/BGSU-RNA/fr3d-python/tree/latest>). It must
   contain `fr3d/__init__.py` and `fr3d/search/FR3D.py`.
2. Set `"data_mode": "run_from_scratch"` in `config/fr3d_config.json`.
3. Install FR3D's dependencies and register it in one step — this also prints a
   readiness report:

   ```text
   rmv_setup FR3D
   ```

4. Run it on a loaded structure:

   ```text
   rmv_fetch 1S72
   rmv_db FR3D
   ```

> **Detailed FR3D instructions** — both modes, query selection, the full config
> reference, first-run cache behavior, and troubleshooting — are in
> [external/FR3D_SETUP.md](external/FR3D_SETUP.md).

---

## RNAMotifScanX setup

RNAMotifScanX (RMSX) is controlled by `data_mode` in `config/rmsx_config.json`.
The default is `preannotated`, which loads precomputed results and needs no
binaries.

### Preannotated mode (default)

1. Download the preannotated bundle from Figshare:
   **[https://doi.org/10.6084/m9.figshare.33826795](https://doi.org/10.6084/m9.figshare.33826795)**
2. Extract `rmsx_preannotated_input_output.tar.gz` inside
   `external/rmsx_preannotated/` so the extracted folder sits there:

   ```text
   external/rmsx_preannotated/rmsx_work_default/
   ```

   RSMViewer can also read the `.tar.gz` archive directly as a fallback; the
   extracted folder is faster.
3. Keep `"data_mode": "preannotated"` and run:

   ```text
   rmv_fetch 1S72
   rmv_db RNA3DMotifAtlas, RNAMotifScanX
   ```

### Run-from-scratch

To generate results with the real scanner instead of reading precomputed logs,
set `"data_mode": "scan_prepared"` (or `"run_from_scratch"`) in
`config/rmsx_config.json` and provide the binaries, then run `rmv_db RNAMotifScanX`.

> **Detailed RMSX instructions** are in [external/RMSX_FROM_SCRATCH.md](external/RMSX_FROM_SCRATCH.md).


---

## Configuration

Both pipelines are configured by small JSON files under `config/`. See
[config/README.md](config/README.md) for every field, allowed value, and how
paths are resolved.

- P-value cutoffs live only under `pvalue_thresholds` in `config/rmsx_config.json`.
  Command-line P-value overrides are intentionally not supported.
- The residue-overlap thresholds are code-defined, not config options. To change
  them, edit `DEFAULT_JACCARD_THRESHOLD` (0.60) and `DEFAULT_CONTAINMENT_THRESHOLD`
  (0.80) in `rsmviewer/database/consolidated_table.py`.

---

## Troubleshooting

**Multiple PDB IDs not loading.** Use comma-separated IDs:
`rmv_fetch 1S72, 1FFK, 3CC2`.

**Atlas/Rfam returns no annotations.** Check network access, then
`rmv_refresh` and re-run `rmv_db`.

**A source `and` query is empty.** The two sources may annotate the motif at
residue sets that neither overlap by ≥ 0.60 Jaccard nor by ≥ 0.80 containment,
so they sit on separate rows. Inspect with `rmv_list`.

**FR3D is unavailable.** Confirm the checkout path in `config/fr3d_config.json`,
run `rmv_setup FR3D`, and check `rmv_fr3d status`. Geometric queries need
`allow_network: true`.

**RMSX returns no motifs.** A result can legitimately contain zero accepted
motifs when every reported P-value exceeds the configured threshold. Check
`data_mode`, the preannotated archive/layout, and `pvalue_thresholds`.

---

## Testing

Run from the project root:

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q rsmviewer tests
RSMVIEWER_ROOT="$PWD" pymol -cq tests/pymol_applications_e2e.py
RSMVIEWER_ROOT="$PWD" pymol -cq tests/pymol_smoke.py
```

`tests/pymol_applications_e2e.py` exercises all six example applications end to
end and prints a PASS/FAIL report.

---

## Documentation

Pick the document that matches what you need:

| I want to… | Read |
| --- | --- |
| Learn RSMViewer step by step | [docs/TUTORIAL.md](docs/TUTORIAL.md) |
| Look up exact command syntax | [Command reference](#command-reference) · in-PyMOL `rmv_help` |
| Understand the query grammar and consolidation | [docs/SUPPLEMENTARY_COLLECTION_AND_TABLE.md](docs/SUPPLEMENTARY_COLLECTION_AND_TABLE.md), [docs/SUPPLEMENT.md](docs/SUPPLEMENT.md) |
| Configure FR3D / RNAMotifScanX | [config/README.md](config/README.md) |
| Understand the design and internals | [DEVELOPED.md](DEVELOPED.md) |
| Read the paper and its supplement | [docs/MANUSCRIPT.md](docs/MANUSCRIPT.md), [docs/SUPPLEMENT.md](docs/SUPPLEMENT.md) |

## Support

- Run `rmv_help` in PyMOL for the full command reference.
- Enable `rmv_debug ON` to print diagnostic messages before reporting an issue.
- When reporting a problem, include your PyMOL version, operating system, the
  exact commands you ran, and the console output.

## Citation

If you use RSMViewer in your research, please cite it. Machine-readable citation
metadata is in [CITATION.cff](CITATION.cff).

> Rahad, R. H., Pranjal, S., Khan, N. S., Zhang, S., & Zhong, C. RSMViewer: A
> PyMOL plugin for RNA structural motif visualization. *Bioinformatics*.
> DOI: [10.6084/m9.figshare.33826777](https://doi.org/10.6084/m9.figshare.33826777)

## License

RSMViewer is released under the MIT License — see [LICENSE](LICENSE). External
analysis software and datasets under `external/` are provided by the user and
retain their own licenses.
