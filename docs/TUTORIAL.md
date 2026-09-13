# RSMViewer Tutorial

A guided walkthrough of the RSMViewer workflow, from installation to the six
example applications and external-pipeline setup. Commands are entered in the
PyMOL command line.

For the detailed collection and consolidated-table rules, see
[SUPPLEMENTARY_COLLECTION_AND_TABLE.md](SUPPLEMENTARY_COLLECTION_AND_TABLE.md).

- [1. Install](#1-install)
- [2. Load structures](#2-load-structures)
- [3. Load annotation sources](#3-load-annotation-sources)
- [4. Query motifs](#4-query-motifs)
- [5. List and view](#5-list-and-view)
- [6. Create objects and superimpose](#6-create-objects-and-superimpose)
- [7. Combine and color groups](#7-combine-and-color-groups)
- [8. Export motifs](#8-export-motifs)
- [9. The six example applications](#9-the-six-example-applications)
- [10. FR3D pipeline](#10-fr3d-pipeline)
- [11. RNAMotifScanX pipeline](#11-rnamotifscanx-pipeline)
- [12. Diagnostics and validation](#12-diagnostics-and-validation)
- [13. Common problems](#13-common-problems)

---

## 1. Install

1. Open PyMOL.
2. Open **Plugin → Plugin Manager → Settings → Add new directory**.
3. Select the directory containing this project and restart PyMOL.
4. Confirm the plugin loaded:

   ```text
   rmv_help
   rmv_db
   ```

The four public sources are `RNA3DMotifAtlas`, `Rfam`, `FR3D`, and
`RNAMotifScanX`. Atlas and Rfam use online **API** retrieval with local caching.
FR3D and RNAMotifScanX use user-provided installations under `external/`.

---

## 2. Load structures

Loading coordinates does **not** retrieve annotations.

```text
rmv_fetch 1S72
rmv_fetch 1FFK
rmv_fetch 1S72, 1FFK, 3CC2, 1J5A, 1RLG
```

Local files are accepted too:

```text
rmv_fetch /absolute/path/to/structure.cif
rmv_fetch /absolute/path/to/structure.pdb
```

---

## 3. Load annotation sources

Select one or more named sources; names are case-insensitive.

```text
rmv_db RNA3DMotifAtlas
rmv_db RNA3DMotifAtlas, Rfam
rmv_db FR3D, RNAMotifScanX
```

`rmv_db` applies the selected sources to every structure already fetched. Use
`rmv_refresh` to bypass the cache and re-fetch.

---

## 4. Query motifs

The selection grammar has four comma-separated clauses:

```text
rmv_select <motif>, <structures>, <sources>, as <group>
```

Examples:

```text
rmv_select SR, 1S72, RNA3DMotifAtlas, as group_SR
rmv_select SR, 1S72 and 1FFK, RNA3DMotifAtlas and Rfam, as group_shared
rmv_select KT, all, RNA3DMotifAtlas and Rfam, as group_KT
rmv_select SR, all, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP
rmv_select SR, all, RNA3DMotifAtlas and not RNAMotifScanX, as group_FN
```

Notes:

- **Motif names** accept abbreviations, separators, and family-index suffixes.
  `SR`, `SARCIN-RICIN`, `sarcin_ricin`, `sarcin-ricin-1`, and `sarcin-ricin-2`
  all resolve to the same family. Distinct families such as `K-TURN` and
  `REVERSE-K-TURN` stay separate.
- **Structures** can be one ID, an `and`/`or` list, or `all`.
- **Sources** form a Boolean expression with `not` > `and` > `or` precedence.
  A source `and` requires both sources on the **same** consolidated row.
- Rows are consolidated by residue-set overlap (Jaccard ≥ 0.60) **or**
  containment (overlap coefficient ≥ 0.80), so a tight motif core and an
  extended annotation of the same motif are treated as one shared row.

Each `rmv_select` saves a **group**: a snapshot of the matching stable motif IDs
plus the original query text.

---

## 5. List and view

List every consolidated row, one motif ID, or a group:

```text
rmv_list
rmv_list group_SR
rmv_list 1S72_00001
```

Highlight residues on the parent structure without creating objects:

```text
rmv_view 1S72_00001
rmv_view group_SR
rmv_view group_SR, color=red
rmv_view group_SR, padding=2
```

Remove a highlight:

```text
rmv_hide group_SR
rmv_hide all
```

Viewing a group does not reorder or resize the saved result.

---

## 6. Create objects and superimpose

Create selectable objects (the camera does not move):

```text
rmv_create_object group_SR
```

Run medoid-based superimposition:

```text
rmv_super group_SR
```

RSMViewer measures pairwise RMSD on temporary copies, reports the
minimum-average-RMSD medoid, and applies the final transformations to the
selected objects. Parent structures are not replaced by motif fragments.

---

## 7. Combine and color groups

Combine saved groups by motif ID, then color them:

```text
rmv_combine group_FP, group_known, as group_combined
rmv_set_color group_FP, red
rmv_set_color group_known, blue
```

The `rmv_set_color` group name follows the command with **no comma** after the
command word.

---

## 8. Export motifs

Export minimal, coordinates-only mmCIF files for a motif ID or group:

```text
rmv_save group_SR cif
rmv_save 1S72_00001 cif
```

Original on-disk coordinates are used where available; unrelated
parent-structure metadata is not copied into the fragments.

---

## 9. The six example applications

### Application 1 — individual motif instance

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR
rmv_list group_SR
rmv_view 1S72_00001
```

### Application 2 — all instances in a family

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR
rmv_view group_SR
```

### Application 3 — multiple families across multiple structures

```text
rmv_fetch 1S72, 1FFK
rmv_db RNA3DMotifAtlas, Rfam
rmv_select SARCIN-RICIN, 1S72 and 1FFK, RNA3DMotifAtlas and Rfam, as group_SR
rmv_select KT, all, RNA3DMotifAtlas and Rfam, as group_KT
rmv_select CL, all, RNA3DMotifAtlas and Rfam, as group_CL
rmv_select EL, all, RNA3DMotifAtlas and Rfam, as group_EL
rmv_view group_SR, group_KT, group_CL, group_EL
```

### Application 4 — superimpose a family

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_select SR, 1S72, RNA3DMotifAtlas, as group_SR
rmv_create_object group_SR
rmv_super group_SR
```

### Application 5 — benchmark a search tool

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas, RNAMotifScanX
rmv_select SR, 1S72, RNA3DMotifAtlas and RNAMotifScanX, as group_TP
rmv_select SR, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP
rmv_select SR, 1S72, RNA3DMotifAtlas and not RNAMotifScanX, as group_FN
```

`group_TP` = predictions supported by the reference; `group_FP` = predictions
without support; `group_FN` = reference motifs the tool missed.

### Application 6 — inspect false positives

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas, RNAMotifScanX
rmv_select SR, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP
rmv_select SR, 1S72, RNA3DMotifAtlas, as group_known
rmv_set_color group_FP, red
rmv_set_color group_known, blue
rmv_combine group_FP, group_known, as group_combined
rmv_create_object group_combined
rmv_super group_combined
```

---

## 10. FR3D pipeline

### Setup

1. Paste the official fr3d-python software into
   `external/fr3d/fr3d-python-latest/` (must contain `fr3d/__init__.py` and
   `fr3d/search/FR3D.py`).
2. Install FR3D's Python dependencies once:

   ```text
   rmv_fr3d setup
   ```

   This installs `numpy`, `scipy`, and `mmcif-pdbx`. Alternatively set
   `python_path` in `config/fr3d_config.json` to an interpreter that has them.
3. Check readiness:

   ```text
   rmv_fr3d status
   ```

### Run

```text
rmv_fetch 1S72
rmv_db FR3D
```

RSMViewer runs FR3D's own default queries from the checkout against the loaded
structure and loads the resulting candidates. FR3D's geometric queries define
their template from a reference PDB, so set `allow_network: true` in
`config/fr3d_config.json` to let FR3D download that reference. Do not pass an
FR3D result file to `rmv_fetch`.

---

## 11. RNAMotifScanX pipeline

RMSX has two modes, set by `data_mode` in `config/rmsx_config.json`.

### Preannotated mode (default)

Place the preannotated bundle at:

```text
external/rmsx_preannotated/rmsx_preannotated_input_output.tar.gz
```

or an extracted directory at:

```text
external/rmsx_preannotated/rmsx_work_default/<pdb_id>/<chain>/
```

Example of a single output file:

```text
external/rmsx_preannotated/rmsx_work_default/1s72/0/sarcin-ricin_consensus.log
```

Then:

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas, RNAMotifScanX
```

The bundle contains both the RMSX **inputs** (`.rmsx.in` / `.rmsx.nch`) and the
precomputed **outputs** (`*_consensus.log`). In preannotated mode RSMViewer
reads the outputs directly.

### From-scratch mode

Set:

```json
"data_mode": "run_from_scratch"
```

and place binaries under `external/rmsx/bin/` (`scan`, `MC-Annotate`, optional
`rnaview`). Even then, RSMViewer first reuses prebuilt inputs (`.rmsx.in`) for
the requested PDB from the preannotated archive/directory if present, skipping
MC-Annotate; it only regenerates inputs from scratch when none are found. It
then runs the RMSX `scan` step. If `data_mode` is left as `preannotated`, the
precomputed outputs from the bundle are shown.

P-value thresholds are configured only under `pvalue_thresholds` in
`config/rmsx_config.json`; the same thresholds apply to preannotated and freshly
generated logs. A result may legitimately contain zero accepted motifs when
every P-value exceeds its threshold.

See [config/README.md](../config/README.md) for the full field reference.

---

## 12. Diagnostics and validation

```text
rmv_db               # list sources and usage (no argument)
rmv_help             # command reference
rmv_refresh          # bypass caches
rmv_debug ON         # verbose diagnostics (off by default)
rmv_reset            # clear visual/session/cache state
```

Automated checks from the project root:

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q rsmviewer tests
RSMVIEWER_ROOT="$PWD" pymol -cq tests/pymol_applications_e2e.py
RSMVIEWER_ROOT="$PWD" pymol -cq tests/pymol_fetch_multi_smoke.py
```

`tests/pymol_applications_e2e.py` runs all six applications and prints a
PASS/FAIL report.

---

## 13. Common problems

**`unknown color '1FFK'` during multiple fetch.** Use the comma syntax shown
above; the earlier parser bug is fixed.

**A source `and` query returns nothing.** The two sources may annotate the motif
at residue sets that neither overlap by ≥ 0.60 Jaccard nor by ≥ 0.80 containment.
Inspect the rows with `rmv_list`.

**RMSX returns zero motifs.** Compare the log P-values against
`config/rmsx_config.json`. This is usually correct filtering, not a load
failure.

**FR3D is not found.** Verify the checkout path in `config/fr3d_config.json`,
run `rmv_fr3d setup`, and check `rmv_fr3d status`.
