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

1. Download or clone the RSMViewer repository and extract it if necessary.
2. Open PyMOL, then **Plugin → Plugin Manager → Settings → Add New Directory**.
3. Select the `rsmviewer` folder inside the extracted `RSMViewer-main` directory
   (the path should end with `RSMViewer-main/rsmviewer`; select the whole folder,
   not `__init__.py`) and restart PyMOL. RSMViewer loads automatically.
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

List every consolidated row, one motif ID, a group, or a whole family:

```text
rmv_list
rmv_list group_SR
rmv_list 1S72_00001
rmv_list SARCIN-RICIN        # every row any source labels as this family
```

Highlight residues on the parent structure without creating objects. The whole
structure (all atoms, including protein and ligands) is set to gray80 first so
the highlighted motif residues stand out:

```text
rmv_view 1S72_00001
rmv_view group_SR
rmv_view group_SR, color=red
rmv_view group_SR, padding=2
```

Several targets can be highlighted together; each keeps its own color instead of
overwriting the previous one:

```text
rmv_view group_SR, group_KT, group_CL
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
selected objects. The superimposed instances inherit a single color — the
group's session color if you set one (`rmv_set_color` / `rmv_color`), otherwise
the default motif-family color — rather than a different color per instance. The
view auto-orients and zooms onto the superimposed section, and only the
superimposed objects stay visible. Parent structures are not replaced by motif
fragments.

---

## 7. Combine and color groups

Combine saved groups by motif ID, then color them:

```text
rmv_combine_groups group_FP, group_known, as group_combined
rmv_set_color group_FP, red
rmv_set_color group_known, blue
```

The `rmv_set_color` group name follows the command with **no comma** after the
command word.

When you set a color on each source group and then combine them, the per-source
colors are preserved for each member. After `rmv_create_object group_combined`,
the `group_FP` objects remain red and the `group_known` objects remain blue. An
explicit color set on the combined group itself instead colors every member
uniformly.

The combined group keeps each database's own columns and original labels (for
example separate `RNAMotifScanX` and `FR3D` columns), not the input group names.
When overlapping instances merge, the surviving row retains every contributing
database's label, including differing family assignments. The original input
groups and raw annotations are left unchanged.

---

## 8. Export motifs

Export minimal, coordinates-only mmCIF files for a motif ID, a group, or all
motifs:

```text
rmv_save group_SR cif
rmv_save 1S72_00001 cif
rmv_save ALL cif
```

Save a high-resolution PNG image of the current view:

```text
rmv_save current
rmv_save current group_SR.png
```

Every save prints its output directory. mmCIF files are written under
`motif_structures/<pdb_id>/`, motif images under `motif_images/<pdb_id>/`, and a
current-view PNG to the filename you provide. Original on-disk coordinates are
used for mmCIF export where available; unrelated parent-structure metadata is
not copied into the fragments.

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
rmv_combine_groups group_FP, group_known, as group_combined
rmv_create_object group_combined
rmv_super group_combined
```

---

## 10. FR3D pipeline

The repository default is `data_mode: "run_from_scratch"` in
`config/fr3d_config.json`: `rmv_db FR3D` executes the official FR3D pipeline
after the setup below. With `data_mode: "cache"`, `rmv_db FR3D` instead loads
cached annotations from an earlier run of the structure from `output/fr3d_runs/`
and does not run FR3D or use the network.

### Setup

1. Paste the official fr3d-python software into
   `external/fr3d/fr3d-python-latest/` (must contain `fr3d/__init__.py` and
   `fr3d/search/FR3D.py`).
2. Run the one-shot setup (finds a Python, installs deps, registers FR3D):

   ```text
   rmv_setup FR3D
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

In `run_from_scratch` mode, RSMViewer runs FR3D's own default queries from the
checkout against the loaded structure and loads the resulting candidates.
FR3D's geometric queries define
their template from a reference PDB, so set `allow_network: true` in
`config/fr3d_config.json` to let FR3D download that reference. Do not pass an
FR3D result file to `rmv_fetch`.

---

## 11. RNAMotifScanX pipeline

RMSX has two modes, set by `data_mode` in `config/rmsx_config.json`.

### Preannotated mode (default)

To give you the most up-to-date RNAMotifScanX annotations, RSMViewer collects the
preannotated data live from our server. You do not download or extract anything
yourself: the first time you request a PDB, RSMViewer downloads
`<preannotated_base_url>/<pdb_lowercase>.tar.gz` (for example
`https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx_work_default/1s72.tar.gz`)
and extracts it to:

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

Later loads of the same PDB use the local copy and need no internet connection
(delete the `<pdb_id>/` folder to fetch the newest published version again).
For each PDB the lookup order is: local folder, then download. The download link
uses the lowercase PDB ID (`<preannotated_base_url>/<pdb>.tar.gz`) and is tried
only for a 4-character PDB ID with no local results; the archive must extract
safely and contain at least one `*_consensus.log`. If the server has no results
for a PDB yet, RSMViewer says so and you can scan it with `run_from_scratch`.

### From-scratch mode

Set `data_mode` in `config/rmsx_config.json` to:

```json
"data_mode": "run_from_scratch"
```

and use the same command as before:

```text
rmv_fetch 1S72
rmv_db RNAMotifScanX
```

**The first time you run it, RSMViewer downloads all the required files and
third-party software by itself** from
[https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz](https://cbb.ittc.ku.edu/RNAMotifScanX_Results/RSMViewer/rmsx.tar.gz)
into `external/rmsx/` and prepares them for your computer, then runs the analysis
and loads the result. PyMOL waits until it finishes (about a minute and a half
for 1S72 after the first-time preparation). On macOS you need the Xcode command
line tools (`xcode-select --install`); on Windows, WSL2 or Docker Desktop; Linux
needs nothing extra. `rmv_rmsx_doctor` shows whether everything is ready. See
[external/rmsx_setup.md](../external/rmsx_setup.md) for the full guide.

P-value thresholds are configured only under `pvalue_thresholds` in
`config/rmsx_config.json`; the same thresholds apply to preannotated and freshly
generated logs, and to single- and multi-source `rmv_db`. The file is re-read on
every `rmv_db`, so an edit takes effect the next time you run it. A result may
legitimately contain zero accepted motifs when every P-value exceeds its
threshold.

See [config/README.md](../config/README.md) for the full field reference.

---

## 12. Diagnostics and validation

```text
rmv_db               # list sources and usage (no argument)
rmv_help             # command reference
rmv_refresh          # bypass caches
rmv_debug ON         # verbose diagnostics (off by default)
rmv_reset            # show details about rmv_reset cache/session (no reset)
rmv_reset cache      # clear only the caches
rmv_reset session    # clear only objects and session state
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

**RMSX says no results are available for a PDB.** The preannotated results are
downloaded from our server, which may not have every PDB yet, or the download
could not reach it (check your internet connection). RSMViewer prints which
case it is and the URL it tried; use `run_from_scratch` to scan the structure
yourself.

**FR3D is not found.** Verify the checkout path in `config/fr3d_config.json`,
run `rmv_setup FR3D`, and check `rmv_fr3d status`.
