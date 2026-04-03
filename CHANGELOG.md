# RSMViewer — Change Log

All notable updates to the RSMViewer project are recorded here.
Entries are listed newest-first.

---

## 2026-04-02 — Padding Refresh Fixes, Combine Suffix Naming, and Jaccard Config

### Bug Fixes: `rmv_show` Padding + Object Lifecycle
1. **Padding object refresh now updates correctly** (`loader.py` `show_motif_instance()`)
  - Re-running the same instance with a new padding value now deletes and recreates the target object
  - Fixes stale geometry issue where `padding=10` persisted after `padding=2`

2. **Conflicting object variant cleanup (base vs padded)** (`loader.py` `show_motif_instance()`)
  - Showing a padded instance now deletes the non-padded variant (and vice versa)
  - Prevents double-render overlap artifacts from simultaneously visible objects

3. **Combine-mode suffix naming fixed for instance/type rendering** (`loader.py`)
  - `show_motif_instance()` and `show_motif_type()` now prefer `filter_suffix` for object naming when provided
  - Ensures combine objects use correct suffix format from active mode (e.g., `_S_8_7`)
  - Fixes incorrect names like `K_TURN_3_1S72_S7` in combine mode

### Feature: User-Configurable Cascade Jaccard Threshold
1. **Optional `jaccard_threshold` in `rmv_db`** (`gui.py`)
  - New command kwarg: `rmv_db 8 7, jaccard_threshold=0.80`
  - Supported formats: decimal (`0.80`), integer percent (`80`), and percent string (`80%`)
  - Values > 1 are treated as percentages and normalized to 0-1
  - Validation added for valid range `(0, 1]`

2. **Threshold wiring into merge pipeline** (`gui.py`)
  - Added GUI state: `self.jaccard_threshold = 0.60` default
  - `_load_combined_motifs()` now passes `self.jaccard_threshold` into `CascadeMerger(...)`
  - Multi-source status output now reports active threshold in percent

### Documentation
1. **`COMBINE_GUIDE.md` updated for threshold configurability**
  - Added optional Jaccard threshold usage in Overview
  - Rewrote Rule 5 as user-configurable (default + override behavior)
  - Updated FAQ with accepted formats, persistence behavior, and examples

### Validation
- Syntax checks passed for modified Python files:
  - `python -m py_compile rsmviewer/gui.py`
  - `python -m py_compile rsmviewer/loader.py`

---

## 2026-04-01 — Combine Mode Refinements & rmv_view Feature

### Bug Fixes: Combine Mode
1. **Stale categories sweep** (`gui.py` `fetch_motif_data_action()`)
   - Before accumulating new motif data, now performs full cleanup of all instances matching `(pdb_id, source_suffix)` across ALL motif types in `existing_loaded`
   - Fixes issue where old motif types would persist in summary when re-loading with fewer types

2. **Object name suffix format** (`gui.py` `_get_source_suffix()`)
   - Combine mode now correctly formats suffix as `_S_{source_ids}` (e.g., `_S_8_7`) instead of `_S{source_ids}` (e.g., `_S8_7`)
   - Example: `K-TURN_ALL_1S72_S_8_7` (was `K-TURN_ALL_1S72_S8_7`)
   - Single-source mode unchanged: `_S7`

3. **Source attribution with shared instances** (`cascade_merger.py`, `loader.py`, `gui.py`)
   - **cascade_merger.py** `_annotate_overlap()`: new method tracks instances discarded during Jaccard overlap detection
     - Adds `_also_found_in` metadata list to surviving ref instances
     - Propagates overlaps for 3+ source cascades
   - **cascade_merger.py** `_pairwise_merge()`: tracks best overlapping instance and calls `_annotate_overlap()` when discarding
   - **loader.py** `_format_source_label()`: new helper method combines `_source_label` + `_also_found_in` for display
     - Example: "NoBIAS + RMSX" for instances with cross-source support
   - **loader.py** `_print_motif_instance_table()`: uses `_format_source_label()` to display source info in instance tables (width: 30 chars)
   - **gui.py** `_print_source_attribution_report()`: three categories of instances:
     - "Unique in NoBIAS"
     - "Unique in RMSX"
     - "Shared (NoBIAS + RMSX)" for overlapping instances

4. **Padding support** (`gui.py`, `loader.py`)
   - **gui.py** `show_motif()`: now checks `**_kwargs['padding']` after positional args (handles PyMOL comma syntax)
   - **loader.py** `show_motif_type()`: when `padding > 0`, creates `_P`-suffixed object copy (e.g., `K-TURN_ALL_1S72_S_8_7_P`)
   - **loader.py** `show_motif_instance()`: appends `_P` to instance object names when padded (e.g., `K-TURN_1_1S72_S_8_7_P`)

### New Feature: `rmv_view` Command
- **Purpose**: zoom to motif regions and create selections (no objects created)
- **`rmv_view motif`**: highlights ALL loaded motif regions on the base PDB structure with unique colors per type (gray80 base + colored residues, no objects)
- **gui.py** `_auto_color_motifs_on_structure()`: method that performs the coloring (called by `rmv_view motif`, NOT by `rmv_load_motif`)
- **loader.py** `view_motif_type()`
  - Zooms to all instances of a motif type
  - Colors only the matching residues; rest of structure in gray80
  - Prints instance table
- **loader.py** `view_motif_instance()`
  - Zooms to specific instance
  - Creates PyMOL *selection* (not object): `sele_<TYPE>_<NO>` (e.g., `sele_K_TURN_1`)
  - Prints instance details
- **gui.py** `view_motif()`: command function with `motif` keyword + instance parsing
- **Command registration**: `cmd.extend('rmv_view', view_motif)`
- **Help**: added `rmv_view motif`, `rmv_view <TYPE>`, `rmv_view <TYPE> <NO>` to `print_help()`

### Documentation
- Added `rmv_view` usage to `print_help()` Next Steps suggestions
- Help command includes: `rmv_view <TYPE>` and `rmv_view <TYPE> <NO>` examples

### Validation
- All modifie files syntax-checked: gui.py, loader.py, cascade_merger.py ✓

---

## 2026-03-29 — Medoid-Based Structural Superimposition

### New Feature: `rmv_super` / `rmv_align`
- **Rewritten file**: `rsmviewer/alignment.py` — full medoid superimposition pipeline
  - `rmv_super <TYPE>` — sequence-independent superimposition (`cmd.super`)
  - `rmv_align <TYPE>` — sequence-dependent superimposition (`cmd.align`)
  - Automatically identifies the **medoid** (instance with minimum average pairwise RMSD)
  - Superimposes all instances onto the medoid
  - Each instance gets a unique colour; medoid highlighted in green
  - Prints formatted report with per-instance RMSD, overall average, and skipped pairs
- **Instance filtering**: `rmv_super KTURN 1,3,5` for specific instances
- **Cross-PDB support**: `rmv_super KTURN pdb=1S72,4V88` — any number of PDBs
- **Motif name aliasing**: `MOTIF_ALIASES` table resolves common variants (KTURN→K-TURN, CLOOP→C-LOOP, etc.) with substring fallback
- **Batch object creation**: suppresses PyMOL GUI updates during bulk creation for performance
- **CSV export**: optional `save_matrix=/path.csv` to export pairwise RMSD matrix

### GUI Updates
- Added `rmv_super`, `rmv_align` to `_RMV_COMMANDS` and `_CMD_SUFFIXES` (typo suggestion system)
- Added `register_alignment_commands()` call in `initialize_gui()`
- Added superimposition example (#7) to `print_help()` Quick Examples
- Added `rmv_super KTURN` hint to Quick Start banner

### Documentation
- **README.md**: added "Structural Superimposition" command table
- **TUTORIAL.md**: added section 9 "Structural Superimposition (Medoid)" with single-PDB, cross-PDB, and alias examples; renumbered sections 10–11
- **MEDOID_SUPERIMPOSITION_PLAN.md**: updated examples section

---

## 2026-03-27 — Major Update: Rename + Alignment

### Project Rename
- Renamed entire package folder: `rna_motif_visualizer/` → `rsmviewer/`
- Updated all module docstrings (26+ files) from "RNA Motif Visualizer" to "RSMViewer"
- Updated display strings: plugin banner (`🧬 RSMViewer 🧬`), GUI help header, status header, colour legend title
- Updated cache directory: `~/.rna_motif_visualizer_cache/` → `~/.rsmviewer_cache/`
- Updated logger names to `rsmviewer.*`
- Updated API User-Agent headers to `RSMViewer/*`
- Updated documentation files: README.md, DEVELOPER.md, TUTORIAL.md
- Cleared all `__pycache__/` directories

### GUI & Help Updates
- Corrected recommended workflow across all files to include `rmv_sources` step
  and full progression: `rmv_fetch` → `rmv_sources` → `rmv_db` → `rmv_load_motif`
  → `rmv_summary` → `rmv_summary <TYPE>` → `rmv_show <TYPE>` → `rmv_show <TYPE> <NO>`
- `print_help()`: added STRUCTURAL ALIGNMENT section
- `_RMV_COMMANDS` list updated

---

## 2026-03-26 — Base-Pair Visualizer

### New Feature: `rmv_pair`
- **New file**: `rsmviewer/pair_visualizer.py`
  - `rmv_pair` command: visualize RNA base pairs with Leontis-Westhof edge labels (W/H/S)
  - `rmv_pair_batch` command: batch visualization from file
  - Input format: `pdbid_chain1_resnum1_chain2_resnum2` (residue names auto-detected)
  - Handles same-chain and cross-chain pairs
  - Solvent removal step added
  - Object naming: `pair_<pdbid>_<chain>_<res1>_<res2>` (same-chain) or `pair_<pdbid>_<ch1>_<res1>_<ch2>_<res2>` (cross-chain)

### Bug Fixes
- Fixed `rmv_load` hanging issue — now shows workflow guide instead of blocking
- Fixed motif object naming conflicts — PDB ID now included in object names
- Fixed parser for numeric chain IDs — uses positional 4-token parsing
- Removed backtick syntax from chain selections in pair_visualizer.py
