# User Annotations Guide: RNAMotifScanX (Source 7) & NoBIAS (Source 8)

This guide explains the folder/file structure for **RNAMotifScanX (RMSX)** and **NoBIAS** user annotation sources, how to set up custom data paths, and the full workflow from loading a structure to visualizing motifs.

---

## Table of Contents

1. [File Format](#file-format)
2. [RNAMotifScanX (Source 7) — Folder Structure](#rnamotifscanx-source-7--folder-structure)
3. [NoBIAS (Source 8) — Folder Structure](#nobias-source-8--folder-structure)
4. [Setting Up a Custom Data Path](#setting-up-a-custom-data-path)
5. [P-Value Thresholds](#p-value-thresholds)
6. [Workflow: RNAMotifScanX (Source 7)](#workflow-rnamotifscanx-source-7)
7. [Workflow: NoBIAS (Source 8)](#workflow-nobias-source-8)

---

## File Format

Both RMSX and NoBIAS use the **same tab-separated format** with a header line:

```
#fragment_ID	aligned_regions	alignment_score	P-value
1S72_0:75-85_89-98_58-60	0:'0'77-4:'0'81,13:'0'93-20:'0'100	167	0.00639525
1S72_0:933-942_1021-1034	1:'0'936-5:'0'940,13:'0'1026-21:'0'1034	148.4	0.00852698
```

| Column | Description |
|--------|-------------|
| `#fragment_ID` | `PDB_CHAIN:range1_range2_range3` — identifies the fragment in the structure |
| `aligned_regions` | `motif_idx:'chain'residue-motif_idx:'chain'residue` pairs (comma-separated) |
| `alignment_score` | Numeric alignment score (higher = better match) |
| `P-value` | Statistical significance (lower = more significant) |

> The header line starting with `#` is **required**. Lines starting with `#` or `No base-stacking` are skipped during parsing.

---

## RNAMotifScanX (Source 7) — Folder Structure

RMSX organizes results into **motif-type subfolders**, each containing multiple result files:

```
RNAMotifScanX/
├── c-loop_consensus/
│   ├── result_0_100_withbs.log    ← HIGHEST PRIORITY (used for parsing)
│   ├── result_0_100.log           ← 2nd priority
│   ├── result_0_withbs.log        ← 3rd priority
│   ├── result_0.log               ← 4th priority
│   ├── result.log                 ← (not auto-selected, but valid)
│   └── result_9.log               ← (not auto-selected, but valid)
├── k-turn_consensus/
│   ├── result_0_100_withbs.log
│   ├── result_0_100.log
│   ├── result_0_withbs.log
│   ├── result_0.log
│   ├── result.log
│   └── result_9.log
├── sarcin-ricin_consensus/
│   └── ... (same structure)
├── e-loop_consensus/
│   └── ... (same structure)
└── reverse-kturn_consensus/
    └── ... (same structure)
```

### Key Rules

1. **Folder naming**: `<motif-name>_consensus` — the `_consensus` suffix is automatically stripped. The folder name determines the motif type:
   - `k-turn_consensus` → **K-TURN**
   - `c-loop_consensus` → **C-LOOP**
   - `sarcin-ricin_consensus` → **SARCIN-RICIN**
   - `e-loop_consensus` → **E-LOOP**
   - `reverse-kturn_consensus` → **REVERSE-K-TURN**

2. **File selection priority**: The parser picks **one file per folder** using this priority order:
   | Priority | Filename | Description |
   |----------|----------|-------------|
   | 1st | `result_0_100_withbs.log` | Top 100 results with base-stacking info |
   | 2nd | `result_0_100.log` | Top 100 results without base-stacking |
   | 3rd | `result_0_withbs.log` | Default results with base-stacking |
   | 4th | `result_0.log` | Default results |

   If none of these exist, the parser falls back to any `result_*.log` file in the folder.

3. **PDB matching**: The parser reads the file content to find lines containing the target PDB ID (e.g., `1S72`). Only matching lines are loaded.

4. **Minimum required**: At least **one motif subfolder** with at least **one `result_*.log` file** containing data for your PDB.

---

## NoBIAS (Source 8) — Folder Structure

NoBIAS uses a **flat file structure** — no subfolders needed:

```
NoBIAS/
├── 1s72_k-turn_nobias.txt
├── 1s72_c-loop_nobias.txt
└── 1s72_sarcin-ricin_nobias.txt
```

### Key Rules

1. **File naming convention**: `<pdb_id>_<motif-name>_nobias.txt`
   - `1s72_k-turn_nobias.txt` → PDB **1S72**, motif type **K-TURN**
   - `1s72_c-loop_nobias.txt` → PDB **1S72**, motif type **C-LOOP**
   - `1s72_sarcin-ricin_nobias.txt` → PDB **1S72**, motif type **SARCIN-RICIN**

2. **PDB ID**: First 4 characters of the filename (before the first `_`). Case-insensitive.

3. **Motif name**: Extracted from between the PDB ID and `_nobias` suffix. Automatically normalized using the same rules as RMSX (e.g., `k-turn` → `K-TURN`).

4. **Minimum required**: At least **one `*_nobias.txt` file** matching the target PDB ID.

---

## Setting Up a Custom Data Path

By default, the plugin looks for data in the bundled directories:
- RMSX: `rsmviewer/database/user_annotations/RNAMotifScanX/`
- NoBIAS: `rsmviewer/database/user_annotations/NoBIAS/`

You can point to **your own data directory** using a path argument:

```
rmv_db 7 /path/to/my/rmsx_data
rmv_db 8 /path/to/my/nobias_data
```

### Custom Path for RMSX (Source 7)

Your custom directory must follow the **same subfolder structure**:

```
/path/to/my/rmsx_data/
├── k-turn_consensus/
│   └── result_0_100_withbs.log    ← must have at least one result_*.log
├── c-loop_consensus/
│   └── result_0_100_withbs.log
└── sarcin-ricin_consensus/
    └── result_0_100.log
```

**Flat file fallback**: If no subfolders are found, the parser will try to load any `.txt` / `.log` / `.tsv` file directly in the directory that contains the target PDB ID. This lets you place RMSX output files directly:

```
/path/to/my/rmsx_data/
├── 1S72_0_kturn.txt
└── 1S72_0_cloop.txt
```

### Custom Path for NoBIAS (Source 8)

Your custom directory just needs flat files following the naming convention:

```
/path/to/my/nobias_data/
├── 1s72_k-turn_nobias.txt
├── 1s72_c-loop_nobias.txt
└── 4v9f_k-turn_nobias.txt
```

**Flat file fallback**: Same as RMSX — if no `*_nobias.txt` files match, any `.txt` / `.log` / `.tsv` file containing the PDB ID will be loaded.

### Single-File Mode

You can also point directly to a single file:

```
rmv_db 7 /path/to/my/result_0_100_withbs.log
rmv_db 8 /path/to/my/1s72_k-turn_nobias.txt
```

---

## P-Value Thresholds

Both sources filter results by motif-specific P-value thresholds. Instances with P-value **above** the threshold are discarded. Default thresholds (same for both RMSX and NoBIAS):

| Motif | Default P-value Threshold |
|-------|--------------------------|
| K-TURN (KINK-TURN) | 0.066 |
| C-LOOP | 0.044 |
| SARCIN-RICIN | 0.040 |
| REVERSE KINK-TURN | 0.018 |
| E-LOOP | 0.018 |
| Other motifs | 0.05 (fallback) |

### Controlling P-value Filtering

```
rmv_db 7 off                         # Disable filtering — show ALL raw results
rmv_db 7 on                          # Re-enable default filtering
rmv_db 7 K-TURN 0.01                 # Custom threshold for K-TURN only
rmv_db 7 K-TURN 0.01 C-LOOP 0.02    # Custom thresholds for multiple motifs

rmv_db 8 off                         # Same commands work for NoBIAS
rmv_db 8 on
rmv_db 8 K-TURN 0.01
```

---

## Workflow: RNAMotifScanX (Source 7)

### Using Bundled Data

```python
# Step 1: Load the PDB structure
rmv_fetch 1S72

# Step 2: Select RMSX as data source
rmv_db 7

# Step 3: (Optional) Turn off P-value filtering to see all raw hits
rmv_db 7 off

# Step 4: Load motif data
rmv_load_motif

# Step 5: View summary of all motif types and counts
rmv_summary

# Step 6: View instances for a specific motif type
rmv_summary K-TURN

# Step 7: Show all instances of a motif type (clustered visualization)
rmv_show K-TURN

# Step 8: Show a specific instance by its number from the summary table
rmv_show K-TURN 1
```

### Using a Custom Data Path

```python
# Step 1: Load the PDB structure
rmv_fetch 1S72

# Step 2: Select RMSX with custom data path
rmv_db 7 /path/to/my/rmsx_results

# Step 3: (Optional) Set custom P-value thresholds
rmv_db 7 K-TURN 0.03 SARCIN-RICIN 0.01

# Step 4: Load motif data
rmv_load_motif

# Step 5-8: Same as above
rmv_summary
rmv_summary K-TURN
rmv_show K-TURN
rmv_show K-TURN 1
```

### Toggling P-values Mid-Session

```python
# Already loaded with filtering ON — now want to see everything:
rmv_db 7 off
rmv_load_motif        # Re-load to apply the change

# Back to default thresholds:
rmv_db 7 on
rmv_load_motif
```

---

## Workflow: NoBIAS (Source 8)

### Using Bundled Data

```python
# Step 1: Load the PDB structure
rmv_fetch 1S72

# Step 2: Select NoBIAS as data source
rmv_db 8

# Step 3: (Optional) Turn off P-value filtering
rmv_db 8 off

# Step 4: Load motif data
rmv_load_motif

# Step 5: View summary of all motif types and counts
rmv_summary

# Step 6: View instances for a specific motif type
rmv_summary K-TURN

# Step 7: Show all instances of a motif type (clustered visualization)
rmv_show K-TURN

# Step 8: Show a specific instance by its number from the summary table
rmv_show K-TURN 1
```

### Using a Custom Data Path

```python
# Step 1: Load the PDB structure
rmv_fetch 4V9F

# Step 2: Select NoBIAS with custom data path
rmv_db 8 /path/to/my/nobias_results

# Step 3: (Optional) Set custom P-value thresholds
rmv_db 8 K-TURN 0.02 C-LOOP 0.01

# Step 4: Load motif data
rmv_load_motif

# Step 5-8: Same as above
rmv_summary
rmv_summary K-TURN
rmv_show K-TURN
rmv_show K-TURN 1
```

### Toggling P-values Mid-Session

```python
# Already loaded with filtering ON — now want to see everything:
rmv_db 8 off
rmv_load_motif        # Re-load to apply the change

# Back to default thresholds:
rmv_db 8 on
rmv_load_motif
```

---

## Quick Reference

| Action | RMSX (Source 7) | NoBIAS (Source 8) |
|--------|----------------|-------------------|
| Select source | `rmv_db 7` | `rmv_db 8` |
| Custom path | `rmv_db 7 /path/to/data` | `rmv_db 8 /path/to/data` |
| Filtering OFF | `rmv_db 7 off` | `rmv_db 8 off` |
| Filtering ON | `rmv_db 7 on` | `rmv_db 8 on` |
| Custom P-values | `rmv_db 7 MOTIF 0.05` | `rmv_db 8 MOTIF 0.05` |
| Load data | `rmv_load_motif` | `rmv_load_motif` |
| Summary | `rmv_summary` | `rmv_summary` |
| Type details | `rmv_summary K-TURN` | `rmv_summary K-TURN` |
| Show all of type | `rmv_show K-TURN` | `rmv_show K-TURN` |
| Show one instance | `rmv_show K-TURN 1` | `rmv_show K-TURN 1` |
