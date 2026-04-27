# RSMViewer: Tutorial

A complete step-by-step guide for using the RSMViewer PyMOL plugin.

---

## Table of Contents

1. [Getting Started](#1-getting-started)
2. [Basic Workflow](#2-basic-workflow)
3. [Exploring Motifs](#3-exploring-motifs)
4. [Working with Data Sources](#4-working-with-data-sources)
5. [Multi-Source Comparison](#5-multi-source-comparison)
6. [User Annotations (FR3D / RMS / RMSX)](#6-user-annotations-fr3d--rms--rmsx)
7. [Customizing Colors](#7-customizing-colors)
8. [Saving Images](#8-saving-images)
9. [Exporting Motif Structures as mmCIF](#9-exporting-motif-structures-as-mmcif)
10. [Structural Superimposition (Medoid)](#10-structural-superimposition-medoid)
11. [Label Chain IDs (Advanced)](#11-label-chain-ids-advanced)
12. [Tips and Tricks](#12-tips-and-tricks)
13. [Quick Reference Card](#13-quick-reference-card)

---

## 1. Getting Started

### Prerequisites

- **PyMOL** (Incentive or Open-Source, version 2.x or later)
- **Python 3.8+** (bundled with PyMOL)
- Internet connection (for BGSU API and Rfam API sources only; offline sources work without internet)

### Installation

1. Open PyMOL
2. Go to **Plugin → Plugin Manager**
3. Click the **Settings** tab
4. Click **Add new directory**
5. Navigate to and select the `rsmviewer` folder
6. Click **OK** and restart PyMOL

After launching PyMOL, you should see:

```
================================================================================
                    🧬 RSMViewer 🧬

 Version 1.0.0
 Last Updated: 02 April 2026
================================================================================
```

Run `rmv_help` at any time to see all available commands.

---

## 2. Basic Workflow

The plugin uses a **3-step workflow**: **Load → Source → Motifs**.

### Step 1 — Load a PDB Structure

```
rmv_fetch 1S72
```

This downloads the structure from RCSB PDB and loads it into PyMOL:

```
[SUCCESS] Loaded structure 1S72 as '1s72'
[INFO] Chain ID mode: auth_asym_id (default)
[INFO] Chains found: 0, 9, A, B, C, D, ...
```

> **Note:** Only the structure is loaded at this point — no motif data yet. This makes loading fast and lets you pick a data source first.

You can also load local files:

```
rmv_fetch /path/to/local_structure.cif
rmv_fetch ~/structures/1s72.pdb
```

### Step 2 — Select a Data Source

```
rmv_db 3
```

This selects BGSU RNA 3D Hub API — the most comprehensive source with ~3000+ structures. Available sources:

| ID | Source | Type | Coverage |
|----|--------|------|----------|
| 1 | RNA 3D Motif Atlas | Local (offline) | 7 motif types (HL, IL, J3–J7) |
| 2 | Rfam | Local (offline) | 19 named motifs |
| 3 | **BGSU RNA 3D Hub** | Online API | **~3000+ PDBs** |
| 4 | Rfam API | Online API | 34 motif families |
| 5 | FR3D | User annotations | Custom |
| 6 | RNAMotifScan (RMS) | User annotations | Custom (P-value filtering) |
| 7 | RNAMotifScanX (RMSX) | User annotations | Custom (P-value filtering) |

Run `rmv_sources` to see detailed information about each source.

### Step 3 — Fetch Motif Data

```
rmv_load_motif
```

This retrieves motif data from the selected source for your loaded PDB:

```
[SUCCESS] Loaded 41 motif types from 1S72 (source: bgsu_api)
```

### Step 4 — View the Summary

```
rmv_summary
```

Displays a table of all detected motif types and instance counts:

```
==================================================
  MOTIF SUMMARY - 1S72
==================================================
  MOTIF TYPE              INSTANCES
--------------------------------------------------
  HAIRPIN LOOP (HL)              35
  GNRA                           23
  INTERNAL LOOP (IL)             19
  3-WAY JUNCTION (J3)            18
--------------------------------------------------
  TOTAL                         223
==================================================
```

### Step 5 — Visualize Motifs

```
rmv_show HL
```

This creates PyMOL objects, colors the hairpin loop residues, and displays them on the gray structure background.

### The Complete Pipeline

```
rmv_fetch 1S72       # Step 1: Load PDB structure
rmv_sources          # Step 2: Check available data sources
rmv_db 3             # Step 3: Select BGSU API
rmv_load_motif       # Step 4: Fetch motif data
rmv_summary          # Step 5: Show motif types & counts
rmv_summary HL       # Step 6: Show hairpin loop instances
rmv_show HL          # Step 7: Render all hairpin loops
rmv_show HL 1        # Step 8: Zoom to specific instance
rmv_view HL 1        # Step 9: Quick highlight without creating objects
```

---

## 3. Exploring Motifs

### View Instances of a Specific Motif Type

```
rmv_summary SARCIN-RICIN
```

Shows a detailed table of all instances of that motif type with chain, residue ranges, and source attribution:

```
====================================================================
  SARCIN-RICIN MOTIF INSTANCES
====================================================================
  NO.    CHAIN      RESIDUE RANGES
--------------------------------------------------------------------
  1      A          A:158-164, A:171-178
  2      A          A:210-216, A:224-229
  3      A          A:290-296, A:355-362
  ...
--------------------------------------------------------------------
```

### View a Single Instance in Detail

```
rmv_summary SARCIN-RICIN 1
```

Prints full residue-level details for instance #1.

### Render a Motif Type on Structure

```
rmv_show HL              # Show all hairpin loops
rmv_show GNRA            # Show all GNRA tetraloops
rmv_show SARCIN-RICIN    # Show all sarcin-ricin motifs
```

### Render a Specific Instance

```
rmv_show GNRA 1          # Zoom to GNRA instance #1
rmv_show GNRA 1,3,5      # Show multiple specific instances
```

### Quick View (Highlight Without Creating Objects)

`rmv_view` is a lightweight alternative to `rmv_show` — it highlights and zooms to motif regions directly on the structure **without** creating separate PyMOL objects:

```
rmv_view HL              # Highlight all hairpin loop residues
rmv_view HL 1            # Highlight and zoom to HL instance 1
rmv_view GNRA 3          # Highlight and zoom to GNRA instance 3
rmv_view HL hide         # Remove HL coloring from structure
rmv_view hide            # Reset all view coloring to gray
```

> **Tip:** Use `rmv_view` for quick exploration and `rmv_show` when you need persistent PyMOL objects for alignment, saving, or side-by-side comparison. Use `rmv_view hide` to undo any view coloring.

### Show All Motif Types

```
rmv_show ALL             # Create PyMOL objects for all motif types
rmv_view all             # Highlight all motif regions on structure (no objects)
rmv_view all hide        # Reset all coloring to gray (undo rmv_view all)
```

### Padding — Expand Visualization Range

Use `padding=N` to expand the displayed residue range by ±N positions for context. This helps see the structural environment around a motif without changing the summary tables or stored data.

```
rmv_show K-TURN, padding=2       # Show K-TURN with ±2 extra residues
rmv_show K-TURN 1, padding=3     # Show instance 1 with ±3 padding
rmv_show HL 1,2,3, padding=5     # Show instances 1,2,3 with ±5 padding
```

Padding also works with superimposition:

```
rmv_super K-TURN, padding=3      # Superimpose with expanded residue ranges
rmv_align SARCIN-RICIN, padding=2
```

> **Note:** Padding only affects **visualization** (coloring and PyMOL object creation). Summary tables, instance details, and stored data are not affected.

### Toggle Visibility

```
rmv_toggle HL off        # Hide hairpin loop objects
rmv_toggle HL on         # Show them again
```

---

## 4. Working with Data Sources

### Local Sources (Offline — No Internet Needed)

Sources 1 and 2 are bundled with the plugin:

```
rmv_db 1               # RNA 3D Motif Atlas — 7 motif types: HL, IL, J3–J7
rmv_load_motif

rmv_db 2               # Rfam — 19 named motifs: GNRA, UNCG, K-turn, T-loop, C-loop, etc.
rmv_load_motif
```

### Online Sources (API — Internet Required)

Sources 3 and 4 query online databases:

```
rmv_db 3               # BGSU RNA 3D Hub — ~3000+ PDBs (most comprehensive source)
rmv_load_motif

rmv_db 4               # Rfam API — 34 Rfam-annotated motif families
rmv_load_motif
```

> **Caching:** API results are cached for 30 days in `~/.rsmviewer_cache/`. Use `rmv_refresh` to bypass the cache and force a fresh fetch.

### Switch Sources Without Reloading PDB

One of the most useful features — the PDB structure stays loaded while you switch sources:

```
rmv_fetch 1S72         # Load once
rmv_db 3               # BGSU
rmv_load_motif         # Fetch from BGSU

rmv_db 1               # Switch to Atlas
rmv_load_motif         # Fetch from Atlas (no PDB re-download)

rmv_db 2               # Switch to Rfam
rmv_load_motif         # Fetch from Rfam
```

### Check Source Information

```
rmv_sources            # List all 7 sources, their types, and coverage
rmv_source info 3      # Detailed information about BGSU source
rmv_source info        # Show currently active source
```

### Source-Tagged PyMOL Objects

When you switch sources and re-run `rmv_show`, each source gets its own PyMOL objects tagged with a suffix:

```
rmv_db 3
rmv_load_motif
rmv_show SARCIN-RICIN      # Creates: SARCIN_RICIN_ALL_S3, SARCIN_RICIN_1_S3, ...

rmv_db 7
rmv_load_motif
rmv_show SARCIN-RICIN      # Creates: SARCIN_RICIN_ALL_S7, SARCIN_RICIN_1_S7, ...
```

This allows direct comparison of objects from different sources within the same PyMOL session.

---

## 5. Multi-Source Comparison

Combine data from multiple sources. The plugin merges results with smart deduplication.

### Basic Two-Source Combine

```
rmv_fetch 1S72
rmv_db 6 7             # Combine RMS [6] + RMSX [7]  (left = highest priority)
rmv_load_motif         # Fetches from both → deduplicates → merged
rmv_summary            # Unified summary from both sources
```

### Three-Source Combine

```
rmv_db 1 3 4           # Atlas + BGSU + Rfam API
rmv_load_motif
rmv_summary
```

### How Merging Works

1. **Fetch** — motifs are fetched from each source independently
2. **Enrich** — generic names (e.g., "HL") are enriched to specific names (e.g., "GNRA") using NR homolog representative lookup
3. **Source stamp** — each instance is tagged with its source (used for attribution)
4. **Within-source dedup** — exact duplicate entries within the same source are removed (same residue set = duplicate)
5. **Cascade merge** — right-to-left merging with Jaccard deduplication (≥60% residue overlap = duplicate; higher-priority version kept)
6. **Cross-source attribution** — overlapping instances are tagged with `_also_found_in` metadata

### Source Priority

Sources listed first have **highest priority**:

```
rmv_db 6 7             # RMS (higher priority) + RMSX (lower priority)
rmv_db 1 3 4           # Atlas (highest) + BGSU + Rfam API (lowest)
```

### Custom Jaccard Threshold

```
rmv_db 6 7, jaccard_threshold=0.80    # 80% overlap required for deduplication
rmv_db 6 7, jaccard_threshold=40      # Values >1.0 treated as percentages
```

### Source Filtering (Combine Mode)

After combining sources, use **source filter keywords** with `rmv_show` to render only instances from a specific source or those shared across sources:

```
rmv_fetch 1S72
rmv_db 6 ~/my_rms_data                # Set RMS data path
rmv_db 7 /path/to/rmsx/output         # Set RMSX data path
rmv_db 6 7                            # Combine RMS (priority) + RMSX
rmv_load_motif

rmv_summary K-TURN                    # View attribution report (unique/shared)

rmv_show K-TURN rms                   # Only RMS-unique K-TURN instances
rmv_show K-TURN rmsx                  # Only RMSX-unique K-TURN instances
rmv_show K-TURN shared                # Only instances found in both sources
rmv_show K-TURN                       # All K-TURN instances (no filter)
```

Supported keywords: `rmsx`, `rms`, `fr3d`, `rfam`, `atlas`, `bgsu`, `shared`.

> **Tip:** After running `rmv_summary <TYPE>` in combine mode, the Next Steps section automatically suggests source filter commands for the active sources.

---

## 6. User Annotations (FR3D / RMS / RMSX)

Use your own motif annotation files from FR3D, RNAMotifScan (RMS), or RNAMotifScanX (RMSX).

### Directory Setup

Place your annotation files in the appropriate subdirectory (make sure the file format matches):

```
rsmviewer/database/user_annotations/
├── fr3d/                    # FR3D output files (CSV)
│   └── <output_files>.csv
├── RNAMotifScan/            # RMS files (organized by motif type)
│   ├── c_loop/
│   │   └── Res_1s72
│   ├── Kturn/
│   │   └── Res_1s72
│   ├── sarcin_ricin/
│   │   └── Res_1s72
│   └── reverse_Kturn/
│       └── Res_1s72
├── RNAMotifScanX/           # RMSX files (consensus subdirectories)
│   ├── c-loop_consensus/
│   │   └── result_0_100_withbs.log
│   ├── k-turn_consensus/
│   │   └── result_0_100_withbs.log
│   └── sarcin-ricin_consensus/
│       └── result_0_100_withbs.log
```

**RMSX file format** (tab-separated):
```
#fragment_ID	aligned_regions	alignment_score	P-value
1S72_0:75-85_89-98_58-60	0:'0'77-4:'0'81,...	167	0.00639525
```

### Load User Annotations

```
rmv_fetch 1S72

rmv_db 5               # FR3D
rmv_load_motif

rmv_db 6               # RNAMotifScan (RMS)
rmv_load_motif

rmv_db 7               # RNAMotifScanX (RMSX)
rmv_load_motif
```

### Custom Data Paths (Per-Source)

You can specify a custom directory for each annotation source independently. Each source remembers its own path:

```
rmv_db 5 /path/to/fr3d/data       # FR3D with custom directory
rmv_db 6 ~/my_rms_data            # RMS with home-relative path
rmv_db 7 /path/to/rmsx/data       # RMSX with custom directory
```

> **Per-source paths:** Setting a custom path for one source does **not** overwrite the path for another. Each source keeps its own path independently.

### P-Value Filtering (RMS / RMSX)

RNAMotifScan and RNAMotifScanX include P-values in their output. The plugin uses these to filter results.

**Default thresholds:**

| Motif | RMS | RMSX |
|-------|-----|------|
| KINK-TURN | 0.07 | 0.066 |
| C-LOOP | 0.04 | 0.044 |
| SARCIN-RICIN | 0.02 | 0.040 |
| REVERSE KINK-TURN | 0.14 | 0.018 |
| E-LOOP | 0.13 | 0.018 |

**Control filtering:**

```
# Turn filtering off (show ALL results including high P-values)
rmv_db 6 off
rmv_db 7 off

# Turn filtering back on (use default thresholds)
rmv_db 6 on
rmv_db 7 on

# Set custom P-value thresholds for specific motif types
rmv_db 6 SARCIN-RICIN 0.01
rmv_db 7 C-LOOP 0.05
```

When filtering is on, only motif instances with P-value ≤ threshold are shown.

### Legacy Loading

The `rmv_user` command provides an alternative loading path:

```
rmv_user fr3d 1S72           # Load FR3D annotations
rmv_user rnamotifscan 1A00   # Load RMS annotations
rmv_user list                # List available annotation files
```

> **Preferred approach:** Use `rmv_db 5/6/7/8` + `rmv_load_motif` instead.

---

## 7. Customizing Colors

### View Current Color Assignments

```
rmv_colors
```

Shows the color legend for all loaded motif types.

### Change a Motif's Color

```
rmv_color HL blue            # Hairpin loops → blue
rmv_color GNRA red           # GNRA tetraloops → red
rmv_color IL marine          # Internal loops → marine blue
```

### Change Background (Non-Motif) Color

```
rmv_bg_color white           # White backbone
rmv_bg_color lightgray       # Light gray backbone
rmv_bg_color gray80          # Default (gray80)
```

### Default Color Assignments

| Motif Type | Default Color |
|-----------|---------------|
| HL (Hairpin Loop) | Red |
| IL (Internal Loop) | Cyan |
| J3 (3-Way Junction) | Yellow |
| J4 (4-Way Junction) | Magenta |
| J5 (5-Way Junction) | Green |
| J6 (6-Way Junction) | Orange |
| J7 (7-Way Junction) | Blue |
| GNRA | Teal Green |
| UNCG | Brown Orange |
| K-TURN | Bright Blue |
| SARCIN-RICIN | Dark Red |
| T-LOOP | Pink |
| C-LOOP | Sky Blue |
| U-TURN | Gold |

Any motif type not in this table gets a unique auto-assigned color from a pool of 20 visually distinct colors. If the pool is exhausted, a deterministic hash-based color is generated.

---

## 8. Saving Images

### Save All Motif Instances

```
rmv_save ALL                  # Default representation (cartoon)
rmv_save ALL sticks           # Sticks representation
rmv_save ALL spheres          # Spheres representation
```

Images are saved to `motif_images/<PDB_ID>/<MOTIF_TYPE>/` inside the plugin directory.

### Save a Specific Motif Type

```
rmv_save HL                   # All hairpin loop instances as cartoon
rmv_save GNRA sticks          # All GNRA instances as sticks
```

### Save a Specific Instance

```
rmv_save GNRA 1               # Instance 1 as cartoon
rmv_save GNRA 1 spheres       # Instance 1 as spheres
```

### Save the Current PyMOL View

```
rmv_save current              # High-res: 2400×1800, ~300 DPI
rmv_save current my_figure    # Save with custom filename
```

### Available Representations

| Representation | Description |
|---------------|-------------|
| `cartoon` | Default — shows secondary structure |
| `sticks` | Ball-and-stick atomic detail |
| `spheres` | Space-filling model |
| `ribbon` | Simplified ribbon trace |
| `lines` | Wire frame |
| `licorice` | Thick stick model |
| `surface` | Molecular surface |
| `cartoon+sticks` | Combined view |

### Output Folder Structure

```
motif_images/
└── 1s72/
    ├── HL/
    │   ├── HL-1-0_55-64.png
    │   ├── HL-2-0_120-135.png
    │   └── ...
    ├── GNRA/
    │   ├── GNRA-1-0_55-64.png
    │   └── ...
    └── ...
```

---

## 9. Exporting Motif Structures as mmCIF

Extract motif instances as standalone mmCIF files with **original coordinates** from the on-disk CIF file (not PyMOL's internal coordinates, which may be slightly modified during loading).

```
rmv_save ALL cif              # Export ALL motif instances as mmCIF
rmv_save HL cif               # Export all hairpin loop instances as mmCIF
rmv_save HL 3 cif             # Export HL instance #3 as mmCIF
```

Structures are saved to `motif_structures/<PDB_ID>/<MOTIF_TYPE>/` inside the plugin directory. Each `.cif` file is a standalone mmCIF that can be loaded independently:

```
# Load an exported motif in a fresh PyMOL session
cmd.load("motif_structures/1s72/HL/HL-3-0_55-64.cif")
```

### Output Folder Structure

```
motif_structures/
└── 1s72/
    ├── HL/
    │   ├── HL-1-0_55-64.cif
    │   ├── HL-2-0_120-135.cif
    │   └── ...
    └── GNRA/
        ├── GNRA-1-0_55-64.cif
        └── ...
```

---

## 10. Structural Superimposition (Medoid)

`rmv_super` and `rmv_align` automatically find the **medoid** (most representative instance — the one with the lowest average RMSD to all others) and superimpose all other instances onto it.

- **`rmv_super`** — uses `cmd.super()` (sequence-independent). Best for comparing the same motif across different organisms/PDBs or when sequences diverge.
- **`rmv_align`** — uses `cmd.align()` (sequence-dependent). Best when sequences are similar.

### Basic Usage

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif
rmv_super KTURN                 # Superimpose all K-TURN instances
```

Each instance gets a unique color (cycling through 15 distinct colors). The medoid is highlighted in **green**.

### Specific Instances

```
rmv_super KTURN 1,3,5           # Only instances 1, 3, 5
```

### Cross-PDB Superimposition

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif

rmv_fetch 4V88
rmv_load_motif

rmv_loaded                      # Check tags: 1S72_S3, 4V88_S3
rmv_super KTURN                 # All K-TURNs across both PDBs
rmv_super KTURN, 1S72_S3, 4V88_S3   # Filter by PDB+source tags
rmv_super KTURN, 1S72_S3        # Only instances from 1S72
```

### Motif Name Aliasing

The commands resolve common name variants automatically:

| User Input | Resolves To |
|------------|-------------|
| `KTURN`, `K_TURN`, `KINK-TURN`, `KINK_TURN` | K-TURN |
| `CLOOP`, `C_LOOP` | C-LOOP |
| `SARCIN`, `SARCINRICIN`, `SARCIN_RICIN` | SARCIN-RICIN |
| `TLOOP`, `T_LOOP` | T-LOOP |
| `ELOOP`, `E_LOOP` | E-LOOP |
| `REVERSE-KTURN`, `REVERSE_KTURN`, `REVERSE_K_TURN`, `REVERSEKTURN` | REVERSE-K-TURN |

### With Padding

```
rmv_super KTURN, padding=3      # Superimpose with ±3 flanking residues
rmv_align SARCIN-RICIN, padding=2
```

### Output

A summary table is printed showing:
- The medoid instance and its average RMSD
- Each instance's color and RMSD to the medoid
- Overall average RMSD and any skipped pairs

---

## 11. Label Chain IDs (Advanced)

### Background

PDB/CIF files contain two types of chain identifiers:

- **auth_asym_id**: Chain IDs assigned by the authors (e.g., `0`, `9`, `A`, `B`). This is the default.
- **label_asym_id**: System-assigned chain IDs (e.g., `A`, `AA`, `BA`, `CA`). Used by some annotation tools.

### When to Use Label Mode

Use `cif_use_auth=0` when:
- Your user annotations (FR3D/RMS/RMSX) use label chain IDs
- You need label_asym_id chains for your analysis
- Your PDB has unusual auth chain IDs that conflict with annotation tools

### Loading with Label Chain IDs

```
rmv_fetch 1S72, cif_use_auth=0
```

The plugin will:
1. Set PyMOL to use `label_asym_id` as the chain property
2. Parse the CIF file's `_atom_site` loop to build an `auth → label` chain mapping
3. Automatically remap motif data (which uses auth chains) to label chains

### Verify Chain Mode

```
rmv_chains
```

Output:

```
Structure: 1S72  |  cif_use_auth = 0 (label_asym_id)  |  Chains: 293
Label chains: A AA AB AC AD ...
```
![chain_id Pymol command ss](images/chain_id.png)

### Example: Auth vs. Label

Default mode (`cif_use_auth=1`):
```
rmv_fetch 1S72
rmv_summary SARCIN-RICIN
# Shows: Chain 0:158-164, 9:75-81
```

Label mode (`cif_use_auth=0`):
```
rmv_fetch 1S72, cif_use_auth=0
rmv_summary SARCIN-RICIN
# Shows: Chain A:158-164, B:75-81
```

---

## 12. Tips and Tricks

### Force API Refresh (Bypass Cache)

```
rmv_refresh              # Re-fetch from current source, ignoring 30-day cache
```

### Check Current Source

```
rmv_source info          # Show currently active source
```

### Quick Full Pipeline

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif
rmv_show HL
```

### Compare Two Sources Side by Side

```
rmv_fetch 1S72
rmv_db 3               # BGSU
rmv_load_motif
rmv_show HL            # Creates: HL_ALL_S3

rmv_db 1               # Switch to Atlas
rmv_load_motif
rmv_show HL            # Creates: HL_ALL_S1 (BGSU objects persist)
```

### Complete Exploration Workflow

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif

rmv_summary                    # What motifs are available?
rmv_summary GNRA               # How many GNRA instances?
rmv_view GNRA 1                # Quick highlight instance 1
rmv_show GNRA                  # Color them on structure
rmv_show GNRA 1                # Zoom to instance 1
rmv_show GNRA 2                # Zoom to instance 2
rmv_show ALL                   # Show all motif types
rmv_save GNRA sticks           # Save all GNRA images as sticks
rmv_save current               # Save current view
rmv_save GNRA cif              # Export GNRA structures as mmCIF
```

### Reset Everything

```
rmv_reset              # Delete all objects, reset plugin to default settings
```

This clears the entire PyMOL session (all objects, selections) and resets the plugin state (source, PDB, colors, filters) — like a fresh start without restarting PyMOL.

### View Highlighting & Reset

Highlight motif regions without creating objects:

```
rmv_view all                      # Color all motif regions on structure
rmv_view K-TURN                   # Highlight only K-TURN regions
rmv_view K-TURN 1                 # Zoom to instance #1

rmv_view K-TURN hide              # Reset K-TURN coloring to gray
rmv_view all hide                 # Reset all coloring to gray
rmv_view hide                     # Same as rmv_view all hide
```

---

## 13. Quick Reference Card

| Step | Command | Purpose |
|------|---------|---------|
| **Load** | `rmv_fetch <PDB>` | Download and load structure |
| | `rmv_fetch /path/to/file.cif` | Load local PDB/mmCIF file |
| **Sources** | `rmv_sources` | List all data sources |
| | `rmv_db <N>` | Select data source (1–7) |
| | `rmv_db <N1> <N2>` | Combine multiple sources |
| | `rmv_db <N> /path` | Set custom data path (5–7) |
| | `rmv_source info` | Show active source info |
| **Motifs** | `rmv_load_motif` | Fetch motif data |
| | `rmv_refresh` | Force API cache refresh |
| **Summary** | `rmv_summary` | View motif counts |
| | `rmv_summary <TYPE>` | View instances of a type |
| | `rmv_summary <TYPE> <N>` | View single instance details |
| **Render** | `rmv_show <TYPE>` | Render all instances of a type |
| | `rmv_show <TYPE> <N>` | Render a specific instance |
| | `rmv_show <TYPE> <N1>,<N2>` | Render multiple specific instances |
| | `rmv_show <TYPE>, padding=N` | Expand view ±N residues |
| | `rmv_show ALL` | Show all motif types with objects |
| **Filter** | `rmv_show <TYPE> rms` | Show source-specific instances (combine) |
| | `rmv_show <TYPE> shared` | Show shared instances (combine) |
| **View** | `rmv_view all` | Highlight all motif regions |
| | `rmv_view <TYPE>` | Highlight a motif type |
| | `rmv_view <TYPE> <N>` | Highlight and zoom to instance |
| | `rmv_view hide` | Reset all view coloring |
| | `rmv_view <TYPE> hide` | Reset coloring for one type |
| | `rmv_view all hide` | Reset all coloring |
| **Super** | `rmv_super <TYPE>` | Superimpose instances (medoid, sequence-independent) |
| | `rmv_super <TYPE> <N1>,<N2>` | Superimpose specific instances |
| | `rmv_super <TYPE>, padding=N` | With expanded context |
| | `rmv_super <TYPE>, <PDB_SRC1>, <PDB_SRC2>` | Cross-PDB superimposition |
| **Align** | `rmv_align <TYPE>` | Sequence-based alignment (medoid) |
| **Save** | `rmv_save ALL` | Save all motif images |
| | `rmv_save <TYPE> <N> sticks` | Save instance with representation |
| | `rmv_save current` | Save current view (high-res) |
| | `rmv_save current my_fig` | Save with custom filename |
| | `rmv_save ALL cif` | Export all motif structures as mmCIF |
| | `rmv_save <TYPE> cif` | Export motif type as mmCIF |
| | `rmv_save <TYPE> <N> cif` | Export specific instance as mmCIF |
| **Info** | `rmv_chains` | Chain ID diagnostics |
| | `rmv_loaded` | List loaded PDBs and sources |
| | `rmv_help` | Full command reference |
| **Color** | `rmv_color <TYPE> <COLOR>` | Change motif color |
| | `rmv_colors` | View color legend |
| | `rmv_bg_color <COLOR>` | Change backbone color |
| **P-value** | `rmv_db <N> off` | Disable filtering (6–7) |
| | `rmv_db <N> on` | Enable filtering (6–7) |
| | `rmv_db <N> <MOTIF> <PVAL>` | Custom P-value threshold |
| **Toggle** | `rmv_toggle <TYPE> on/off` | Toggle visibility |
| **Reset** | `rmv_reset` | Delete all objects & reset plugin |

---

*RSMViewer v1.0.0 — CBB LAB KU @Dr. Zhong*
