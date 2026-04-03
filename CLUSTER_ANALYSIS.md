# Cluster Analysis — RSMViewer

## Overview

The **Cluster Analysis** feature (Source 8) lets you load groups of structurally
related RNA motif instances from a CSV file, visualize them side-by-side in
PyMOL with distinct colours, and structurally align any pair of objects to
obtain RMSD values.

---

## CSV File Format

Each line defines one **cluster**.  Fields are comma-separated.

```
cluster_name,PDBID_CHAIN:range1_range2_...,PDBID_CHAIN:range1_range2_...,...
```

| Field | Description |
|-------|-------------|
| `cluster_name` | Unique label for the cluster (e.g. `ML1_2`) |
| `PDBID_CHAIN` | 4-char PDB ID + underscore + chain ID (e.g. `1S72_0`, `4QCN_A`) |
| `range` | Residue range in `start-end` format (e.g. `2552-2555`) |

- Multiple residue ranges per entry are separated by `_` (underscore).
- A trailing comma at the end of a line is tolerated.
- Lines starting with `#` are treated as comments.
- The same PDB + chain can appear more than once in a cluster (duplicates get
  an auto-incrementing index suffix).

### Example CSV

```csv
ML1_2,1S72_0:2552-2555_2580-2582_2596-2602,3J7Y_A:3004-3007_3032-3034_3048-3054,4QCN_A:2529-2532_2557-2559_2573-2579
ML8_1,3FO4_A:45-54_72-75_21-25,3LA5_A:45-54_72-75_21-25,4TZX_X:45-54_72-75_21-25
```

### Bundled CSV

A default CSV (`ML_train_Motif_input.csv`) is bundled at:

```
rsmviewer/database/cluster_analysis/ML_train_Motif_input.csv
```

It contains 5 clusters (39 total entries): `ML1_2`, `ML2_1`, `ML3_1`,
`ML8_1`, `H`.

---

## Commands

| Command | Description |
|---------|-------------|
| `rmv_db 8` | Load the bundled cluster CSV |
| `rmv_db 8 /path/to/file.csv` | Load a custom cluster CSV |
| `rmv_cluster` | List all loaded clusters |
| `rmv_cluster <NAME>` | Visualize a cluster — fetches PDBs, creates per-entry objects with unique colours, prints summary table, zooms |
| `rmv_cluster_clear` | Remove all `cluster_*` objects from the PyMOL session |
| `align <mobile>, <target>` | PyMOL's native sequence-dependent alignment (reports RMSD) |
| `super <mobile>, <target>` | PyMOL's native sequence-independent superposition (reports RMSD) |

---

## Workflow

### 1. Load cluster data

```
rmv_db 8
```

Output lists the available clusters and their entry counts.

To use your own CSV instead:

```
rmv_db 8 /path/to/my_clusters.csv
```

### 2. List clusters

```
rmv_cluster
```

### 3. Visualize a cluster

```
rmv_cluster ML1_2
```

This will:
1. Fetch each PDB structure (skips already-loaded ones).
2. Remove solvent from fetched structures.
3. Create a PyMOL object for each cluster entry containing only the motif
   residues (named `cluster_ML1_2_1S72_0`, `cluster_ML1_2_3J7Y_A`, etc.).
4. Assign a distinct colour to each object (cycles through 15 colours).
5. Print a summary table with object names, PDB IDs, chains, and residue
   ranges.
6. Zoom to the first object.

### 4. Align two objects

Use PyMOL's built-in alignment commands directly:

Sequence-dependent alignment:

```
align cluster_ML1_2_1S72_0, cluster_ML1_2_3J7Y_A
```

Sequence-independent superposition:

```
super cluster_ML1_2_1S72_0, cluster_ML1_2_3J7Y_A
```

Both commands print RMSD, atom count, and alignment details in the PyMOL
console.

### 5. Clean up

```
rmv_cluster_clear
```

Removes all objects whose names start with `cluster_`.

---

## Object Naming Convention

```
cluster_<CLUSTER_NAME>_<PDB_ID>_<CHAIN>[_<INDEX>]
```

- The `_<INDEX>` suffix is only added when the same PDB + chain appears
  multiple times in a cluster (e.g. `cluster_ML2_1_4QCN_A` and
  `cluster_ML2_1_4QCN_A_1`).

---

## Colour Palette

Objects are coloured in order using this 15-colour palette (cycles if more
than 15 entries):

red, blue, green, orange, magenta, cyan, yellow, salmon, lime, slate,
hotpink, teal, olive, violet, deepteal

---

## Tips

- You can use PyMOL's `align` / `super` commands with **any** PyMOL objects — not just
  cluster objects.  For example, align motif objects from different sources:
  ```
  super SARCIN_RICIN_3_1S72_S3, SARCIN_RICIN_3_1S72_S7
  ```
- After alignment, the mobile object is moved onto the target.  Use PyMOL's
  `undo` to revert.
- To visualize multiple clusters at once, simply call `rmv_cluster` for each
  cluster name.  Use `rmv_cluster_clear` to remove them all.
