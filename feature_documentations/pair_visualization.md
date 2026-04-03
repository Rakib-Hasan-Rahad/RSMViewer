# Base-Pair Visualization

RSMViewer provides base-pair analysis with **Leontis-Westhof (LW)** classification labels
for motif instances, implemented in `pair_visualizer.py`.

---

## Commands

| Command | Purpose |
|---------|---------|
| `rmv_pair <TYPE> <N>` | Show base-pair interactions for instance N |
| `rmv_pair_batch <TYPE>` | Batch base-pair analysis for all instances |

---

## Usage

### Single Instance

```
rmv_pair SARCIN-RICIN 1       # Show base pairs for instance 1
rmv_pair KTURN 3              # Show base pairs for K-TURN instance 3
rmv_pair GNRA 1               # Show base pairs for GNRA instance 1
```

### Batch Analysis

```
rmv_pair_batch SARCIN-RICIN   # Analyze all sarcin-ricin instances
rmv_pair_batch GNRA           # Analyze all GNRA instances
```

---

## How It Works

The pair visualizer (`pair_visualizer.py`) is function-based (no class).

### Step 1 — Parse Pair Descriptor

```python
def parse_pair_descriptor(descriptor: str) -> dict
    # Parses: pdbid_chain1_resnum1_chain2_resnum2
```

Extracts the two paired residues' chain and residue number from the pair descriptor.

### Step 2 — Detect Nucleotide Identity

```python
def _get_resname_from_structure(obj_name, chain, resi) -> str
    # Auto-detects nucleotide identity (A, G, C, U) from loaded structure
```

### Step 3 — Place Edge Labels

For each nucleotide, the visualizer identifies the Watson-Crick (W), Hoogsteen (H),
and Sugar (S) edges using atom coordinates:

```python
EDGE_ATOMS = {
    "A": {"W": ["N1","C2","N6"],        "H": ["C5","C6","N7","C8"],  "S": ["N3","C2"]},
    "G": {"W": ["O6","N1","N2"],        "H": ["C5","N7","C8"],       "S": ["N3","C2","N2"]},
    "C": {"W": ["N3","O2","N4"],        "H": ["C5","C6"],            "S": ["N1","C2","O2"]},
    "U": {"W": ["O4","N3","O2"],        "H": ["C5","C6"],            "S": ["N1","C2","O2"]},
}
```

For each edge, a **pseudoatom** is placed at the center of mass of the defining atoms,
labeled "W", "H", or "S".

```python
def _place_edge_label(pair_obj, chain, resi, resname, label_prefix)
    # Places W, H, S pseudoatom labels at centroid of defining atoms
```

### Step 4 — Draw Interactions

Dashed lines are drawn between paired bases in the PyMOL viewport. Each pair is labeled
with its LW classification (e.g., `cWW` for cis Watson-Crick/Watson-Crick, `tHS` for
trans Hoogsteen/Sugar).

---

## Leontis-Westhof Classification

The LW system classifies base pairs by two properties:

### Edge Types

| Edge | Symbol | Atoms Involved |
|------|--------|---------------|
| **Watson-Crick** | W | Standard hydrogen-bonding face |
| **Hoogsteen** | H | Major groove face |
| **Sugar** | S | Minor groove (2'-OH) face |

### Orientation

| Prefix | Meaning |
|--------|---------|
| **c** | cis — glycosidic bonds on same side |
| **t** | trans — glycosidic bonds on opposite sides |

### Common LW Labels

| Label | Full Name |
|-------|-----------|
| `cWW` | cis Watson-Crick / Watson-Crick (canonical) |
| `tWW` | trans Watson-Crick / Watson-Crick |
| `cWH` | cis Watson-Crick / Hoogsteen |
| `tWH` | trans Watson-Crick / Hoogsteen |
| `cWS` | cis Watson-Crick / Sugar |
| `tWS` | trans Watson-Crick / Sugar |
| `cHH` | cis Hoogsteen / Hoogsteen |
| `tHH` | trans Hoogsteen / Hoogsteen |
| `cHS` | cis Hoogsteen / Sugar |
| `tHS` | trans Hoogsteen / Sugar |
| `cSS` | cis Sugar / Sugar |
| `tSS` | trans Sugar / Sugar |

---

## Typical Workflow

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif

rmv_summary SARCIN-RICIN         # Check available instances
rmv_show SARCIN-RICIN 1          # Visualize structure
rmv_pair SARCIN-RICIN 1          # Overlay base-pair labels
rmv_pair_batch SARCIN-RICIN      # Batch analysis for all instances
```

---

## Related Commands

| Command | Purpose |
|---------|---------|
| `rmv_show <TYPE> <N>` | Visualize instance before pair analysis |
| `rmv_view <TYPE> <N>` | Quick highlight before pair analysis |
| `rmv_summary <TYPE> <N>` | Check residue details |
| `rmv_super <TYPE>` | Superimpose instances for structural comparison |
| `rmv_save <TYPE> <N>` | Save image after pair visualization |
