# Structural Superimposition (Medoid)

RSMViewer implements automatic medoid-based structural superimposition for comparing multiple
instances of the same motif type. The medoid (most representative instance) is identified
and all other instances are aligned onto it.

---

## Commands

| Command | Method | Best For |
|---------|--------|----------|
| `rmv_super <TYPE>` | `cmd.super()` | Sequence-independent comparison |
| `rmv_align <TYPE>` | `cmd.align()` | Sequence-similar comparison |

---

## Usage

### Basic Superimposition

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif
rmv_super KTURN              # Superimpose all K-TURN instances
```

### Specific Instances

```
rmv_super KTURN 1,3,5        # Only instances 1, 3, 5
```

### Cross-PDB Superimposition

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif

rmv_fetch 4V88
rmv_load_motif

rmv_super KTURN                  # All K-TURNs across both PDBs
rmv_super KTURN pdb=1S72,4V88    # Explicit PDB selection
rmv_super KTURN pdb=1S72         # Only instances from 1S72
```

### With Padding

```
rmv_super K-TURN, padding=3      # Extended residue context
rmv_align SARCIN-RICIN, padding=2
```

---

## How It Works

The superimposition pipeline is implemented in `alignment.py`.

### Step 1 — Pairwise RMSD Matrix

```python
def compute_pairwise_rmsd(objects, method="super"):
    """Compute NxN RMSD matrix using temporary copies."""
```

For each pair of instances `(i, j)`:
1. Creates a temporary PyMOL copy of instance `i`
2. Runs `cmd.super` (or `cmd.align`) against instance `j`
3. Records the RMSD value
4. Failed pairs get `float('inf')` (skipped in medoid calculation)

### Step 2 — Medoid Selection

```python
def find_medoid(matrix):
    """Return (medoid_index, avg_rmsd_list).
    Medoid = row with minimum average RMSD (excluding inf pairs)."""
```

For each instance:
- Computes the average RMSD to all other instances (excluding infinite/failed pairs)
- The instance with the **smallest average RMSD** is selected as the medoid

The medoid is the most "central" or representative instance in the structural ensemble.

### Step 3 — Align All to Medoid

```python
def superimpose_onto_medoid(objects, medoid_idx, method="super"):
    """Move every non-medoid object onto the medoid in place."""
```

All non-medoid instances are structurally aligned to the medoid using the selected method.

### Step 4 — Coloring

- **Medoid:** Always colored **green**
- **Other instances:** Cycle through `SUPER_COLORS` (red, blue, orange, magenta, cyan, ...)
- Each instance gets a unique color for easy identification

---

## Output

A summary table is printed:

```
==========================================================
  SUPERIMPOSITION SUMMARY — K-TURN (5 instances)
==========================================================
  INSTANCE    PDB     RMSD       COLOR      NOTE
----------------------------------------------------------
  3           1S72    0.000 Å    green      MEDOID
  1           1S72    1.234 Å    red
  2           1S72    1.567 Å    blue
  4           4V88    2.012 Å    orange
  5           4V88    2.345 Å    magenta
----------------------------------------------------------
  Average RMSD: 1.790 Å
  Skipped pairs: 0
==========================================================
```

---

## Motif Name Aliasing

The superimposition commands resolve common name variants automatically:

| Input | Matches |
|-------|---------|
| `KTURN` | K-TURN |
| `KINK-TURN` | K-TURN |
| `SARCIN` | SARCIN-RICIN |

---

## super vs align

| Feature | `rmv_super` | `rmv_align` |
|---------|------------|------------|
| PyMOL function | `cmd.super()` | `cmd.align()` |
| Sequence dependence | No | Yes |
| Best for | Different organisms, divergent sequences | Same species, similar sequences |
| Structural comparison | Pure 3D overlay | Sequence-guided overlay |

---

## Related Commands

| Command | Purpose |
|---------|---------|
| `rmv_show <TYPE>` | Visualize before superimposing |
| `rmv_summary <TYPE>` | Check instance count |
| `rmv_pair <TYPE> <N>` | Analyze base pairs of an instance |
| `rmv_save <TYPE>` | Save images after superimposition |
