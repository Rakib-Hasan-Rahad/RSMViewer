# Multi-Source Combine Pipeline

RSMViewer can combine motif data from multiple sources and merge them with smart
deduplication using the **cascade merge** pipeline.

---

## Quick Start

```
rmv_fetch 1S72
rmv_db 1 3                  # Combine Atlas [1] + BGSU [3]
rmv_load_motif              # Fetches from both → deduplicates → merged
rmv_summary                 # Unified results from both sources
```

### Three-Source Combine

```
rmv_db 1 3 4                # Atlas + BGSU + Rfam API
rmv_load_motif
```

### Custom Jaccard Threshold

```
rmv_db 1 3 jaccard_threshold=0.80   # Stricter deduplication
```

---

## Pipeline Overview

The combine pipeline has 5 stages:

```
Fetch → Enrich → Source Stamp → Within-Source Dedup → Cascade Merge
```

### Stage 1 — Fetch

Each source independently fetches motifs for the PDB:

```
Atlas → [HL(35), IL(19), J3(18), ...]
BGSU  → [HL(42), IL(23), J3(21), GNRA(15), SARCIN-RICIN(8), ...]
```

### Stage 2 — Enrich (Homolog Enrichment)

Generic motif names are enriched using NR homolog representative lookup via
`homolog_enricher.py`. For example:

```
Atlas HL_345 → matches BGSU GNRA_123 (same residues) → rename to GNRA
```

Uses `nrlist_4.24_all.csv` for the BGSU NR representative set.

### Stage 3 — Source Stamping

Each motif instance is tagged with `_source_id` and `_source_label` metadata.
This enables the SOURCE column in instance tables and source attribution reporting.

### Stage 4 — Within-Source Deduplication

Exact duplicates within the same source are removed. Two instances are duplicates
if they have the **same set of `(chain, residue_number)` pairs**.

### Stage 5 — Cascade Merge

Implemented in `database/cascade_merger.py`.

---

## Cascade Merge Algorithm

### Core Class

```python
class CascadeMerger:
    def __init__(self, jaccard_threshold: float = 0.60)

    def merge_sources(self, ordered_sources, source_labels=None):
        """Merge multiple source results with Jaccard deduplication."""
```

### Right-to-Left Merge

Sources are merged in **right-to-left order** (later sources have higher priority):

```
rmv_db 1 3 4    →    result = merge(merge(source_4, source_3), source_1)
```

The right-most source (source 4 in this example) has the highest priority.

### Pairwise Merge Logic

For each pair of sources `(reference, updater)` where reference has higher priority:

1. Reference keeps **all** its motif instances
2. For each updater motif instance:
   - Extract **base category** via `_get_base_category()` (strips trailing `_1`, `_2`, maps variants using `MOTIF_CANONICAL_MAP`)
   - Check Jaccard overlap **only within the same base category**
   - If Jaccard ≥ threshold (default 0.60) → **discard** updater motif, annotate reference with `_also_found_in`
   - Otherwise → **keep** updater motif (it's unique)

### Jaccard Similarity

```python
def _jaccard(set_a: Set, set_b: Set) -> float:
    """J(A, B) = |A ∩ B| / |A ∪ B|"""
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0
```

Residue sets use `(chain, residue_number)` tuples for chain-aware comparison.

### Example

Given two instances:
- **BGSU GNRA #3:** residues `{(A, 160), (A, 161), (A, 162), (A, 163)}`
- **Atlas HL #12:** residues `{(A, 160), (A, 161), (A, 162), (A, 163), (A, 164)}`

Jaccard = |{160,161,162,163}| / |{160,161,162,163,164}| = 4/5 = 0.80

Since 0.80 ≥ 0.60 (threshold), these are considered the **same motif**. The BGSU version
(higher priority) is kept, and the Atlas version is discarded with an `_also_found_in` annotation.

---

## Motif Canonical Map

The `MOTIF_CANONICAL_MAP` in cascade_merger.py maps variant names to their base categories
for deduplication purposes:

```
K-TURN-1, K-TURN-2, pK-TURN → K-TURN
SARCIN-RICIN-1, SARCIN-RICIN-2 → SARCIN-RICIN
RIGHT-ANGLE-2, RIGHT-ANGLE-3 → RIGHT-ANGLE
```

This ensures that `K-TURN-1` from Rfam and `K-TURN` from BGSU are compared against
each other during deduplication.

---

## Source-Tagged PyMOL Objects

When combining or switching sources, PyMOL objects include a source suffix:

```
rmv_db 3
rmv_show HL            # Creates: HL_ALL_S3, HL_1_S3, ...

rmv_db 1 3             # Combine mode
rmv_show HL            # Creates: HL_ALL_S1_3, HL_1_S1_3, ...
```

---

## Adjusting the Threshold

| Threshold | Effect |
|-----------|--------|
| 0.40 | Aggressive dedup — more motifs treated as duplicates |
| **0.60** | **Default** — balanced deduplication |
| 0.80 | Conservative — only very similar motifs deduplicated |
| 1.00 | Only exact duplicates removed |

```
rmv_db 1 3 jaccard_threshold=0.80
rmv_load_motif
```

---

## Related Commands

| Command | Purpose |
|---------|---------|
| `rmv_db <N1> <N2> [N3...]` | Select sources to combine |
| `rmv_load_motif` | Fetch and merge from all selected sources |
| `rmv_summary` | View merged results |
| `rmv_show <TYPE>` | Visualize merged motifs |
| `rmv_source info` | Show currently active source(s) |
| `rmv_loaded` | List loaded PDBs and sources |

---

## Further Reading

See [COMBINE_GUIDE.md](../COMBINE_GUIDE.md) for a detailed walkthrough of the combine
pipeline with code references and worked examples.
