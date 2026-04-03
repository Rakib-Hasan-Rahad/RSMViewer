# `rmv_super` / `rmv_align` — Medoid-Based Superimposition Plan

## Goal

Given multiple instances of the same RNA motif (within one PDB or across
multiple PDBs), automatically identify the **medoid** — the instance whose
average pairwise RMSD to all other instances is lowest — and superimpose every
other instance onto it.  This produces a "consensus view" centred on the most
representative structure. Each superimposed instance gets a unique colour.

---

## Key Concepts

| Term | Definition |
|------|-----------|
| **Motif instance** | A single occurrence of a motif type (e.g., one K-TURN in 1S72). Already represented as `MotifInstance` in `base_provider.py`. |
| **Pairwise RMSD** | The structural distance between two instances, computed by PyMOL's `cmd.super()` or `cmd.align()`. |
| **Medoid** | The instance *i* that minimises $\frac{1}{N-1}\sum_{j \neq i} \text{RMSD}(i,j)$.  Unlike a centroid (which is a virtual average coordinate), the medoid is an actual member of the set. |
| **Distance matrix** | An $N \times N$ symmetric matrix of pairwise RMSDs. |

---

## Commands

Two separate commands, identical pipeline, different PyMOL backend:

| Command | PyMOL function | Best for |
|---------|---------------|----------|
| `rmv_super` | `cmd.super()` — sequence-independent | Comparing same motif across different organisms/PDBs |
| `rmv_align` | `cmd.align()` — sequence-dependent | Comparing motifs with known sequence similarity |

### Syntax

```
rmv_super <motif_type> [instance_list] [pdb=PDB1,PDB2,...]
rmv_align <motif_type> [instance_list] [pdb=PDB1,PDB2,...]
```

### Examples

```
rmv_super KTURN                         All K-TURN instances (sequence-independent)
rmv_super KTURN 1,3,5                   Only instances 1, 3, 5
rmv_super KTURN pdb=1S72                Only instances from 1S72
rmv_super KTURN pdb=1S72,4V88           Instances from these 2 PDBs
rmv_super KTURN pdb=1S72,4V88,4V9F      3+ PDBs supported
rmv_super IL                            All internal loops
```

For manual pairwise alignment of any two PyMOL objects, use PyMOL natively:
```
super obj1, obj2
align obj1, obj2
```

---

## Workflow (User Perspective)

### Single-PDB workflow

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif
rmv_super KTURN                     # Medoid superimposition of all K-TURNs
rmv_super KTURN 1,3,5               # Or specific instances
```

### Cross-PDB workflow

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif                      # Creates KTURN_ALL_1S72_S3, etc.

rmv_fetch 4V88
rmv_load_motif                      # Creates KTURN_ALL_4V88_S3, etc.

rmv_super KTURN                     # ALL K-TURN instances across both PDBs
rmv_super KTURN pdb=1S72,4V88       # Explicit: only these 2 PDBs
rmv_super KTURN pdb=1S72            # Only instances from 1S72
rmv_super KTURN pdb=1S72,4V88,4V9F  # 3+ PDBs
```

The `pdb=` argument is **optional**. Without it, all loaded instances of
that motif type are included regardless of PDB. When provided, it accepts
a comma-separated list of any number of PDB IDs.

---

## Motif Name Aliasing Problem

### The Problem

Different data sources use different names for the same motif:

| Source | K-turn example | Sarcin-ricin example |
|--------|---------------|---------------------|
| Atlas (1) / BGSU API (3) | Not in their vocabulary (they use HL, IL, J3–J7) | Not in their vocabulary |
| Rfam (2) | `k-turn-1`, `k-turn-2` | `sarcin-ricin-1`, `sarcin-ricin-2` |
| FR3D (5) | Listed in CSV as-is | Listed in CSV as-is |
| RMS (6) | `Kturn`, `kturn` | `sarcin_ricin` |
| RMSX (7) | `K-TURN` (normalised by `_KNOWN_MOTIFS`) | `SARCIN-RICIN` |

After loading, motif types are stored as **uppercase keys** in
`loaded_motifs` dict (e.g., `KTURN`, `K-TURN`, `KINK-TURN`).
These are different keys even though they represent the same structural motif.

### Current State

The RMSX converter already has a `_KNOWN_MOTIFS` mapping that normalises
some variations:

```python
_KNOWN_MOTIFS = {
    'KTURN': 'K-TURN', 'K-TURN': 'K-TURN', 'KINK-TURN': 'K-TURN',
    'CLOOP': 'C-LOOP', 'C-LOOP': 'C-LOOP', 'C_LOOP': 'C-LOOP',
    'SARCIN': 'SARCIN-RICIN', 'SARCIN-RICIN': 'SARCIN-RICIN',
    ...
}
```

But this only applies during RMSX file parsing. It does **not** help at
`rmv_super` time when the user types a name that doesn't exactly match
the key in `loaded_motifs`.

### Proposed Solution: Alias Lookup at `rmv_super`/`rmv_align` Time

When the user types `rmv_super KTURN`, and `KTURN` is not a key in
`loaded_motifs`, do a fuzzy/alias match:

1. **Exact match** (case-insensitive) → use it.
2. **Canonical alias table** — a global mapping of known motif name
   variations to a canonical form:
   ```python
   MOTIF_ALIASES = {
       'KTURN': 'K-TURN', 'KINK-TURN': 'K-TURN', 'KINK_TURN': 'K-TURN',
       'K_TURN': 'K-TURN',
       'CLOOP': 'C-LOOP', 'C_LOOP': 'C-LOOP',
       'SARCIN': 'SARCIN-RICIN', 'SARCINRICIN': 'SARCIN-RICIN',
       'SARCIN_RICIN': 'SARCIN-RICIN',
       'REVERSE-KTURN': 'REVERSE-K-TURN', 'REVERSE_KTURN': 'REVERSE-K-TURN',
       'ELOOP': 'E-LOOP', 'E_LOOP': 'E-LOOP',
       'TLOOP': 'T-LOOP', 'T_LOOP': 'T-LOOP',
       ...
   }
   ```
3. **Substring match fallback** — if neither exact nor alias works, check
   if the user's input is a substring of any loaded key (or vice versa).
   E.g., `KTURN` matches `K-TURN`, `SARCIN` matches `SARCIN-RICIN`.
4. **Report ambiguity** — if multiple keys match, list them and ask the
   user to be more specific.

This alias resolution should be applied in `rmv_super`, `rmv_align`, and
ideally also in `rmv_show` and `rmv_summary` for consistency.

### Cross-Source Matching for Cross-PDB

When doing `rmv_super KTURN pdb=1S72,4V88` where 1S72 was loaded from
RMSX (motif named `K-TURN`) and 4V88 from BGSU API (where K-turns may
appear as part of `IL` or custom names), the alias table ensures:

- User types `KTURN` → resolves to both `K-TURN` and `KTURN` keys
- All matching instances from all PDBs are collected into one pool
- The medoid pipeline runs on the combined set

---

## Architecture

### Reuse `rsmviewer/alignment.py` (repurposed)

```
alignment.py
├── _validate_objects(*names)              # (existing) check object existence
├── _print_alignment_result(...)           # (existing) format RMSD table
│
├── MOTIF_ALIASES                          # NEW — canonical name mapping
├── resolve_motif_name(user_input, loaded) # NEW — alias + fuzzy lookup
├── compute_pairwise_rmsd(objects, method) # NEW — pairwise RMSD matrix
├── find_medoid(rmsd_matrix)               # NEW — index of medoid
├── superimpose_onto_medoid(objects, idx, method) # NEW
├── color_superimposed(objects, medoid_idx)# NEW — unique colours
├── print_medoid_report(...)               # NEW — formatted summary
│
└── register_alignment_commands()          # NEW — registers rmv_super + rmv_align
```

### Integration Points

| Existing component | How it's used |
|--------------------|---------------|
| `loader.py` | `_create_single_instance_object()` creates per-instance objects. Medoid command reuses this for batch creation. |
| `gui.py` | Register `rmv_super` and `rmv_align` in `_RMV_COMMANDS` and `_CMD_SUFFIXES`, call `register_alignment_commands()`. |
| `plugin.py` | Import and call registration. |
| `colors.py` | Use existing palette for unique colouring. |

---

## Performance Strategy: Object Creation

### The Problem

`rmv_show KTURN` creates **one combined object** (`KTURN_ALL_1S72_S3`) for all
K-TURN instances fused together. Individual instance objects (e.g.,
`KTURN_1_1S72_S3`, `KTURN_2_1S72_S3`) are only created on-demand via
`rmv_show KTURN 1`. Creating many objects at once in PyMOL is slow.

### The Solution: Two-Phase Pipeline

**Phase 1 — Batch object creation (optimised):**

1. Gather all instance selections as strings first (pure Python, instant).
2. Use `cmd.feedback("disable", "all", "actions")` to suppress PyMOL's
   per-object GUI updates during creation.
3. Create all N objects in a tight loop.
4. Re-enable feedback, then do a single `cmd.rebuild()` + `cmd.refresh()`.
5. Print a progress indicator: `Creating objects... 5/23`

```python
def _batch_create_instance_objects(motif_type, instances, structure_name, source_suffix):
    """Create all instance objects with suppressed GUI updates."""
    cmd.feedback("disable", "all", "actions")
    try:
        created = []
        for i, detail in enumerate(instances, 1):
            obj_name = _create_one(motif_type, i, detail, structure_name, source_suffix)
            if obj_name:
                created.append(obj_name)
            if i % 10 == 0:
                print(f"  Creating objects... {i}/{len(instances)}")
        return created
    finally:
        cmd.feedback("enable", "all", "actions")
        cmd.rebuild()
        cmd.refresh()
```

**Phase 2 — RMSD computation + superimposition** (uses temporary copies).

User sees "Creating objects..." first, then "Computing RMSD matrix..." — clear
progress indication.

**Timing estimate:** Object creation ~0.1–0.5s per object. For 20 instances,
2–10 seconds. RMSD matrix for 20 instances = $\binom{20}{2} = 190$
superimpositions at ~0.01s each ≈ 2 seconds.

---

## Detailed Algorithm

### Step 1 — Resolve Motif Name + Collect Target Objects

```python
def _collect_motif_objects(motif_type, indices=None, pdb_filter=None):
    """
    Get PyMOL object names for a motif type's rendered instances.
    
    Object naming: MOTIF_NO_PDB_SOURCE (e.g., KTURN_1_1S72_S3)
    
    Steps:
    1. Resolve user's motif name via alias table + fuzzy match
    2. Find all matching keys in loaded_motifs
    3. Find/create individual instance objects
    4. Apply pdb= and index filters
    """
    all_objects = cmd.get_object_list()
    
    # Resolve aliases — may match multiple keys (e.g., "KTURN" → ["K-TURN"])
    resolved_keys = resolve_motif_name(motif_type, loaded_motifs)
    
    matching = []
    for key in resolved_keys:
        prefix = key + "_"
        key_matches = [o for o in all_objects 
                       if o.upper().startswith(prefix) and "_ALL_" not in o.upper()]
        matching.extend(key_matches)
    
    if pdb_filter:
        pdb_set = {p.upper() for p in pdb_filter}
        matching = [o for o in matching if _extract_pdb(o) in pdb_set]
    
    if indices:
        matching = [o for o in matching if _extract_index(o) in indices]
    
    return sorted(matching)
```

If no individual objects exist yet (user only ran `rmv_show KTURN` which
creates the combined object), the command **auto-creates** them using the
batch method above.

### Step 2 — Pairwise RMSD Matrix

`cmd.super()` physically moves the mobile object. Use temporary copies:

```python
def compute_pairwise_rmsd(objects, method="super"):
    """Compute NxN RMSD matrix using temporary copies."""
    n = len(objects)
    matrix = [[0.0] * n for _ in range(n)]
    skipped = []
    
    align_fn = cmd.super if method == "super" else cmd.align
    
    # Create temporary reference copies (prefixed with _ to hide from panel)
    temps = []
    for i, obj in enumerate(objects):
        tmp = f"_medoid_ref_{i}"
        cmd.create(tmp, obj)
        temps.append(tmp)
    
    total = n * (n - 1) // 2
    done = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            work = "_medoid_work"
            cmd.create(work, temps[i])
            try:
                result = align_fn(work, temps[j])
                rmsd = result[0]
                matrix[i][j] = rmsd
                matrix[j][i] = rmsd
            except Exception:
                matrix[i][j] = float('inf')
                matrix[j][i] = float('inf')
                skipped.append((objects[i], objects[j]))
            cmd.delete(work)
            done += 1
            if done % 50 == 0:
                print(f"  Computing RMSD... {done}/{total}")
    
    for tmp in temps:
        cmd.delete(tmp)
    
    return matrix, skipped
```

### Step 3 — Find Medoid

```python
def find_medoid(matrix):
    """Row index with minimum average RMSD (excluding ∞ pairs)."""
    n = len(matrix)
    avg_rmsd = []
    for i in range(n):
        finite = [matrix[i][j] for j in range(n)
                  if j != i and matrix[i][j] != float('inf')]
        avg_rmsd.append(sum(finite) / len(finite) if finite else float('inf'))
    return avg_rmsd.index(min(avg_rmsd)), avg_rmsd
```

### Step 4 — Superimpose All onto Medoid

```python
def superimpose_onto_medoid(objects, medoid_idx, method="super"):
    """Superimpose every non-medoid object onto the medoid (in place)."""
    medoid_obj = objects[medoid_idx]
    align_fn = cmd.super if method == "super" else cmd.align
    results = []
    for i, obj in enumerate(objects):
        if i == medoid_idx:
            results.append((obj, 0.0, True))
            continue
        try:
            result = align_fn(obj, medoid_obj)
            results.append((obj, result[0], True))
        except Exception:
            results.append((obj, float('inf'), False))
    return results
```

### Step 5 — Unique Colouring

Each instance gets a unique colour. Medoid is highlighted in green.

```python
SUPER_COLORS = [
    "red", "blue", "orange", "magenta", "cyan",
    "yellow", "salmon", "lime", "slate", "hotpink",
    "teal", "olive", "violet", "deepteal", "white",
]
MEDOID_COLOR = "green"

def color_superimposed(objects, medoid_idx):
    for i, obj in enumerate(objects):
        if i == medoid_idx:
            cmd.color(MEDOID_COLOR, obj)
        else:
            cmd.color(SUPER_COLORS[i % len(SUPER_COLORS)], obj)
```

### Step 6 — Report

```
  ┌────────────────────────────────────────────────────────────────┐
  │  rmv_super — KTURN  (5 instances)                             │
  ├────────────────────────────────────────────────────────────────┤
  │  Medoid:  KTURN_3_1S72_S3  (avg RMSD: 1.234 Å)               │
  ├────────────────────────────────────────────────────────────────┤
  │  #   Instance                  Color     RMSD to medoid       │
  │  1   KTURN_1_1S72_S3           red       1.052 Å              │
  │  2   KTURN_2_1S72_S3           blue      1.417 Å              │
  │  3   KTURN_3_1S72_S3           GREEN     — (medoid)           │
  │  4   KTURN_4_4V88_S3           orange    1.198 Å              │
  │  5   KTURN_5_4V88_S3           magenta   1.271 Å              │
  ├────────────────────────────────────────────────────────────────┤
  │  Overall avg RMSD: 1.235 Å                                   │
  │  Skipped pairs: 0                                             │
  └────────────────────────────────────────────────────────────────┘
```

If any pairs were skipped (superimposition failed):
```
  ⚠  1 pair(s) skipped (superimposition failed):
     KTURN_1_1S72_S3 ↔ KTURN_4_4V88_S3
```

---

## Optional Flags

| Flag | Default | Description |
|------|---------|-------------|
| `pdb=PDB1,PDB2,...` | all loaded PDBs | Filter instances by PDB ID. Accepts any number of comma-separated PDB IDs. |
| `color_medoid=<color>` | `green` | Colour for the medoid object |
| `save_matrix=/path.csv` | — | Export the pairwise RMSD matrix as CSV |

---

## Cross-PDB Support

### Default: All Loaded PDBs

```
rmv_super KTURN
```

Collects every KTURN instance across all loaded PDBs. The algorithm is
PDB-agnostic — `cmd.super()` operates on atomic coordinates regardless
of source PDB.

### Filter by PDB (optional)

```
rmv_super KTURN pdb=1S72                # Only 1S72
rmv_super KTURN pdb=1S72,4V88           # Only these 2
rmv_super KTURN pdb=1S72,4V88,4V9F      # 3+ PDBs
```

The `pdb=` argument accepts a comma-separated list of any number of PDB IDs.

### Comparing Motifs Across PDBs

To compare all K-TURNs from PDB A vs PDB B:

```
rmv_fetch 1S72
rmv_db 3
rmv_load_motif

rmv_fetch 4V88
rmv_load_motif

rmv_super KTURN pdb=1S72,4V88
```

The medoid could come from either PDB — it's the instance most structurally
central overall. The report shows which PDB each instance belongs to (visible
in the object name: `KTURN_1_1S72_S3` vs `KTURN_3_4V88_S3`).

The RMSD matrix CSV export (`save_matrix=`) lets you extract cross-PDB vs
within-PDB distances for further analysis.

---

## Edge Cases

| Case | Handling |
|------|----------|
| Only 1 instance | Print: "Need ≥ 2 instances for medoid analysis." |
| 2 instances | Superimpose directly, report RMSD. Both have same avg RMSD. |
| Instances not rendered yet | Auto-create individual objects using batch method (see Performance Strategy). |
| `cmd.super()` fails for a pair | Set RMSD = ∞, exclude from averages. Print warning listing every skipped pair. |
| All pairs fail for one instance | That instance gets avg RMSD = ∞, cannot be medoid. Warn user. |
| Tied medoids | Pick first; note the tie in the report. |
| Motif name mismatch across sources | Alias table resolves variants (see Motif Name Aliasing). |

### Large Sets (N > 100)

For N = 100: $\binom{100}{2} = 4950$ superimpositions ≈ 50 seconds + object
creation time.

**Strategy:**
1. Print warning: `"100 instances → 4950 pairwise comparisons. This may take ~1 min. Continue? [Y/n]"`
2. Show continuous progress: `Computing RMSD... 500/4950`
3. If N > 200, suggest narrowing with instance indices or `pdb=` filter.
4. **Future optimisation:** Subsample K random instances for approximate medoid
   selection, then superimpose all N onto it exactly.

---

## Files to Create / Modify

| File | Action |
|------|--------|
| `rsmviewer/alignment.py` | **REWRITE** — Medoid pipeline, alias table, `register_alignment_commands()` for both `rmv_super` + `rmv_align` |
| `rsmviewer/gui.py` | **EDIT** — Add `'rmv_super'`, `'rmv_align'` to `_RMV_COMMANDS` and `_CMD_SUFFIXES`, call `register_alignment_commands()`, add help section |
| `rsmviewer/plugin.py` | **EDIT** — Add commands to Quick Start if desired |
| `README.md` | **EDIT** — Add `rmv_super` / `rmv_align` to command table |
| `TUTORIAL.md` | **EDIT** — Add usage example |

---

## Dependencies

- **PyMOL** (`cmd.super`, `cmd.align`, `cmd.create`, `cmd.delete`,
  `cmd.get_object_list`, `cmd.color`, `cmd.feedback`, `cmd.rebuild`) —
  already available
- **No external libraries** needed. Pure Python nested lists for the
  distance matrix.

---

## Testing Strategy

1. **Unit test `find_medoid()`** — known distance matrix → verify correct index.
2. **Single-PDB** — `rmv_fetch 1S72`, `rmv_show KTURN`, `rmv_super KTURN`.
   Verify: objects created, uniquely coloured, medoid green, RMSD table printed.
3. **`rmv_align`** — same as above but uses `cmd.align()` internally.
4. **Cross-PDB** — fetch 1S72 + 4V88, `rmv_super KTURN`. Verify instances
   from both PDBs included.
5. **PDB filter** — `rmv_super KTURN pdb=1S72`. Verify only 1S72 instances.
6. **Multi-PDB filter** — `rmv_super KTURN pdb=1S72,4V88,4V9F`.
7. **Alias resolution** — load RMSX (K-TURN), type `rmv_super KTURN`. Verify it
   resolves to K-TURN instances.
8. **Failed pair** — verify warning printed, pair excluded from averages.
9. **Edge cases** — 1 instance, 2 instances, non-existent motif type.

---

## Summary

| Aspect | Detail |
|--------|--------|
| Commands | `rmv_super` (sequence-independent) / `rmv_align` (sequence-dependent) |
| Syntax | `rmv_super <MOTIF_TYPE> [indices] [pdb=X,Y,Z]` |
| File | `rsmviewer/alignment.py` (repurposed) |
| Algorithm | Pairwise RMSD → medoid (min avg) → superimpose all onto medoid |
| Performance | Batch object creation with suppressed GUI updates; progress indicators |
| Colouring | Each instance unique colour; medoid highlighted green |
| Cross-PDB | Works automatically; optional `pdb=` filter for scoping (any number of PDBs) |
| Name aliasing | Alias table + fuzzy matching resolves motif name variants across sources |
| Failed pairs | Excluded from averages, warning printed with details |
| Large sets | Warning + progress bar; suggest filtering if N > 200 |
| Dependencies | None (pure PyMOL + Python) |
