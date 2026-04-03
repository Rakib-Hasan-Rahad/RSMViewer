# RSMViewer — Multi-Source Combine Guide

A detailed walkthrough of the multi-source combine pipeline: how data from
different sources is fetched, enriched, deduplicated, and merged into a single
unified motif set.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Merge Rules — Complete Reference](#2-merge-rules--complete-reference)
3. [Pipeline Steps at a Glance](#3-pipeline-steps-at-a-glance)
4. [Example 1 — NoBIAS + RMSX (Source 8 + 7)](#4-example-1--nobias--rmsx-source-8--7)
5. [Example 2 — Rfam + RMSX (Source 2 + 7) — Numbered Suffixes](#5-example-2--rfam--rmsx-source-2--7--numbered-suffixes)
6. [Example 3 — Atlas + BGSU (Source 1 + 3)](#6-example-3--atlas--bgsu-source-1--3)
7. [Example 4 — Three Sources (Source 1 + 3 + 7)](#7-example-4--three-sources-source-1--3--7)
8. [Detailed Code Walkthrough](#8-detailed-code-walkthrough)
9. [Key Data Structures](#9-key-data-structures)
10. [FAQ / Troubleshooting](#10-faq--troubleshooting)

---

## 1. Overview

When you run:

```
rmv_db 8 7
rmv_fetch 1S72
rmv_load_motif
```

RSMViewer enters **combine mode**. Instead of loading motifs from a single
source, it:

1. Fetches data from **each** source independently
2. Enriches generic names (for sources 1, 3, 5) via homolog lookup
3. Stamps each instance with its source origin (source ID + label)
4. Removes exact duplicates **within** each source
5. **Category-aware cascade merge** — compares motifs only within the same
   biological category, using Jaccard similarity to detect overlaps

The result is a single unified motif set where overlapping entries are
resolved by **priority order** (left = highest, right = lowest), while
biologically distinct motif types are always preserved.

### Priority Order

The order you type the source IDs determines priority:

```
rmv_db 8 7      # Source 8 = HIGHEST priority, Source 7 = LOWEST
rmv_db 7 8      # Source 7 = HIGHEST priority, Source 8 = LOWEST
rmv_db 1 3 7    # Source 1 > Source 3 > Source 7
```

If the same motif region is found in both sources under the same category,
the higher-priority version is kept and the lower-priority duplicate is
discarded.

### Priority Also Controls Naming Convention

The **priority source's naming convention** is preserved in the final output.
If source 2 (Rfam) calls a motif `SARCIN_RICIN_1` and source 7 (RMSX) calls
the same region `SARCIN-RICIN`, then:

- `rmv_db 2 7` → output key is `SARCIN_RICIN_1` (Rfam's name, highest priority)
- `rmv_db 7 2` → output key is `SARCIN-RICIN` (RMSX's name, highest priority)

### Configurable Jaccard Similarity Threshold (Optional)

You can optionally adjust the **Jaccard similarity threshold** for overlap
detection. This controls how strict the merger is when deciding if two motifs
are "the same" region.

Default: **0.60** (60%)

```
rmv_db 8 7, jaccard_threshold=0.80   # Stricter: only very similar motifs merge
rmv_db 1 3, jaccard_threshold=0.50   # Looser: more likely to merge nearby regions
rmv_db 2 7                            # Uses default 0.60
```

Supported formats:
- `jaccard_threshold=0.80` — decimal (0.0 to 1.0)
- `jaccard_threshold=80` — percentage (values > 1 treated as %, so 80 → 0.80)
- `jaccard_threshold=80%` — with % sign

Higher thresholds (e.g., 0.80) are **stricter** → fewer motifs merge → more
instances kept. Lower thresholds (e.g., 0.40) are **looser** → more motifs
merge → fewer instances kept.

---

## 2. Merge Rules — Complete Reference

This section lists every rule the cascade merger follows, in order of
application.

### Rule 1 — Source Priority (Left = Highest)

The first source ID in the `rmv_db` command has the highest priority. When
two motifs overlap, the higher-priority source's version is kept.

```
rmv_db 2 7     →  Source 2 (Rfam) wins ties over Source 7 (RMSX)
rmv_db 7 2     →  Source 7 (RMSX) wins ties over Source 2 (Rfam)
```

### Rule 2 — Right-to-Left Cascade

Merging proceeds from right to left. With three sources `[A, B, C]`:

```
Step 1: intermediate = pairwise_merge(B, C)    # B has priority over C
Step 2: final       = pairwise_merge(A, intermediate)  # A has priority over all
```

Each pairwise merge: all ref (higher-priority) motifs are kept unconditionally;
each updater (lower-priority) motif is tested against ref motifs of the same
category.

### Rule 3 — Category-Aware Comparison (Same Type Only)

Overlap is checked **only between motifs of the same base category**. A
`K-TURN` will never be discarded because it spatially overlaps a `C-LOOP`
from a higher-priority source. Both are preserved because they represent
different biological structures.

**Why:** Different motif types annotate different structural roles. Even if
they happen to share residues in 3D space, discarding one would lose the
biological information that those residues participate in both roles.

### Rule 4 — Base Category Extraction

Before comparison, every motif type name is reduced to a **base category**
using two steps:

**Step A — Strip trailing numeric suffix (`_N` or `-N`):**

Some tools append sub-family numbers (e.g., `SARCIN_RICIN_1`,
`SARCIN_RICIN_2`, `K_TURN_1`, `K_TURN_2`). These all represent the same
biological family, so the suffix is stripped before category comparison.

The regex `[-_]\d+$` matches a hyphen or underscore followed by one or more
digits at the end of the name. It requires a separator (`-` or `_`), so
names like `J3` or `J7` are **not** affected (no separator before the digit).

```
SARCIN_RICIN_1   → strip "_1"  → SARCIN_RICIN
SARCIN_RICIN_2   → strip "_2"  → SARCIN_RICIN
K_TURN_1         → strip "_1"  → K_TURN
K_TURN_2         → strip "_2"  → K_TURN
SARCIN-RICIN     → no suffix   → SARCIN-RICIN
K-TURN           → no suffix   → K-TURN
GNRA             → no suffix   → GNRA
J3               → no match    → J3  (no separator before digit)
J7               → no match    → J7
```

**Step B — Harmonize via `MOTIF_CANONICAL_MAP`:**

After stripping, variant spellings are mapped to a canonical form:

```
SARCIN_RICIN  → SARCIN-RICIN    (canonical)
K_TURN        → K-TURN           (canonical)
KINK-TURN     → K-TURN
KINK_TURN     → K-TURN
CLOOP         → C-LOOP
C_LOOP        → C-LOOP
REVERSE_KTURN → REVERSE-K-TURN
E_LOOP        → E-LOOP
T_LOOP        → T-LOOP
GNRA          → GNRA             (no mapping = unchanged, uppercased)
```

**Combined examples:**

| Original Name     | After Strip | After Harmonize | Base Category |
|-------------------|-------------|-----------------|---------------|
| `SARCIN_RICIN_1`  | `SARCIN_RICIN` | `SARCIN-RICIN` | **SARCIN-RICIN** |
| `SARCIN_RICIN_2`  | `SARCIN_RICIN` | `SARCIN-RICIN` | **SARCIN-RICIN** |
| `SARCIN-RICIN`    | `SARCIN-RICIN` | `SARCIN-RICIN` | **SARCIN-RICIN** |
| `K_TURN_1`        | `K_TURN`       | `K-TURN`       | **K-TURN**       |
| `K_TURN_2`        | `K_TURN`       | `K-TURN`       | **K-TURN**       |
| `K-TURN`          | `K-TURN`       | `K-TURN`       | **K-TURN**       |
| `KINK-TURN`       | `KINK-TURN`    | `K-TURN`       | **K-TURN**       |
| `T_LOOP`          | `T_LOOP`       | `T-LOOP`       | **T-LOOP**       |
| `GNRA`            | `GNRA`         | `GNRA`         | **GNRA**         |
| `J3`              | `J3`           | `J3`           | **J3**           |

### Rule 5 — Jaccard Similarity Threshold (User-Configurable)

Two motif instances are compared by their residue sets — sets of
`(chain, residue_number)` pairs:

```
J(A, B) = |A ∩ B| / |A ∪ B|
```

- **J = 1.0** → identical residue sets
- **J = 0.0** → no shared residues (completely different locations)
- **J ≥ threshold** → considered overlapping → lower-priority version is **discarded**
- **J < threshold** → considered independent → both are **kept**

#### Default Threshold

The **default threshold is 0.60** (60%). Set in the `MotifVisualizerGUI.__init__()`:

```python
self.jaccard_threshold = 0.60  # Instance variable on MotifVisualizerGUI
```

Used during merge:

```python
merger = CascadeMerger(jaccard_threshold=self.jaccard_threshold)
```

#### User Override

You can override the threshold on the command line:

```
rmv_db 8 7, jaccard_threshold=0.80   # One-time override for this combine
rmv_db 1 3, jaccard_threshold=75     # Looser (75%)
```

The new threshold **persists** for subsequent combines until explicitly changed
to a different value. To reset to default:

```
rmv_db 1 3, jaccard_threshold=0.60   # Reset to default
```

### Rule 6 — Original Keys Preserved

The cascade merger does **not** rename motif-type keys. Each source's motifs
enter (and survive) the merge under their **original** uppercased key.

In a pairwise merge of ref (priority) vs. updater:
- **Ref motifs** are added under their own keys unconditionally.
- **Surviving updater motifs** (no overlap) are added under their own keys.

This means the final output may contain **multiple keys** that map to the
same base category when contributions come from different sources:

```
rmv_db 2 7  with 1S72

Source 2 (Rfam, priority):    SARCIN_RICIN_1 (2 instances), SARCIN_RICIN_2 (3 instances)
Source 7 (RMSX):              SARCIN-RICIN (11 instances)

After merge — keys in output:
  SARCIN_RICIN_1:  2 from Rfam (all kept, highest priority)
  SARCIN_RICIN_2:  3 from Rfam (all kept, highest priority)
  SARCIN-RICIN:    X from RMSX (only those NOT overlapping any Rfam SARCIN-RICIN)
```

### Rule 7 — Source Attribution (Overlap Annotation)

When a lower-priority motif is discarded due to overlap, the surviving
higher-priority instance is annotated with `_also_found_in` metadata. This
lets the display layer show cross-source confirmation (e.g., "NoBIAS + RMSX").

For three or more sources, these annotations are propagated through the
cascade: if source C's motif is discarded in favor of B's, and then B's is
discarded in favor of A's, the final surviving instance from A carries
attribution from both B and C.

### Rule 8 — Within-Source Exact Deduplication (Before Merge)

Before the cascade merge, each source is deduplicated independently using
**exact** residue-set matching (frozenset equality). Two instances with the
identical `{(chain, residue_number)}` set within the same source and same
motif type are considered duplicates — only the first is kept.

This differs from the cascade merge, which uses **fuzzy** Jaccard matching
across sources.

### Rule 9 — Enrichment for Generic Sources (Before Merge)

Sources 1 (Atlas), 3 (BGSU), and 5 (FR3D) produce **generic** motif-type
names like `HL`, `IL`, `J3`. Before merging, these are enriched via homolog
representative lookup to obtain specific names like `GNRA`, `K-TURN`, etc.

This enrichment runs **before** source stamping and deduplication, so by the
time the cascade merge runs, all sources have their best possible names.

### Summary — How It All Fits Together

```
For each updater motif instance:
  1. Compute base_category = strip_suffix(name) → harmonize(name)
  2. Find all ref motifs whose base_category matches
  3. If no ref motifs share that category → KEEP (unique category)
  4. Compute Jaccard against each same-category ref motif
  5. If Jaccard ≥ 0.60 with ANY → DISCARD (annotate surviving ref)
  6. If Jaccard < 0.60 with ALL → KEEP (unique instance within category)
```

---

## 3. Pipeline Steps at a Glance

```
User: rmv_db 2 7  →  rmv_load_motif

  ┌──────────────────────────────────────────────────────┐
  │ Step 1: Fetch raw motifs from each source            │
  │   _fetch_from_single_source(pdb, 2)  → Rfam data    │
  │   _fetch_from_single_source(pdb, 7)  → RMSX data    │
  ├──────────────────────────────────────────────────────┤
  │ Step 2: Enrich generic sources (1, 3, 5 only)       │
  │   Sources 2 and 7 are NOT generic → skip            │
  ├──────────────────────────────────────────────────────┤
  │ Step 2.5: Source Stamping                            │
  │   Tag every instance with _source_id + _source_label│
  ├──────────────────────────────────────────────────────┤
  │ Step 2.7: Within-Source Deduplication                │
  │   Remove exact duplicate residue sets per source     │
  ├──────────────────────────────────────────────────────┤
  │ Step 3: Cascade Merge (category-aware)               │
  │   Right-to-left: merge(Src2, Src7)                   │
  │   Base category grouping → Jaccard within category   │
  │   Jaccard ≥ 0.60 → discard lower-priority duplicate  │
  │   Original keys preserved from each source           │
  └──────────────────────────────────────────────────────┘

  Result: unified dict {motif_type: [MotifInstance, ...]}
```

---

## 4. Example 1 — NoBIAS + RMSX (Source 8 + 7)

The simplest combine scenario: two user-annotation sources with specific
motif names and identical naming conventions.

### 4.1 Setup

```python
rmv_db 8 /path/to/nobias/output
rmv_db 7 /path/to/rmsx/output
rmv_db 8 7                            # Combine: 8 (high) + 7 (low)
rmv_fetch 1S72
rmv_load_motif
```

### 4.2 What Happens

Both sources use the same naming convention (`K-TURN`, `SARCIN-RICIN`, etc.),
so base category extraction is straightforward — names map to themselves.

**Step 1 — Fetch:**
```
NoBIAS (8): K-TURN ×7, SARCIN-RICIN ×2
RMSX   (7): K-TURN ×5, C-LOOP ×2, SARCIN-RICIN ×3
```

**Steps 2–2.7 — Enrichment skipped, stamping, dedup.**

**Step 3 — Cascade merge (pairwise: NoBIAS ref, RMSX updater):**

All 9 NoBIAS motifs go into the merged set unconditionally.

For each RMSX motif:

| RMSX Instance | Base Category | Same-Cat Ref Motifs | Best Jaccard | Decision |
|--------------|---------------|---------------------|--------------|----------|
| K-TURN #1    | K-TURN        | NoBIAS K-TURN ×7    | 0.85         | DISCARD (overlap with NoBIAS) |
| K-TURN #2    | K-TURN        | NoBIAS K-TURN ×7    | 0.78         | DISCARD |
| K-TURN #3    | K-TURN        | NoBIAS K-TURN ×7    | 0.92         | DISCARD |
| K-TURN #4    | K-TURN        | NoBIAS K-TURN ×7    | 0.80         | DISCARD |
| K-TURN #5    | K-TURN        | NoBIAS K-TURN ×7    | 0.10         | **KEEP** (unique location) |
| C-LOOP #1    | C-LOOP        | *none* (NoBIAS has no C-LOOP) | — | **KEEP** (unique category) |
| C-LOOP #2    | C-LOOP        | *none*              | —            | **KEEP** |
| SARCIN-RICIN #1 | SARCIN-RICIN | NoBIAS SARCIN-RICIN ×2 | 0.88  | DISCARD |
| SARCIN-RICIN #2 | SARCIN-RICIN | NoBIAS SARCIN-RICIN ×2 | 0.75  | DISCARD |
| SARCIN-RICIN #3 | SARCIN-RICIN | NoBIAS SARCIN-RICIN ×2 | 0.05  | **KEEP** (unique location) |

**Final result:** 13 motifs in 3 categories
- K-TURN: 7 (NoBIAS) + 1 (RMSX unique) = **8**
- SARCIN-RICIN: 2 (NoBIAS) + 1 (RMSX unique) = **3**
- C-LOOP: 0 (NoBIAS) + 2 (RMSX unique) = **2**

> **Key point:** The two RMSX C-LOOP instances were **never compared** against
> NoBIAS K-TURN or SARCIN-RICIN — they are a different base category and are
> kept unconditionally.

---

## 5. Example 2 — Rfam + RMSX (Source 2 + 7) — Numbered Suffixes

This example shows how **trailing numeric suffixes** are handled: sources
that use sub-family numbering (like Rfam's `SARCIN_RICIN_1`, `K_TURN_2`)
are correctly compared against sources that use flat names (like RMSX's
`SARCIN-RICIN`, `K-TURN`).

### 5.1 Setup

```python
rmv_db 2 7          # Rfam (high priority) + RMSX (low priority)
rmv_fetch 1S72
rmv_load_motif
```

### 5.2 Source Data

```
Rfam (Source 2):
  GNRA                 14
  T_LOOP                6
  SARCIN_RICIN_2        3
  C_LOOP                3
  U_TURN                3
  SARCIN_RICIN_1        2
  UNCG                  1
  K_TURN_1              1
  TWIST_UP              1
  K_TURN_2              1
  ──────────────────────
  Total: 35 motifs in 10 categories

RMSX (Source 7):
  REVERSE-K-TURN       12
  SARCIN-RICIN         11
  K-TURN                9
  E-LOOP                3
  C-LOOP                2
  ──────────────────────
  Total: 37 motifs in 5 categories
```

### 5.3 Base Category Extraction

Before comparison, every motif type is reduced to its base category:

| Source | Original Key       | Strip Suffix | Harmonize       | Base Category       |
|--------|--------------------|--------------|-----------------|---------------------|
| Rfam   | `SARCIN_RICIN_1`   | `SARCIN_RICIN` | `SARCIN-RICIN` | **SARCIN-RICIN**   |
| Rfam   | `SARCIN_RICIN_2`   | `SARCIN_RICIN` | `SARCIN-RICIN` | **SARCIN-RICIN**   |
| RMSX   | `SARCIN-RICIN`     | `SARCIN-RICIN` | `SARCIN-RICIN` | **SARCIN-RICIN**   |
| Rfam   | `K_TURN_1`         | `K_TURN`       | `K-TURN`       | **K-TURN**         |
| Rfam   | `K_TURN_2`         | `K_TURN`       | `K-TURN`       | **K-TURN**         |
| RMSX   | `K-TURN`           | `K-TURN`       | `K-TURN`       | **K-TURN**         |
| Rfam   | `C_LOOP`           | `C_LOOP`       | `C-LOOP`       | **C-LOOP**         |
| RMSX   | `C-LOOP`           | `C-LOOP`       | `C-LOOP`       | **C-LOOP**         |
| Rfam   | `T_LOOP`           | `T_LOOP`       | `T-LOOP`       | **T-LOOP**         |
| Rfam   | `GNRA`             | `GNRA`         | `GNRA`         | **GNRA**           |
| RMSX   | `E-LOOP`           | `E-LOOP`       | `E-LOOP`       | **E-LOOP**         |
| RMSX   | `REVERSE-K-TURN`   | `REVERSE-K-TURN` | `REVERSE-K-TURN` | **REVERSE-K-TURN** |

### 5.4 Cascade Merge (Rfam ref, RMSX updater)

All 35 Rfam motifs go into the merged set under their original keys.

For each RMSX motif, only **same-base-category** ref motifs are considered:

**RMSX `SARCIN-RICIN` (11 instances) vs Rfam base category `SARCIN-RICIN`:**

The Rfam ref pool for category SARCIN-RICIN includes:
- `SARCIN_RICIN_1` (2 instances)
- `SARCIN_RICIN_2` (3 instances)

→ Total 5 ref instances to compare against.

| RMSX Instance      | Compared Against               | Best Jaccard | Decision |
|--------------------|---------------------------------|--------------|----------|
| SARCIN-RICIN #1    | Rfam SARCIN_RICIN_1 + _2 (5)   | 0.90         | DISCARD  |
| SARCIN-RICIN #2    | same pool                       | 0.85         | DISCARD  |
| SARCIN-RICIN #3    | same pool                       | 0.78         | DISCARD  |
| SARCIN-RICIN #4    | same pool                       | 0.65         | DISCARD  |
| SARCIN-RICIN #5    | same pool                       | 0.72         | DISCARD  |
| ...                | ...                             | ...          | ...      |
| SARCIN-RICIN #9    | same pool                       | 0.15         | **KEEP** |
| SARCIN-RICIN #10   | same pool                       | 0.08         | **KEEP** |
| SARCIN-RICIN #11   | same pool                       | 0.20         | **KEEP** |

**RMSX `K-TURN` (9 instances) vs Rfam base category `K-TURN`:**

The Rfam ref pool for K-TURN includes:
- `K_TURN_1` (1 instance)
- `K_TURN_2` (1 instance)

→ Total 2 ref instances to compare against.

| RMSX Instance | Compared Against           | Best Jaccard | Decision |
|--------------|----------------------------|--------------|----------|
| K-TURN #1    | Rfam K_TURN_1 + _2 (2)     | 0.88         | DISCARD  |
| K-TURN #2    | same pool                  | 0.75         | DISCARD  |
| K-TURN #3    | same pool                  | 0.10         | **KEEP** |
| ...          | ...                        | ...          | ...      |
| K-TURN #9    | same pool                  | 0.05         | **KEEP** |

**RMSX `C-LOOP` (2 instances) vs Rfam base category `C-LOOP`:**

Rfam has `C_LOOP` (3 instances) → base = C-LOOP.

| RMSX Instance | Compared Against       | Best Jaccard | Decision |
|--------------|------------------------|--------------|----------|
| C-LOOP #1    | Rfam C_LOOP (3)        | 0.82         | DISCARD  |
| C-LOOP #2    | same pool              | 0.70         | DISCARD  |

**RMSX `E-LOOP` (3 instances) vs Rfam base category `E-LOOP`:**

Rfam has no E-LOOP → no same-category ref motifs → all 3 are **KEPT**.

**RMSX `REVERSE-K-TURN` (12 instances) vs Rfam base category `REVERSE-K-TURN`:**

Rfam has no REVERSE-K-TURN → all 12 are **KEPT**.

### 5.5 Final Merged Output

```
From Rfam (priority, preserved as-is):
  GNRA                 14
  T_LOOP                6
  SARCIN_RICIN_2        3
  C_LOOP                3
  U_TURN                3
  SARCIN_RICIN_1        2
  UNCG                  1
  K_TURN_1              1
  K_TURN_2              1
  TWIST_UP              1

From RMSX (survivors, under RMSX's original keys):
  REVERSE-K-TURN       12   (unique category — all kept)
  SARCIN-RICIN          3   (unique locations not in Rfam)
  K-TURN                7   (unique locations not in Rfam)
  E-LOOP                3   (unique category — all kept)
```

> **Note the naming:** The output has BOTH `SARCIN_RICIN_1` (from Rfam) AND
> `SARCIN-RICIN` (from RMSX) as separate keys, because original keys are
> preserved. They were compared against each other during the merge (same
> base category), but surviving entries keep their source's name.

### 5.6 Flipping Priority: `rmv_db 7 2`

If you flip the order, RMSX is priority and Rfam is the updater:

```python
rmv_db 7 2          # RMSX (high priority) + Rfam (low priority)
```

Now **all 37 RMSX motifs** go into the merged set first. Then each of
Rfam's 35 motifs is checked:

- Rfam `SARCIN_RICIN_1` (2 instances) are compared against RMSX
  `SARCIN-RICIN` (11 instances) — same base category.
- Rfam `K_TURN_1` (1 instance) is compared against RMSX `K-TURN` (9
  instances) — same base category.
- Rfam `GNRA` (14 instances) has no RMSX counterpart → all **KEPT**.

The output keys would be `SARCIN-RICIN`, `K-TURN` (from RMSX) plus `GNRA`,
`T_LOOP`, `SARCIN_RICIN_1` survivors (from Rfam, if any don't overlap).

---

## 6. Example 3 — Atlas + BGSU (Source 1 + 3)

This scenario shows **homolog enrichment**, which is the key difference from
the user-source examples above.

### 6.1 Setup

```python
rmv_db 1 3          # Atlas (high priority) + BGSU (low priority)
rmv_fetch 1S72
rmv_load_motif
```

### 6.2 What Happens Inside

#### Step 1 — Fetch

```python
raw_sources[1] = atlas_provider.get_motifs_for_pdb('1S72')
# Returns generic names: {'HL': [...], 'IL': [...], 'J3': [...]}

raw_sources[3] = bgsu_api_provider.get_motifs_for_pdb('1S72')
# Returns: {'GNRA': [...], 'K-TURN': [...], 'HL': [...], ...}
```

Atlas returns **generic** names (`HL`, `IL`, `J3`, etc.) because its bundled
JSON data doesn't include semantic annotations. BGSU API returns a mix of
generic and specific names.

#### Step 2 — Enrichment

Both sources 1 and 3 are in `GENERIC_SOURCES = {1, 3, 5}`, so both get
enriched via the HomologEnricher:

```
1. Query: 1S72 → NR list → Representative: 4V9F
2. Fetch 4V9F annotations from BGSU (has full semantic annotations)
3. Build lookup: motif_group_id → specific_category
     e.g., HL_75660.8 → GNRA
           HL_75661.2 → T-LOOP
4. For each generic instance in 1S72:
     - Get its motif_group from metadata
     - Look up in representative's table → "GNRA"
     - Rename: HL → GNRA (residues unchanged)
```

After enrichment, both sources have specific names where possible. Whatever
remains generic (`HL`, `IL`) could not be resolved via the representative.

#### Steps 2.5–2.7 — Source stamping + within-source dedup.

#### Step 3 — Cascade Merge

```python
merged = merger.merge_sources([Atlas, BGSU], ['RNA 3D Atlas', 'BGSU RNA 3D Hub'])
```

Category-aware comparison: Atlas `GNRA` instances are compared against BGSU
`GNRA` instances (same base category). Atlas `HL` against BGSU `HL`.
Different categories are never compared — BGSU `K-TURN` is not compared
against Atlas `HL`, even if they share residues.

### 6.3 Why This Matters

Atlas is bundled and works offline, but it has limited coverage (759 PDBs)
and uses generic names. BGSU has much broader coverage (~3000+ PDBs) and
sometimes provides additional motifs that Atlas misses. Combining them gives
you the best of both worlds:

- Atlas motifs are used where available (they're usually validated)
- BGSU fills in gaps with additional motifs not present in Atlas
- Enrichment ensures generic HL/IL names are replaced with GNRA/K-TURN etc.

---

## 7. Example 4 — Three Sources (Source 1 + 3 + 7)

This demonstrates the **right-to-left cascade** with more than two sources.

### 7.1 Setup

```python
rmv_db 7 /path/to/rmsx/data
rmv_db 1 3 7                    # Atlas > BGSU > RMSX (priority order)
rmv_fetch 1S72
rmv_load_motif
```

### 7.2 Cascade Order

With three sources `[1, 3, 7]`, the merge proceeds **right-to-left**:

```
Step A: intermediate = pairwise_merge(Src3_BGSU, Src7_RMSX)
            ref = BGSU (priority)     updater = RMSX (lower)
            → intermediate result

Step B: final = pairwise_merge(Src1_Atlas, intermediate)
            ref = Atlas (priority)    updater = intermediate
            → final merged result
```

### 7.3 Concrete Trace

Suppose:

| Source | Motifs |
|--------|--------|
| Atlas (1) | GNRA ×2, K-TURN ×1, HL ×3 |
| BGSU (3) | GNRA ×2, SARCIN-RICIN ×1, HL ×5 |
| RMSX (7) | K-TURN ×2, C-LOOP ×1, SARCIN-RICIN ×1 |

**Step A — Merge BGSU (ref) + RMSX (updater):**

| RMSX Motif | Base Category | Same-Cat Refs in BGSU | Best Jaccard | Decision |
|-----------|---------------|----------------------|--------------|----------|
| K-TURN #1 | K-TURN | *none* (BGSU has no K-TURN) | — | **KEEP** |
| K-TURN #2 | K-TURN | *none* | — | **KEEP** |
| C-LOOP #1 | C-LOOP | *none* (BGSU has no C-LOOP) | — | **KEEP** |
| SARCIN-RICIN #1 | SARCIN-RICIN | BGSU SARCIN-RICIN ×1 | 0.85 | DISCARD |

Intermediate = BGSU motifs + RMSX K-TURN ×2 + RMSX C-LOOP ×1
(10 motifs total, RMSX SARCIN-RICIN #1 was discarded)

**Step B — Merge Atlas (ref) + intermediate (updater):**

| Intermediate Motif | Base Category | Same-Cat Atlas Refs | Best Jaccard | Decision |
|-------------------|---------------|---------------------|--------------|----------|
| BGSU GNRA #1 | GNRA | Atlas GNRA ×2 | 0.92 | DISCARD |
| BGSU GNRA #2 | GNRA | Atlas GNRA ×2 | 0.88 | DISCARD |
| BGSU SARCIN-RICIN #1 | SARCIN-RICIN | *none* (Atlas has no SARCIN-RICIN) | — | **KEEP** |
| BGSU HL #1 | HL | Atlas HL ×3 | 0.95 | DISCARD |
| BGSU HL #2 | HL | Atlas HL ×3 | 0.90 | DISCARD |
| BGSU HL #3 | HL | Atlas HL ×3 | 0.87 | DISCARD |
| BGSU HL #4 | HL | Atlas HL ×3 | 0.02 | **KEEP** |
| BGSU HL #5 | HL | Atlas HL ×3 | 0.05 | **KEEP** |
| RMSX K-TURN #1 | K-TURN | Atlas K-TURN ×1 | 0.78 | DISCARD |
| RMSX K-TURN #2 | K-TURN | Atlas K-TURN ×1 | 0.10 | **KEEP** |
| RMSX C-LOOP #1 | C-LOOP | *none* (Atlas has no C-LOOP) | — | **KEEP** |

**Final result:**
Atlas originals (GNRA ×2, K-TURN ×1, HL ×3) + survivors (BGSU SARCIN-RICIN,
BGSU HL ×2, RMSX K-TURN, RMSX C-LOOP) = **10 motifs** from 3 sources.

> **Key difference from old behavior:** In the old cross-type system, RMSX
> C-LOOP #1 could have been discarded if it happened to share residues with
> an Atlas HL instance. Now it is unconditionally kept because C-LOOP and HL
> are different base categories.

---

## 8. Detailed Code Walkthrough

### 8.1 Entry Point — `_handle_multi_source()` (gui.py)

**File:** `rsmviewer/gui.py`

Called when `rmv_db` receives more than one source ID. Sets up combine mode:

```python
self.combined_source_ids = source_ids       # e.g., [8, 7]
self.current_source_mode = 'combine'
self.current_source_id = '_'.join(str(s) for s in source_ids)  # "8_7"
```

The `current_source_id` string is used as a suffix for PyMOL object names
(e.g., `HL_ALL_S_8_7`) so objects from combine mode are visually distinct.

### 8.2 Fetching — `_fetch_from_single_source()` (gui.py)

**File:** `rsmviewer/gui.py`

Routes each source ID to the correct provider:

| Source Type | Provider |
|------------|----------|
| `local` (1, 2) | `source_selector.providers['atlas']` or `['rfam']` |
| `web` (3, 4) | `source_selector.providers['bgsu_api']` or `['rfam_api']` |
| `user` (5-8) | `UserAnnotationProvider` with tool-specific config |

For user sources, per-source custom paths are applied:

```python
_udp = self.user_data_paths.get(source_id)
if _udp:
    tool_name_map = {
        'rmsx': 'RNAMotifScanX',
        'rms': 'RNAMotifScan',
        'fr3d': 'fr3d',
        'nobias': 'NoBIAS',
    }
    internal_tool = tool_name_map.get(tool.lower(), tool)
    provider.override_tool_dirs[internal_tool] = Path(_udp)
```

> **Per-source paths:** Each source stores its custom path independently in
> `self.user_data_paths` (a dict keyed by source ID int). Setting a path for
> source 8 does NOT affect source 7's path.

### 8.3 Enrichment — `HomologEnricher.enrich()` (homolog_enricher.py)

**File:** `rsmviewer/database/homolog_enricher.py`

Only runs for sources in `GENERIC_SOURCES = {1, 3, 5}` (Atlas, BGSU, FR3D).
These sources return generic names like `HL`, `IL`, `J3`.

The enricher uses the BGSU Non-Redundant (NR) representative set:

```
Query PDB (e.g., 1S72)
  → NR list lookup
  → Representative PDB (e.g., 4V9F)
  → Fetch 4V9F's full annotations from BGSU
  → Build motif_group → annotation lookup
  → For each generic instance in 1S72:
      Check its motif_group ID
      If the representative has a specific name for that group: rename
```

If a match is found, a new `MotifInstance` is created with:
- The specific name (e.g., `GNRA`)
- The original residues (unchanged)
- Extra metadata: `enriched_from: "HL"`, `enrichment_source: "homolog"`

### 8.4 Source Stamping — Step 2.5 (gui.py)

**File:** `rsmviewer/gui.py`

Every instance is tagged before merging:

```python
for sid in source_ids:
    label = SOURCE_ID_MAP[sid]['name']
    for mtype, instances in raw_sources[sid].items():
        for inst in instances:
            if inst.metadata is None:
                inst.metadata = {}
            inst.metadata['_source_id'] = sid
            inst.metadata['_source_label'] = label
```

This metadata persists through the cascade merge and is used by:
- `rmv_show TYPE N` — shows `SOURCE` column in instance tables
- `rmv_summary` — prints source attribution report

### 8.5 Within-Source Deduplication — Step 2.7 (gui.py)

**File:** `rsmviewer/gui.py`

Removes exact duplicates **within** a single source. Two instances are
considered duplicates if they have the identical set of `(chain, residue_number)`
pairs.

**How it differs from the cascade merge:**
- **Within-source dedup**: exact match (frozenset equality) — same residues = duplicate
- **Cascade merge**: fuzzy match (Jaccard ≥ 0.60) — overlapping but not identical

### 8.6 Cascade Merge — `CascadeMerger` (cascade_merger.py)

**File:** `rsmviewer/database/cascade_merger.py`

#### Algorithm

```
Given sources [S1, S2, S3] (S1 = highest priority):
  result = pairwise_merge(S2, S3)    # S2 wins ties
  final  = pairwise_merge(S1, result) # S1 wins ties
```

#### Base Category Extraction (`_get_base_category`)

```python
_TRAILING_NUM_RE = re.compile(r'[-_]\d+$')

def _get_base_category(name: str) -> str:
    key = name.strip().upper()
    stripped = _TRAILING_NUM_RE.sub('', key)       # Remove trailing _1, -2, etc.
    return MOTIF_CANONICAL_MAP.get(stripped, stripped)  # Harmonize spelling
```

#### Category-Aware Pairwise Merge

The core merge function groups ref motifs by base category and only compares
updater motifs against same-category refs:

```python
def _pairwise_merge(self, ref, updater, ...):
    # 1. All ref motifs go into merged under their original (uppercased) keys
    merged = {}
    for mtype, instances in ref.items():
        merged[mtype.upper()] = list(instances)

    # 2. Group ref residue sets by BASE CATEGORY
    ref_by_cat = {}
    for mtype, instances in ref.items():
        cat = _get_base_category(mtype)          # e.g., SARCIN_RICIN_1 → SARCIN-RICIN
        for inst in instances:
            rset = _get_residue_set(inst)
            ref_by_cat.setdefault(cat, []).append((mtype.upper(), inst, rset))

    # 3. For each updater motif, check same-category refs ONLY
    for u_type, u_instances in updater.items():
        u_cat = _get_base_category(u_type)       # e.g., SARCIN-RICIN → SARCIN-RICIN
        same_cat_refs = ref_by_cat.get(u_cat, [])  # Only same-category!

        for u_inst in u_instances:
            u_rset = _get_residue_set(u_inst)

            is_overlap = False
            for _r_key, r_inst, r_rset in same_cat_refs:
                if _jaccard(u_rset, r_rset) >= 0.60:
                    is_overlap = True
                    _annotate_overlap(r_inst, u_inst)  # Mark cross-source support
                    break

            if not is_overlap:
                # Unique — add under updater's original key
                merged.setdefault(u_type.upper(), []).append(u_inst)
```

#### Jaccard Similarity

```python
def _jaccard(set_a, set_b):
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)
```

- J = 1.0 → identical residue sets
- J = 0.0 → no shared residues
- J ≥ 0.60 → considered overlapping (default threshold)

#### Source Attribution

When a lower-priority motif is discarded, the surviving ref instance's
metadata is updated:

```python
def _annotate_overlap(ref_inst, discarded_inst):
    also = ref_inst.metadata.get('_also_found_in', [])
    d_label = discarded_inst.metadata.get('_source_label', '')
    if d_label and d_label not in also:
        also.append(d_label)
    # Propagate from discarded's own _also_found_in (3+ source cascade)
    for label in discarded_inst.metadata.get('_also_found_in', []):
        if label and label not in also:
            also.append(label)
    ref_inst.metadata['_also_found_in'] = also
```

---

## 9. Key Data Structures

### MotifInstance

Each motif instance carries:

```python
class MotifInstance:
    instance_id: str          # e.g., "HL_75660.8"
    motif_id: str             # e.g., "GNRA"
    pdb_id: str               # e.g., "1S72"
    residues: List[ResidueSpec]
    annotation: str           # Display name
    metadata: dict            # Extra info including:
        # _source_id: int       (e.g., 8)
        # _source_label: str    (e.g., "NoBIAS")
        # _also_found_in: list  (e.g., ["RNAMotifScanX (RMSX)"])
        # enriched_from: str    (e.g., "HL", if enriched)
        # motif_group: str      (e.g., "HL_75660.8")
```

### ResidueSpec

```python
class ResidueSpec:
    chain: str                # e.g., "A"
    residue_number: int       # e.g., 100
    residue_name: str         # e.g., "G" (optional)
    insertion_code: str       # e.g., "" (optional)
```

### SOURCE_ID_MAP

Centralized registry in `rsmviewer/database/config.py`:

| ID | Name | Type | Tool |
|----|------|------|------|
| 1 | RNA 3D Atlas | local | — |
| 2 | Rfam | local | — |
| 3 | BGSU RNA 3D Hub | web | — |
| 4 | Rfam API | web | — |
| 5 | FR3D Annotations | user | fr3d |
| 6 | RNAMotifScan (RMS) | user | rms |
| 7 | RNAMotifScanX (RMSX) | user | rmsx |
| 8 | NoBIAS | user | nobias |

### MOTIF_CANONICAL_MAP

Centralized in `rsmviewer/database/cascade_merger.py`:

| Variant | Canonical |
|---------|-----------|
| `KTURN`, `K_TURN`, `KINK-TURN`, `KINK_TURN`, `KINKTURN` | `K-TURN` |
| `CLOOP`, `C_LOOP` | `C-LOOP` |
| `SARCIN`, `SARCINRICIN`, `SARCIN_RICIN` | `SARCIN-RICIN` |
| `REVERSE-KTURN`, `REVERSE_KTURN`, `REVERSEKTURN`, `REVERSE_K_TURN`, `REVERSE_KINK_TURN`, `REVERSE-KINK-TURN` | `REVERSE-K-TURN` |
| `ELOOP`, `E_LOOP` | `E-LOOP` |
| `TLOOP`, `T_LOOP` | `T-LOOP` |

Names not in this map pass through unchanged (uppercased).

### Naming Conventions by Source

| Source | Example Names |
|--------|--------------|
| Atlas (1) | `HL`, `IL`, `J3`, `J4`, `J5`, `J6`, `J7` (generic, enriched) |
| Rfam (2) | `GNRA`, `T_LOOP`, `SARCIN_RICIN_1`, `SARCIN_RICIN_2`, `K_TURN_1`, `K_TURN_2`, `C_LOOP`, `TWIST_UP` |
| BGSU (3) | Mix of generic (`HL`, `IL`) and specific (`GNRA`, `K-TURN`) |
| Rfam API (4) | Similar to Rfam |
| FR3D (5) | Generic (`HL`, `IL`) (enriched) |
| RMS (6) | `KINK-TURN`, `C-LOOP`, `SARCIN-RICIN`, `REVERSE-KINK-TURN` |
| RMSX (7) | `K-TURN`, `SARCIN-RICIN`, `C-LOOP`, `REVERSE-K-TURN`, `E-LOOP` |
| NoBIAS (8) | Same conventions as RMSX |

### Per-Source Custom Paths

```python
self.user_data_paths: Dict[int, str] = {}
# Example state:
# {
#     7: "/home/user/rmsx_output",
#     8: "/home/user/nobias_output",
# }
```

Each source stores its path independently. `rmv_db 7 /path/to/data` sets
`user_data_paths[7]` without affecting any other source.

---

## 10. FAQ / Troubleshooting

### Q: What Jaccard threshold is used? Can I change it?

**0.60** by default (60% residue overlap). You can override it on the command
line:

```
rmv_db 8 7, jaccard_threshold=0.80   # Use 80% threshold for this combine
rmv_db 1 3, jaccard_threshold=50     # Use 50% threshold
rmv_db 2 7                            # Use default 60%
```

Supported formats: `0.80` (decimal), `80` (percentage), `80%` (with %).

The threshold **persists** until you change it again. Higher values are
**stricter** (fewer merges); lower values are **looser** (more merges).

### Q: Does the order of source IDs matter?

**Yes, critically.** The first source ID has highest priority. If both sources
have a motif covering the same residues (in the same category), the first
source's version is kept.

### Q: Are different motif types ever compared against each other?

**No.** Comparison is **strictly within the same base category**. A C-LOOP
will never be compared against a K-TURN, even if they cover the exact same
residues. Both are preserved because they represent different biological
structures.

### Q: What about SARCIN_RICIN_1 vs SARCIN-RICIN?

They are the **same base category** (`SARCIN-RICIN`). The trailing `_1` is
stripped, then `SARCIN_RICIN` is harmonized to `SARCIN-RICIN`. So yes, they
are compared against each other during the merge.

### Q: What about K_TURN_1 vs K-TURN?

Same principle. `K_TURN_1` → strip `_1` → `K_TURN` → harmonize → `K-TURN`.
They share the base category `K-TURN` and are compared.

### Q: Why does the output have both SARCIN_RICIN_1 and SARCIN-RICIN as keys?

Because original keys are preserved. If Rfam's `SARCIN_RICIN_1` (2 instances)
is the priority source and 3 RMSX `SARCIN-RICIN` instances survive the merge,
the output will have both keys. They were compared (same base category) but
surviving entries keep their own source's naming.

### Q: Does J3, J4, J5 etc. get its suffix stripped?

**No.** The regex requires a separator (`_` or `-`) before the digits:
`[-_]\d+$`. `J3` has no separator, so it stays `J3`. Similarly `J7` stays
`J7`. They are treated as distinct categories.

### Q: Can I combine local + web + user sources?

Yes. For example: `rmv_db 1 3 8` combines Atlas (local) + BGSU (web) +
NoBIAS (user).

### Q: What if a source returns no data?

That source is simply skipped. The pipeline continues with whatever sources
returned results. If **no** source returns data, an error is printed.

### Q: How do I see which source each motif came from?

Use `rmv_show TYPE N` — the instance table includes a `SOURCE` column.
Or use `rmv_summary` which prints a source attribution report at the end.

### Q: I see duplicates in the table — why?

Check if the duplicates are from the **same** source or **different** sources:
- Same source: the within-source dedup should catch these (exact residue match).
  If they have slightly different residue sets, they are considered distinct.
- Different sources: the cascade merge only removes overlaps with Jaccard ≥ 0.60.
  Two motifs with less than 60% residue overlap are considered independent.

### Q: Does enrichment modify residues?

No. Enrichment only changes the **name** (e.g., `HL` → `GNRA`). The residue
coordinates are preserved exactly as they came from the original source.

### Q: Do P-value settings carry over into combine mode?

Yes. Configure P-values **before** entering combine mode:

```
rmv_db 6 K-TURN 0.05       # Set custom P-value for RMS
rmv_db 8 on                # Enable NoBIAS filtering
rmv_db 8 7 6               # Now combine — settings are applied per-source
```

---

*RSMViewer — CBB LAB @Rakib Hasan Rahad*
