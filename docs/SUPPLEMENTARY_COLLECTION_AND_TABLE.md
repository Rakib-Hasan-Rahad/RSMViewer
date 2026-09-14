# Supplementary Documentation: Collection and Consolidated Annotation Tables

This document explains how RSMViewer collects annotations, builds its
structure-local table, preserves source labels and hierarchy, and evaluates
SQL-like queries.

It also explains the difference between the number printed during source
loading and the number of rows printed by `rmv_list`.

## 1. The two counts have different meanings

A typical Atlas load may report:

```text
Found 223 motifs in 1S72
Loaded 223 motif(s) in 44 families from 1S72
```

Later, a selected or combined group may contain fewer physical rows because
overlapping annotations are consolidated within that explicit scope. This is
expected and does not mean that raw annotations were lost.

The counts mean:

| Count | Meaning |
| --- | --- |
| `223 motifs` | Raw motif instances returned by the selected provider and accepted by the loader. |
| `44 families` | Distinct display-family keys in the loaded summary. |
| Source-table rows | Raw source records retained by `rmv_db`. |
| Selected/combined rows | Distinct physical residue-set clusters after scoped consolidation. |
| `L1`, `L2`, `L3` | Source hierarchy labels stored on a selected or combined row. |

Several raw annotations can describe the same physical residue fragment. They
remain separate after `rmv_db`; they occupy one consolidated row only when the
user selects a family or explicitly combines groups. Their source labels are
retained in `source_annotations` and `source_hierarchy`.

For example, these two raw annotations may become one row:

```text
Atlas annotation A: residues {0:248, 0:249, 0:260, 0:261, 0:262, 0:263, 0:264, 0:265}
Atlas annotation B: residues {0:248, 0:249, 0:260, 0:261, 0:262, 0:263, 0:264, 0:265}
```

The table keeps one motif ID and preserves the labels/hierarchy associated with
that residue set. `rmv_list` is a physical-fragment view, not a raw-provider
record dump.

## 2. Collection pipeline

For each fetched structure, RSMViewer follows this sequence:

```text
Structure coordinates
    -> provider retrieval
    -> provider conversion to MotifInstance objects
    -> source-label normalization for matching only
    -> residue normalization into source-specific raw records
    -> separate source tables (no overlap merging at rmv_db)
    -> family selection or explicit group combination
    -> residue-based consolidation within the requested scope
    -> source labels and hierarchy attached to rows
    -> stable motif IDs and query groups
    -> PyMOL operations
```

The source-specific labels are not discarded during normalization. A canonical
family is used to make matching reliable, while the original source label is
retained for display and provenance.

### 2.1 Structure loading

```text
rmv_fetch 1S72
```

`rmv_fetch` loads coordinates into PyMOL and records the structure in the
session. It does not load motif annotations. Each structure receives its own
annotation table. Loading another structure does not overwrite the first
structure's table.

Multiple structures can be fetched together:

```text
rmv_fetch 1S72, 4V88
```

Internally, the tables are indexed by structure ID:

```text
annotation_tables["1S72"]
annotation_tables["4V88"]
```

### 2.2 Source loading

```text
rmv_db RNA3DMotifAtlas
```

The selected provider returns motif instances. For multiple sources:

```text
rmv_db RNA3DMotifAtlas, Rfam
rmv_db RNA3DMotifAtlas, RNAMotifScanX
rmv_db FR3D, RNAMotifScanX
```

The public source names are:

```text
RNA3DMotifAtlas
Rfam
FR3D
RNAMotifScanX
```

Atlas and Rfam are retrieved through their APIs. FR3D and RNAMotifScanX use the
configured external pipelines or preannotated data.

### 2.3 Provider records

Every provider result is converted to a `MotifInstance`. Conceptually, a
record contains:

```text
source
structure_id
motif_type
instance_id
residues
annotation
metadata
```

The important field for consolidation is `residues`. A motif may contain
residues from multiple chains and non-contiguous residue ranges.

## 3. Residue identity and normalization

A residue is represented internally by:

```text
(chain, residue_number, insertion_code, model)
```

For example:

```text
("0", 248, "", 1)
("0", 249, "", 1)
("9", 75, "", 1)
```

The table removes duplicate residue entries and sorts them deterministically.
Chain IDs are preserved as strings, so numeric-looking chains such as `0` and
`9` remain distinct from alphabetic chains such as `A`.

A motif may therefore be displayed compactly as:

```text
0:248-249,260-265
```

Single residues and ranges are both valid:

```text
0:125,128-129
0:2136,2237-2239
```

The parser and selection builder support both forms.

## 4. Motif family names and aliases

RSMViewer uses one centralized family-definition table in
`rsmviewer/database/motif_aliases.py`.

It distinguishes two concepts:

### Original source label

This is retained exactly as returned by the source whenever possible. Examples:

```text
Kink-turn
Sarcin-Ricin
Hairpin Loop (HL)
Internal Loop (IL)
```

### Canonical query family

This is used only for matching equivalent names across sources:

```text
Kink-turn       -> K-TURN
KINK-TURN       -> K-TURN
Kink_Turn       -> K-TURN
KT              -> K-TURN

Sarcin-Ricin    -> SARCIN-RICIN
Sarcin-ricin    -> SARCIN-RICIN
SR              -> SARCIN-RICIN
```

Thus this query works regardless of the source's spelling:

```text
rmv_select Kink-Turn, 1S72, RNA3DMotifAtlas, as group_KT
```

The same query family can match `Kink-turn`, `KINK-TURN`, `K-TURN`, and `KT`.

`K-TURN` and `REVERSE-K-TURN` are intentionally different families. The
matching logic checks the more specific reverse form before the plain K-turn
form when reading FR3D free-text names.

## 5. Single structure, single source

Example:

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
```

The provider may return 223 raw instances. The table then processes each
instance in deterministic order:

```text
for each motif instance:
    normalize its residue set
    find a matching existing row for structure 1S72
    if no row matches:
        create a new row
    otherwise:
        attach the source label and hierarchy to the existing row
```

The result is a table such as:

```text
MOTIF_ID       RESIDUES                 RNA3DMotifAtlas.L1  L2
1S72_00001     0:21-27,516-522         Triple sheared      -
1S72_00010     0:77-81,93-100         Kink-turn           -
1S72_00016     0:158-164,171-178      Sarcin-Ricin        -
1S72_00024     0:248-249,260-265      Internal Loop (IL)   Kink-turn
```

The row `1S72_00024` is one physical residue row with two hierarchy levels,
not two independent stable motif IDs.

## 6. Row matching rules

For every existing row belonging to the same structure, RSMViewer calculates:

### Jaccard similarity

$$
J(A,B)=\frac{|A\cap B|}{|A\cup B|}
$$

Default threshold:

```text
DEFAULT_JACCARD_THRESHOLD = 0.60
```

### Containment / overlap coefficient

$$
C(A,B)=\frac{|A\cap B|}{\min(|A|,|B|)}
$$

Default threshold:

```text
DEFAULT_CONTAINMENT_THRESHOLD = 0.80
```

### Match condition

Two residue sets match when either condition is true:

```text
Jaccard(A, B) >= 0.60
OR
Containment(A, B) >= 0.80
```

Example:

```text
Atlas core:       {1, 2, 3}
Rfam extension:   {1, 2, 3, 4, 5, 6}
```

```text
Jaccard     = 3 / 6 = 0.50
Containment = 3 / 3 = 1.00
```

They are consolidated into one row because containment is at least `0.80`.
This handles sources that describe the same motif at different granularity.

If multiple existing rows qualify, the table chooses the strongest candidate
by `(Jaccard, containment)` and breaks an exact tie deterministically by motif
ID.

Rows from different structures never match, even if their residue numbers are
identical:

```text
1S72 chain 0 residue 100 != 4V88 chain 0 residue 100
```

## 7. Source labels and hierarchy storage

Each row contains dictionaries keyed by canonical public source name:

```python
row.source_annotations = {
    "RNA3DMotifAtlas": ("Sarcin-Ricin",)
}

row.source_hierarchy = {
    "RNA3DMotifAtlas": ("Internal Loop (IL)", "Kink-turn")
}
```

Labels are accumulated rather than overwritten. If two source records match the
same row, each source receives its own entry:

```python
row.source_annotations = {
    "RNA3DMotifAtlas": ("Sarcin-Ricin",),
    "Rfam": ("sarcin-ricin-1",),
}
```

If one source provides a generic parent and a specific child, those become
hierarchy levels:

```text
RNA3DMotifAtlas.L1 = Internal Loop (IL)
RNA3DMotifAtlas.L2 = Kink-turn
```

A label with a spelling difference such as `Sarcin-Ricin` versus
`Sarcin-ricin` is normalized for queries, but the source label remains
available for provenance/display.

## 8. Stable motif IDs

After a new physical row is created, rows for that structure are sorted by
normalized residue set and reindexed:

```text
1S72_00001
1S72_00002
1S72_00003
```

The ID identifies the consolidated physical fragment, not a raw provider
record. Adding another annotation to an existing row does not create another
ID.

IDs are structure-local. The same sequence number in another structure is a
different motif:

```text
1S72_00010
4V88_00010
```

## 9. Multiple structures, one source

Example:

```text
rmv_fetch 1S72, 4V88
rmv_db RNA3DMotifAtlas
```

RSMViewer loads each structure independently:

```text
Table 1: 1S72 rows -> 1S72_00001, 1S72_00002, ...
Table 2: 4V88 rows -> 4V88_00001, 4V88_00002, ...
```

A query can select both:

```text
rmv_select KT, 1S72 and 4V88, RNA3DMotifAtlas, as group_KT
```

The resulting group contains IDs from both structures. They do not overwrite
each other and can be converted into separate PyMOL objects:

```text
rmv_create_object group_KT
rmv_super group_KT
```

## 10. One structure, multiple sources

Example:

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas, Rfam
```

The sources are retained separately in the structure-local annotation table.
`rmv_db` does not perform containment filtering, Jaccard merging, or
cross-source consolidation. Source-specific family tables report raw counts
independently. Consolidation is deferred to `rmv_select` and
`rmv_combine_groups`:

```text
Source: RNA3DMotifAtlas
SELECTABLE NAME       ANNOTATION NAME       COUNT
Sarcin-Ricin          Sarcin-Ricin          8

Source: Rfam
SELECTABLE NAME       ANNOTATION NAME       COUNT
Sarcin-Ricin          sarcin-ricin-1 /      5
                       sarcin-ricin-2
```

A shared-source query requires that **both sources label the same row as the
queried motif**:

```text
rmv_select SR, 1S72, RNA3DMotifAtlas and Rfam, as group_shared_SR
```

The source predicate is evaluated per source: a source name is true for a row
only when that source labels the row as the queried motif. A row where Atlas
calls it Sarcin-Ricin but Rfam calls it something else does not satisfy
`RNA3DMotifAtlas and Rfam`. This is what makes benchmarking (Application 5)
exact: `not RNA3DMotifAtlas and RNAMotifScanX` selects rows RMSX calls SR that
Atlas does not, even when Atlas labels that row a different family.

## 11. Multiple structures, multiple sources

Example:

```text
rmv_fetch 1S72, 4V88
rmv_db RNA3DMotifAtlas, Rfam
rmv_select SR, 1S72 and 4V88, RNA3DMotifAtlas and Rfam, as group_shared
```

The query is evaluated as:

```text
for each table in {1S72, 4V88}:
    keep rows whose structure is 1S72 or 4V88
    keep rows where both RNA3DMotifAtlas and Rfam label the row as SR
```

The resulting group stores stable IDs only. Viewing, object creation, and
superimposition resolve those IDs back to their owning structure table.

## 12. SQL-like query semantics

RSMViewer uses a SQL-like command grammar; `rmv_select` is not raw SQL and users
do not write SQL statements directly.

```text
rmv_select <motif>, <structures>, <sources>, as <group>
```

### Motif predicate

Canonical family matching applies across source spellings:

```text
rmv_select K-TURN, all, RNA3DMotifAtlas, as group_KT
rmv_select Kink-turn, all, RNA3DMotifAtlas, as group_KT
rmv_select KT, all, RNA3DMotifAtlas, as group_KT
```

### Structure predicate

```text
rmv_select SR, 1S72, RNA3DMotifAtlas, as group_one
rmv_select SR, 1S72 and 4V88, RNA3DMotifAtlas, as group_two
rmv_select SR, all, RNA3DMotifAtlas, as group_all
```

### Source predicate

Operator precedence is:

```text
not > and > or
```

Examples:

```text
# Shared annotations on the same consolidated row
rmv_select SR, all, RNA3DMotifAtlas and Rfam, as group_shared

# Atlas-only rows
rmv_select SR, all, RNA3DMotifAtlas and not Rfam, as group_atlas_only

# RMSX rows without Atlas support
rmv_select SR, all, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP
```

The source predicate checks which sources label the row as the queried motif.
It does not query raw API responses after the table has been built.

## 13. FR3D query-name matching

FR3D query files may use family names in their filename or JSON `name` field:

```json
"name": "geometric_3_sarcin3geometric"
```

or:

```text
geometric_5_sarcin_ricin.json
```

RSMViewer uses a guarded free-text family scan for FR3D names:

```text
sarcin / sarcin-ricin       -> SARCIN-RICIN
kink-turn / kturn           -> K-TURN
reverse-kturn               -> REVERSE-K-TURN
```

The specific reverse forms are checked before the general K-turn forms. Thus:

```text
KT query + reverse_kturn FR3D name = no match
KT query + kink_turn FR3D name     = match
SR query + sarcin_ricin FR3D name  = match
```

Other sources use exact canonical label matching rather than FR3D free-text
inference.

## 14. SQLite hierarchy cache

The file displayed by `rmv_list` is:

```text
rsmviewer/database/motif_hierarchy_cache.sqlite3
```

It stores display-oriented hierarchy and alias information. The authoritative
rows used by `rmv_list` and `rmv_select` are held in the in-memory
`ConsolidatedAnnotationTable` for each structure during the PyMOL session.

The SQLite cache records fields equivalent to:

```text
pdb_id
residue_key
custom_id
source_key
source_label
motif_label
hierarchy_level
rank
cached_at
```

The cache is not a second independent annotation table. It supports hierarchy
display, stable lookup, and saved aliases. `rmv_reset` clears the in-memory
state and resets the cache connection/data.

## 15. Why nested rows appear in `rmv_list`

A row such as:

```text
1S72_00024  0:248-249,260-265  Internal Loop (IL)  Kink-turn  -
```

means:

- one physical residue set;
- source hierarchy level 1: `Internal Loop (IL)`;
- source hierarchy level 2: `Kink-turn`.

It does not represent a missing record. It represents preserved classification
information from the same source.

Likewise, a row with several source labels represents one residue cluster with
multiple annotations, not a lost or overwritten source record.

## 16. Interpreting the supplied 1S72 output

The supplied output reports:

```text
223 raw motifs
44 families
```

The subsequent listing contains fewer stable rows because consolidation merged
raw records with matching or nested residue sets. The hierarchy examples in the
listing are valid and show that levels are being preserved:

```text
1S72_00024  Internal Loop (IL)  Kink-turn
1S72_00030  3-way Junction       4-way Junction  Hairpin Loop
1S72_00102  4-way Junction       Major groove platform
1S72_00126  Intercalated tWH      Internal Loop (IL)
1S72_00161  3-way Junction       6-way Junction
1S72_00182  Hairpin Loop (HL)     LSU A loop
```

The apparent discrepancy is therefore a raw-instance-versus-physical-row
comparison. To audit a session programmatically, compare:

```python
raw_count = sum(
    len(instances)
    for instances in provider_motifs.values()
)
row_count = len(gui.annotation_tables["1S72"].rows)
```

Then inspect each row's `source_annotations`, `source_hierarchy`, and
`residue_set` rather than comparing only printed family counts.

## 17. Known spelling corrections

The Atlas output contains the label `Ribsomal LSU H95` in the supplied session.
If that spelling comes directly from the API, RSMViewer preserves it as source
provenance. The canonical family table accepts both:

```text
Ribsomal LSU H95
Ribosomal LSU H95
```

This allows queries to remain robust without silently rewriting the provider's
original label.

## 18. Practical audit commands

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_list
rmv_select KT, 1S72, RNA3DMotifAtlas, as group_KT
rmv_list group_KT
rmv_create_object group_KT
rmv_super group_KT
```

For source comparison:

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas, Rfam
rmv_select SR, 1S72, RNA3DMotifAtlas and Rfam, as group_shared
rmv_list group_shared
```

The table, not the raw provider count, is the basis for all subsequent
selection, visualization, object creation, superimposition, and export.

## 19. Concrete audit example: why the numbers differ

The following is a real 1S72 combined-source check using:

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas, Rfam
```

The consolidated family table reported these selected values:

```text
Family                  Rows containing the family
Hairpin Loop (HL)       32
GNRA                    24
Bulged                  22
Internal Loop (IL)      18
Kink-turn                8
Sarcin-Ricin             8
C-loop                   6
```

The important audit is:

```text
rmv_select SR, 1S72, RNA3DMotifAtlas or Rfam, as probe_sr
```

Result:

```text
rmv_select SR ...       8 motif IDs
Family table SR         8 rows
```

The same check for Kink-turn also returned eight motif IDs and eight table
rows. This confirms that the family table is counting consolidated rows using
the same labels and canonical aliases as `rmv_select`.

These numbers must not be added together as though they were disjoint. For
example, one physical row can contain both `Sarcin-Ricin` and another family
label. That row contributes once to the Sarcin-Ricin count and once to the
other family count. Therefore:

```text
sum of family counts >= number of consolidated rows
```

is normal. It does not imply duplicate motif IDs. The row count answers
“How many distinct physical residue fragments are present?”, while a family
count answers “How many physical fragments can be selected by this family
query?”

The same distinction explains the earlier single-source example:

```text
223 raw provider instances
210 consolidated physical rows
```

The difference of 13 means that at least 13 raw records were absorbed into
existing residue-set rows. They remain visible through source labels and
hierarchy, but they do not receive additional stable IDs. The correct audit is
therefore to compare raw instances with `len(annotation_tables["1S72"].rows)`
and to compare a family-table value with the corresponding `rmv_select`
result, rather than comparing a raw total with a family total.

## 20. Suspicious hierarchy worth verifying: `1FFK_00043`

The following row was observed in the real combined `1FFK` Atlas+Rfam load:

```text
MOTIF_ID     RESIDUES ON CHAIN 0
1FFK_00043   0:567-572,585-591
```

The source labels attached to this one consolidated physical row were:

```text
RNA3DMotifAtlas:
        annotation = Sarcin-Ricin
        hierarchy  = Sarcin-Ricin

Rfam:
        annotations = C-loop, sarcin-ricin-1
        hierarchy   = L1 C-loop -> L2 sarcin-ricin-1
```

In other words, both providers point to the same chain-0 residue set, but they
classify it differently. Atlas gives the row one broad `Sarcin-Ricin` label.
Rfam contributes two labels on the same row: `C-loop` and
`sarcin-ricin-1`. Even within Rfam, the family name is therefore not identical
at every level: the parent is `C-loop`, while the child is the indexed name
`sarcin-ricin-1`.

This is worth verifying against the original Rfam annotation. It may be a
legitimate nested Rfam classification in which a C-loop is also assigned to a
specific sarcin-ricin family, but the current table alone cannot establish
whether `C-loop -> sarcin-ricin-1` is the intended biological hierarchy or an
annotation/name-mapping artifact.

The row is retained as data, not corrected automatically. That is deliberate:
changing or removing either Rfam label would alter source provenance and could
change valid queries. Until the source annotation is confirmed, treat this as
one physical motif with multiple source classifications, not as two separate
motif instances. For example, an inspection should show:

```text
rmv_list 1FFK_00043

1FFK_00043  0:567-572,585-591
    RNA3DMotifAtlas: Sarcin-Ricin
    Rfam: C-loop; sarcin-ricin-1
```

## 21. Worked example: every number in a combined load

This section reproduces a full `RNA3DMotifAtlas, Rfam` combined load for `1S72`
and `1FFK` and explains, step by step, how each printed number is produced. All
values below are the actual current outputs.

```text
rmv_fetch 1S72, 1FFK
rmv_db RNA3DMotifAtlas, Rfam
```

### 21.1 Step 1 — raw annotations fetched per source

Each provider is queried independently. The counts printed on the
`... raw annotations in N categories` lines are the raw provider totals before
any consolidation.

| Structure | RNA3DMotifAtlas raw | Rfam raw |
| --- | --- | --- |
| 1S72 | 223 (44 categories) | 35 (10 categories) |
| 1FFK | 214 (44 categories) | 10 (6 categories) |

These are counts of raw provider instances, not table rows. A single physical
motif can be described by several raw instances.

### 21.2 Step 2 — deferred consolidation

The raw source records remain separate after `rmv_db`. Two residue sets are
merged only within the scope requested by `rmv_select`, or when the user
explicitly calls `rmv_combine_groups`, using

```text
Jaccard(A, B) >= 0.60   OR   containment(A, B) >= 0.80
```

The selected or combined group contains physical residue rows. Its source
columns retain the original database labels, including different family
assignments for the same merged region. The original source records and input
groups remain unchanged.

### 21.3 Step 3 — family counts in the source tables

The `COUNT` column is source-specific. It counts, for each canonical family,
how many raw rows from that source carry at least one label that maps to that
family. The selectable name is canonicalized for matching, while the annotation
column retains the source's original wording. Selection may subsequently merge
rows within its requested family and source scope, so a source-table count is
not a combined-group row count.

A row is counted once per family it belongs to. Because one row can carry more
than one family label, a row can be counted under several families.

Worked 1S72 examples:

```text
Sarcin-Ricin = 8    rows whose Atlas 'Sarcin-Ricin' or Rfam
                    'sarcin-ricin-1'/'sarcin-ricin-2' label is present
Kink-turn    = 8
C-loop       = 6
GNRA         = 24
```

Cross-check against selection (union semantics):

```text
rmv_select SR,   1S72, RNA3DMotifAtlas or Rfam   -> 8   (table shows 8)  ✓
rmv_select KT,   1S72, RNA3DMotifAtlas or Rfam   -> 8   (table shows 8)  ✓
rmv_select CL,   1S72, RNA3DMotifAtlas or Rfam   -> 6   (table shows 6)  ✓
rmv_select GNRA, 1S72, RNA3DMotifAtlas or Rfam   -> 24  (table shows 24) ✓
```

The `and` form is a stricter, row-level predicate and returns fewer rows because
it requires both sources on the same row:

```text
rmv_select SR, 1S72, RNA3DMotifAtlas and Rfam    -> 5   (subset of the 8)
```

### 21.4 Step 4 — the source-specific `Total` line

`Total` at the bottom of each source table is that source's raw annotation
total. It is intentionally not a cross-source total and not a count of rows in
a selected or combined group. A family count can be smaller or larger than a
physical-row count because one raw row may carry multiple hierarchy labels.

### 21.5 Display-name normalization

The `SELECTABLE NAME` column shows one canonical copyable name per family. Three
names are normalized for selection, because their raw source spelling is
inconsistent while `rmv_select` treats them as one family. The `ANNOTATION NAME`
column retains the original source wording:

```text
Ribsomal LSU H95     -> Ribosomal LSU H95   (source typo)
right_angle-3        -> Right-angle         (Rfam index suffix)
twist_up             -> Twist-up            (Rfam short name)
```

The stored `source_annotations` and displayed annotation values keep the
original spelling for provenance/display; only the selectable family label is
normalized. Queries such as
`rmv_select Right-angle, ...` and `rmv_select twist_up, ...` both resolve to the
same canonical family.

### 21.6 Reading the numbers together

For each PDB, the load log reports the raw total for each source. The following
source table reports the same source-specific total and breaks it down by
selectable family. A later `rmv_select` or `rmv_combine_groups` result is a
separate, scope-specific physical-fragment view and should not be added to the
raw source totals.


2. "43 families" — that's correct, not a bug
It's genuinely what BGSU Atlas annotates for 1S72 (the large 23S/5S rRNA subunit). The table lists 43 distinct canonical families totaling 223 motif instances. The counter groups every source spelling into its canonical family:

HAIRPIN LOOP (HL) + OTHER HL → one family (Hairpin-Loop)
RIBSOMAL LSU H95 (source typo) → Ribosomal-Lsu-H95