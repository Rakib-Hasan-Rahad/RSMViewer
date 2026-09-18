# Supplementary Information

**RSMViewer: A PyMOL plugin for RNA structural motif visualization**

Rakib Hasan Rahad, Smriti Pranjal, Nabila Shahnaz Khan, Shaojie Zhang, and
Cuncong Zhong

> This supplement expands on the implementation and algorithmic details of
> RSMViewer that are summarized in the main manuscript
> ([MANUSCRIPT.md](MANUSCRIPT.md)). Every statement below reflects the released
> source code. Cross-references to the main text are marked as *(main text)*.
> End-user command usage is in [../README.md](../README.md) and
> [TUTORIAL.md](TUTORIAL.md); the maintainer reference is in
> [../DEVELOPED.md](../DEVELOPED.md); the collection-and-table rules are in
> [SUPPLEMENTARY_COLLECTION_AND_TABLE.md](SUPPLEMENTARY_COLLECTION_AND_TABLE.md).

---

## Contents

- [S1. Scope](#s1-scope)
- [S2. Residue-set motif representation](#s2-residue-set-motif-representation)
- [S3. Annotation consolidation](#s3-annotation-consolidation)
- [S4. Motif-family name normalization](#s4-motif-family-name-normalization)
- [S5. SQL-like selection grammar and semantics](#s5-sql-like-selection-grammar-and-semantics)
- [S6. Medoid-based superimposition](#s6-medoid-based-superimposition)
- [S7. Annotation source integration](#s7-annotation-source-integration)
- [S8. Caching architecture](#s8-caching-architecture)
- [S9. Worked example: benchmarking (Application 5)](#s9-worked-example-benchmarking-application-5)
- [S10. Visualization, objects, and export](#s10-visualization-objects-and-export)
- [S11. Configuration](#s11-configuration)
- [S12. Reproducibility and testing](#s12-reproducibility-and-testing)
- [S13. Design choices and limitations](#s13-design-choices-and-limitations)
- [S14. Command summary](#s14-command-summary)

---

## S1. Scope

The main manuscript describes RSMViewer at the level of concept and workflow: a
residue-set motif representation, redundancy filtering by the Jaccard index, a
two-dimensional consolidated table with per-source hierarchy columns, a
medoid-based superimposition heuristic, and six example applications. This
supplement provides the additional detail needed to reproduce and extend those
results, including (i) the exact merging criteria used in the implementation,
(ii) the precise per-source semantics of the SQL-like source predicate that make
the benchmarking applications (Applications 5–6) well defined, (iii) the data
layouts and execution modes of the external FR3D and RNAMotifScanX pipelines,
(iv) the multi-tier caching architecture, and (v) a fully worked numerical
benchmark that corresponds to Application 5.

---

## S2. Residue-set motif representation

*(main text: "Motif Representation and Annotation Consolidation")* Each motif
instance is represented internally as a **set of residues**, with no constraint
on the number of RNA strands or molecular chains. A residue is identified by the
tuple

$$r = (\text{chain},\; \text{residue number},\; \text{insertion code},\; \text{model}).$$

A motif instance is therefore a set $R = \{r_1, r_2, \dots\}$. This
representation is agnostic to how the motif was defined — by sequence, secondary
structure, noncanonical base-interaction pattern, or 3D geometry — and it admits
multi-strand and multi-chain motifs without special handling.

Three data structures realize this model in the code base:

| Structure | Role |
| --- | --- |
| `ResidueSpec` (`database/base_provider.py`) | One residue: chain, number, nucleotide, insertion code, model. |
| `MotifInstance` (`database/base_provider.py`) | One provider result: instance ID, source label, structure ID, residue list, annotation text, metadata. |
| `AnnotationRow` (`database/consolidated_table.py`) | One consolidated fragment: stable `motif_id`, structure ID, normalized residue set, per-source labels, per-source hierarchy, provenance. |

Residue sets are normalized (sorted, duplicate-free) before comparison so that
set identity is independent of provider ordering.

---

## S3. Annotation consolidation

### S3.1 Redundancy filtering (merging criterion)

*(main text)* The main text states the Jaccard criterion: two residue sets with
a Jaccard index below 60% are retained as distinct fragments; otherwise they are
merged. The Jaccard index of two residue sets $A$ and $B$ is

$$J(A,B) = \frac{|A \cap B|}{|A \cup B|}.$$

**Implementation refinement.** In addition to the Jaccard criterion, the
implementation applies a *containment* (overlap-coefficient) criterion so that a
tight motif core annotated by one source and an extended annotation of the same
motif by another source are recognized as the **same** fragment even when their
Jaccard index is moderate. The overlap coefficient is

$$C(A,B) = \frac{|A \cap B|}{\min(|A|, |B|)}.$$

Two residue sets are merged into one consolidated row when **either** criterion
holds:

$$J(A,B) \ge \tau_J \quad \textbf{or} \quad C(A,B) \ge \tau_C,$$

with the code-defined constants

$$\tau_J = 0.60 \;\; (\texttt{DEFAULT\_JACCARD\_THRESHOLD}), \qquad
\tau_C = 0.80 \;\; (\texttt{DEFAULT\_CONTAINMENT\_THRESHOLD}),$$

both in `database/consolidated_table.py`. The containment criterion captures
nested annotations that Jaccard alone would score below threshold; the Jaccard
criterion still retains neighboring instances that share only closing canonical
base pairs (as noted in the main text). When several existing rows are eligible,
`_find_matching_row` selects the strongest match and breaks ties by `motif_id`,
which keeps consolidation deterministic. These thresholds are intentionally
code-defined rather than user-configurable so that stable motif IDs remain
reproducible across sessions.

### S3.2 Deterministic stable motif identifiers

After redundancy filtering, each residue set receives a unique internal
identifier of the form `<PDBID>_00001`. Identifiers are assigned deterministically:
`_reindex_structure` sorts a structure's rows by `(structure_id, residue_set)`
and assigns sequential five-digit suffixes. Because the ordering is by residue
set and not by provider load order or cache sequence, the same structure and
source selection always produces the same IDs (e.g. `1S72_00016`). Stable IDs
are the currency of all downstream commands (`rmv_view`, `rmv_create_object`,
`rmv_save`, group snapshots).

### S3.3 The consolidated table and hierarchy columns

*(main text)* Each structure has exactly one `ConsolidatedAnnotationTable`
indexed by residue set. Every row stores, per source, the source's own label(s)
(`source_annotations`) and, when the source uses a hierarchical classification,
each level as a separate column (`source_hierarchy`). For example, RNAMotifScanX
rows may expose two levels, printed as `RNAMotifScanX.L1` and `RNAMotifScanX.L2`
(e.g. `REVERSE-K-TURN` at L1 and `SARCIN-RICIN` at L2), whereas a flat source
such as the RNA 3D Motif Atlas occupies a single column. This preserves the
semantics of every source and avoids collapsing disagreements into a single
consensus label, which is essential for the cross-source comparisons in
Applications 3, 5, and 6.

### S3.4 Multi-source residue-based merge

When two or more sources are combined (e.g. `rmv_db RNA3DMotifAtlas,RNAMotifScanX`),
`database/residue_merger.py` (`ResidueMerger`) performs a name-agnostic,
residue-set-based merge. Overlap is decided **only** by residues/chains (exact
match, subset/superset, or Jaccard ≥ threshold) — never by the motif label — so
that two sources describing the same physical region are placed on one row even
when they disagree on the family name. Sources are processed right-to-left, so
the leftmost (highest-precedence) source in the `rmv_db` list wins ties. Each
source's own label and hierarchy are retained on the merged row and persisted in
the display cache (S8).

---

## S4. Motif-family name normalization

Different sources spell the same family differently (`SR`, `Sarcin-Ricin`,
`sarcin-ricin-1`, `sarcin_ricin`). `database/motif_aliases.py` is the single
source of truth for family-name handling:

- `canonical_motif(value)` maps abbreviations, separators, and family-index
  suffixes to one canonical family name; unknown families fall back to a
  normalized self-matching form.
- `labels_match_motif(query, labels, free_text_labels)` returns whether any
  source label denotes the queried family. Matching is exact on the canonical
  family, with a guarded substring fallback for compound/novel labels
  (e.g. `"sarcin-ricin core region"`). Distinct families such as `K-TURN` and
  `REVERSE-K-TURN` never cross-match. FR3D query names, which embed the family
  loosely in a filename or JSON `name`, are matched by keyword through
  `free_text_labels`.

Because the selection grammar and provider ingestion both delegate here, adding
one alias updates the whole system consistently.

---

## S5. SQL-like selection grammar and semantics

### S5.1 Grammar

`rmv_select` takes four comma-separated clauses parsed by
`database/query_parser.py`:

```text
rmv_select <motif>, <structures>, <sources>, as <group>
```

producing a `QueryExpression(motif, structures, sources, group, text)`.

### S5.2 Structure predicate

`<structures>` is a single PDB ID, an `and`/`or` list of IDs, or the keyword
`all`. A row passes when its structure matches the expression.

### S5.3 Source predicate — per-source semantics

`<sources>` is a Boolean expression over source names with operator precedence

$$\texttt{not} \;>\; \texttt{and} \;>\; \texttt{or}.$$

The predicate is evaluated **per source and per motif**: for a given row and a
given queried family $m$, a source name $s$ evaluates to **true** only when
source $s$ labels that row as family $m$ (via `labels_match_motif` on that
source's own labels/hierarchy). Formally, let

$$\mathrm{Src}_m(\text{row}) = \{\, s : s \text{ labels the row as } m \,\}.$$

The row is selected iff the Boolean expression evaluates to true on the set
$\mathrm{Src}_m(\text{row})$, and iff that set is non-empty (so a purely negated
predicate cannot select rows unrelated to $m$). Consequently:

- `A and B` selects rows that **both** $A$ and $B$ label as $m$.
- `not A and B` selects rows that $B$ labels as $m$ but $A$ does **not** — even
  when $A$ annotates the same residues under a *different* family.
- `A and not B` selects rows that $A$ labels as $m$ but $B$ does not.

This per-source rule is what makes Applications 5 and 6 well defined: a row is a
true positive only when the reference *and* the tool both call it the queried
family, a false positive when the tool calls it the family but the reference
does not (regardless of what other label the reference assigns to that region),
and a false negative when the reference calls it the family but the tool does
not. A source that annotates the same residue region under a different family is
therefore **not** counted as agreeing on the queried family.

### S5.4 Query groups

`rmv_select ... as <group>` stores a snapshot: the matching stable motif IDs
plus the original query text. Groups are immutable snapshots that drive
`rmv_list`, `rmv_view`, `rmv_create_object`, `rmv_super`/`rmv_align`,
`rmv_combine_groups`, `rmv_set_color`, and `rmv_save`. `rmv_combine_groups` unions the motif
IDs of several groups and records each member's origin group so that per-source
colors set with `rmv_set_color` are preserved after the combine (S10).

---

## S6. Medoid-based superimposition

*(main text: "Large-scale Motif Superimposition")* PyMOL's `super`/`align`
operate on two objects at a time. RSMViewer extends this to $n$ instances with a
medoid heuristic (`alignment.py`):

1. **Pairwise RMSD.** For every unordered pair $(i, j)$, superimpose temporary
   copies and record the RMSD $d_{ij}$. Failed pairs are recorded as skipped.
2. **Medoid selection.** Compute each instance's average RMSD to the others
   (ignoring failed pairs) and choose the medoid $k = \arg\min_i \frac{1}{n-1}
   \sum_{j \ne i} d_{ij}$.
3. **Superimpose onto the medoid.** Move every non-medoid instance onto $k$.
4. **Report and color.** Print a table (medoid, per-instance color and RMSD to
   the medoid, overall average RMSD, skipped pairs) and color the medoid green
   with distinct colors for the rest.

`rmv_super` uses `cmd.super` (sequence-independent) and `rmv_align` uses
`cmd.align` (sequence-dependent). The pairwise stage is $O(n^2)$ in alignment
calls, which is why the main text notes that dozens of instances superimpose
within a few seconds; parent structures are never replaced by motif fragments.

---

## S7. Annotation source integration

RSMViewer exposes four **public** sources — `RNA3DMotifAtlas`, `Rfam`, `FR3D`,
`RNAMotifScanX` — whose names are case-insensitive. Public names are deliberately
separate from internal adapter identity; no public command exposes adapter names
or numeric IDs.

### S7.1 RNA 3D Motif Atlas and Rfam (API + caching)

The RNA 3D Motif Atlas and Rfam are retrieved from their respective public
**APIs**. Responses are cached locally on first retrieval
(`database/cache_manager.py`), so repeated queries for the same structure are
served offline; `rmv_refresh` bypasses the cache and re-fetches. Atlas returns
generic secondary-structure keys that are re-categorized into semantic families
during ingestion, and both sources' labels and hierarchies are stored per row.

### S7.2 FR3D

By default, FR3D is loaded in `cache` mode (`data_mode: "cache"` in
`config/fr3d_config.json`), which serves previously generated FR3D results for
the structure from `output/fr3d_runs/` and does not execute FR3D or make network
requests. When `data_mode` is set to `"run_from_scratch"`, FR3D is executed
through a **user-provided** official fr3d-python checkout placed
under `external/fr3d/fr3d-python-latest/`; RSMViewer never vendors or modifies
that checkout. The runner (`tools/fr3d_search_runner.py`) resolves the checkout
and a Python interpreter that can import `numpy`, `scipy`, `mmcif-pdbx`, and
`fr3d`; stages the loaded structure as a local mmCIF target so FR3D annotates the
exact coordinates in memory; runs FR3D's own default queries; and ingests the
resulting CSV output through `FR3DConverter`. Bundled geometric queries define
their template from an external reference PDB and therefore require
`allow_network: true` in `config/fr3d_config.json`; a query whose reference is
unreachable is skipped with a message while the remaining queries proceed.
Typical usage:

```text
rmv_setup FR3D            # one-time: install FR3D's Python deps + register
rmv_fr3d status           # verify the checkout and interpreter
rmv_fetch 1S72
rmv_db FR3D               # run FR3D's default queries and ingest results
```

### S7.3 RNAMotifScanX

RNAMotifScanX (RMSX) has two modes, selected by `data_mode` in
`config/rmsx_config.json` and driven by `tools/rmsx_runner.py`.

**Preannotated mode (default).** RSMViewer reads precomputed RMSX results
distributed as a preannotated bundle under `external/rmsx_preannotated/`
(available for download from Figshare:
[https://doi.org/10.6084/m9.figshare.33826795](https://doi.org/10.6084/m9.figshare.33826795)).
The bundle is keyed by PDB ID and chain and contains both the RMSX inputs and the
precomputed alignment outputs:

```text
rmsx_work_default/
└── <pdb_id_lowercase>/                e.g. 1s72/
    ├── _prep_main/                     inputs used to build the targets
    │   ├── <PDB>.pdb                   coordinates
    │   ├── <PDB>.pdb.mca               MC-Annotate output
    │   └── <PDB>_<chain>.rmsx.in/.nch  RMSX target inputs
    └── <chain>/                        one folder per scanned chain, e.g. 0/
        ├── sarcin-ricin_consensus.log  alignment OUTPUT (one per family)
        ├── k-turn_consensus.log
        ├── c-loop_consensus.log
        ├── e-loop_consensus.log
        └── reverse-kturn_consensus.log
```

For each requested PDB, RSMViewer collects the matching `*_consensus.log` files
across **all** chains of a family, concatenating them into a single result file
so that no chain's hits overwrite another's, then applies the per-family P-value
threshold and consolidates. Because enumerating members of the large gzip
archive requires decompressing the whole stream, the extracted logs are cached
per PDB (S8) so later loads and later PyMOL sessions are fast. The preferred
input is an *extracted* `rmsx_work_default/` directory; a compressed
`.tar.gz` archive is a fallback that is extracted once on first use.

**From-scratch mode.** With `"data_mode": "run_from_scratch"` and working
binaries under `external/rmsx/bin/` (`scan`, `MC-Annotate`, optional `rnaview`),
RSMViewer executes the standalone RMSX pipeline. Even then, it first reuses
prebuilt `.rmsx.in`/`.nch` inputs for the requested PDB from the preannotated
directory or archive (skipping MC-Annotate) and only regenerates inputs from
scratch when none are found; it then runs the RMSX `scan` step to produce fresh
alignment logs.

**P-value thresholds.** Per-family acceptance thresholds live only under
`pvalue_thresholds` in `config/rmsx_config.json`; command-line overrides are
intentionally unsupported. The same thresholds apply to preannotated and freshly
generated logs. A result may legitimately contain zero accepted motifs when
every reported P-value exceeds its threshold — this is correct filtering, not a
load failure. The RNAMotifScanX software and binaries are **not** distributed
with RSMViewer.

---

## S8. Caching architecture

RSMViewer uses three complementary caches, all cleared by `rmv_reset cache`:

1. **API response cache** (`database/cache_manager.py`) — provider API responses
   stored outside the install directory with provenance and expiry; bypassed by
   `rmv_refresh`.
2. **SQLite display cache** (`database/motif_hierarchy_cache.py`) — a text-keyed
   store of per-source labels/hierarchies keyed by canonical `source_key`
   values, used to render `.L1`/`.L2` columns and saved-query metadata.
3. **Preannotated RMSX extraction cache** (`tools/rmsx_runner.py`) — per-PDB
   extracted `*_consensus.log` files at
   `output/rmsx_results/.preannotated_cache/<pdb_id>/`, stamped with the source
   archive/directory identity (path, modification time, size) so a changed
   source invalidates the cache automatically.

A snapshot of the SQLite display cache and the preannotated extraction cache is
distributed with the repository so a fresh clone works immediately; both are
regenerated on demand. `rmv_reset` requires an explicit subcommand to actually
reset anything. With no argument it only prints details about the two
subcommands below and performs no reset:

- `rmv_reset cache` — clears the three caches above (data and file, including
  `-wal`/`-shm` for the SQLite cache); loaded objects, query groups, and other
  session state are left untouched.
- `rmv_reset session` — deletes all PyMOL objects and resets session state
  (loaded structures, query groups, source selections, motif loader, custom
  colors); caches on disk are left untouched.

---

## S9. Worked example: benchmarking (Application 5)

This section provides a fully worked numerical benchmark corresponding to
Application 5 of the main text, using the RNA 3D Motif Atlas as the ground-truth
reference for RNAMotifScanX on the large ribosomal subunit structure 1S72 and the
SARCIN-RICIN (SR) family. The exact counts depend on data versions and the
configured P-value thresholds; the values below were reproduced on a clean
session with **only** the two benchmarked sources loaded (RSMViewer 2.0.0,
cached RNA 3D Motif Atlas API data and the preannotated RNAMotifScanX bundle,
$\tau_J = 0.60$, $\tau_C = 0.80$). Load only these two sources: adding a third
source (e.g. Rfam) changes the residue-set consolidation and therefore shifts
the TP/FP/FN partition, even though the per-source SR marginals stay the same.

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas,RNAMotifScanX
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas and RNAMotifScanX, as group_TP
rmv_select SARCIN-RICIN, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas and not RNAMotifScanX, as group_FN
```

Applying the per-source semantics of S5.3 partitions the SR-associated rows into
three disjoint groups that together equal the union
`RNA3DMotifAtlas or RNAMotifScanX`:

| Group | Selection | Meaning | Count |
| --- | --- | --- | ---: |
| TP | `RNA3DMotifAtlas and RNAMotifScanX` | both sources label the row SR | 1 |
| FP | `not RNA3DMotifAtlas and RNAMotifScanX` | RMSX labels SR; Atlas does not | 11 |
| FN | `RNA3DMotifAtlas and not RNAMotifScanX` | Atlas labels SR; RMSX does not | 7 |
| Union | `RNA3DMotifAtlas or RNAMotifScanX` | any source labels SR | 19 |

The partition is disjoint and exhaustive ($1 + 11 + 7 = 19$), and the marginals
are internally consistent: Atlas-SR rows $= \mathrm{TP} + \mathrm{FN} = 8$ and
RMSX-SR rows $= \mathrm{TP} + \mathrm{FP} = 12$, giving the union
$8 + 12 - 1 = 19$. Treating the Atlas as ground truth,

$$\text{Precision} = \frac{\mathrm{TP}}{\mathrm{TP}+\mathrm{FP}} = \frac{1}{12} \approx 8.3\%,$$

$$\text{Recall} = \frac{\mathrm{TP}}{\mathrm{TP}+\mathrm{FN}} = \frac{1}{8} = 12.5\%,$$

$$F_1 = \frac{2\,\mathrm{TP}}{2\,\mathrm{TP}+\mathrm{FP}+\mathrm{FN}} = \frac{2}{20} = 10.0\%.$$

**Interpretation and caveat.** Two subtleties, both consequences of the
representation and thresholds described above, matter when interpreting such a
benchmark:

1. A source predicate is per-source and per-family (S5.3). For instance, an
   RMSX SR row is a false positive whether the Atlas assigns that region a
   *different* family (e.g. `3-way Junction`) or no annotation at all; likewise
   an Atlas SR row that RMSX labels `REVERSE-K-TURN` (or leaves unannotated) is
   a false negative — a same-residue annotation under a different family does
   **not** count as agreement.
2. Several FP/FN pairs describe the *same physical region* with slightly
   different residue boundaries reported by the two tools. Whether an Atlas SR
   row and an RMSX SR row collapse into one TP row is decided by the
   **cross-source residue merge** (S3.4: strict subset/superset, or
   Jaccard $\ge \tau_J = 0.60$) — the containment coefficient $\tau_C$ used for
   within-source consolidation (S3.1) is intentionally *not* applied across
   sources. In the 1S72 SR benchmark the seven overlapping FP/FN pairs have
   residue-set Jaccard indices of $0.44$–$0.53$, just below $\tau_J = 0.60$
   (one pair reaches containment $0.82$ but, lacking a subset relationship and
   with Jaccard $< 0.60$, still does not merge). Reported precision and recall
   should therefore be read as a strict, residue-identity lower bound; the
   biological agreement is likely higher, and such near-duplicate pairs are
   exactly the cases Application 6 is designed to inspect visually.

`rmv_list SARCIN-RICIN` lists all 19 SR-associated rows with their per-source
labels, and `rmv_list group_TP` / `group_FP` / `group_FN` list each partition.

---

## S10. Visualization, objects, and export

- **Highlighting.** `rmv_view <ID|group>` recolors motif residues on the parent
  structure without creating objects; it accepts `color=` and `padding=` and can
  take several targets at once. On a single structure, residues shared by
  overlapping rows take the color applied last.
- **Objects.** `rmv_create_object <ID|group>` builds one selectable
  `motif_<id>` object per row, restricted to polymeric nucleic-acid atoms and
  rendered as a cartoon backbone. When a group was produced by `rmv_combine_groups`,
  each object is colored by the source group it came from, so per-source colors
  set with `rmv_set_color` (e.g. red for false positives, blue for known
  instances in Application 6) are preserved; an explicit color set on the
  combined group overrides this and colors all members uniformly.
- **Export.** `rmv_save <group|ID|ALL> cif` writes minimal, coordinates-only
  mmCIF files using original on-disk coordinates where available;
  `rmv_save current [file.png]` saves the current view as a high-resolution PNG;
  `rmv_save ALL [representation]` saves per-motif images. Every save prints its
  output directory (`motif_structures/<pdb_id>/`, `motif_images/<pdb_id>/`, or
  the provided PNG path).

---

## S11. Configuration

Both external pipelines are configured by small, human-readable JSON files under
`config/` (see `config/README.md` for every field and how relative paths are
resolved):

- `config/fr3d_config.json` — FR3D checkout path, Python interpreter, data mode,
  query directory, `allow_network`, and output directory.
- `config/rmsx_config.json` — RMSX `data_mode`, executable paths, the
  preannotated archive/directory, motif families, output directory, and
  `pvalue_thresholds`.

The residue-overlap thresholds ($\tau_J = 0.60$, $\tau_C = 0.80$) are code-defined
in `database/consolidated_table.py`, not configuration options, to keep stable
motif IDs reproducible.

---

## S12. Reproducibility and testing

From the project root:

```bash
python3 -m compileall -q rsmviewer tests
python3 -m unittest discover -s tests -v
RSMVIEWER_ROOT="$PWD" pymol -cq tests/pymol_applications_e2e.py
```

`tests/pymol_applications_e2e.py` executes all six example applications in a real
headless PyMOL session and writes a PASS/FAIL report. Focused unit tests cover
the merge criteria and stable IDs (`test_consolidated_table.py`), family-name
normalization including the `K-TURN`/`REVERSE-K-TURN` separation
(`test_motif_aliases.py`), the selection grammar and operator precedence
(`test_query_parser.py`), and multi-structure fetching
(`pymol_fetch_multi_smoke.py`).

---

## S13. Design choices and limitations

- **Annotation is preserved, not reconciled.** RSMViewer never collapses
  disagreeing sources into a single consensus label; each source keeps its own
  label and hierarchy. This is a deliberate choice so that agreement,
  source-specific annotation, and conflict all remain queryable.
- **Merging is residue-based, not name-based.** Whether two annotations describe
  the same physical fragment is decided purely from residues (S3.1, S3.4). Two
  annotations of the same region under different family names therefore share a
  row and are distinguished by their per-source labels, not split apart.
- **Strict residue-identity benchmarking.** As shown in S9, near-duplicate
  residue sets that fall just below the merge thresholds are counted separately.
  Benchmark precision/recall are best interpreted as strict lower bounds.
- **External software is user-provided.** FR3D and RNAMotifScanX binaries and
  large offline datasets are not distributed with RSMViewer; the repository
  ships only the directory structure, documentation, and a small sample so that
  users can install the software and paste large data locally.

---

## S14. Command summary

| Category | Commands |
| --- | --- |
| Structures & sources | `rmv_fetch`, `rmv_db`, `rmv_source`, `rmv_refresh` |
| Query & results | `rmv_select`, `rmv_list`, `rmv_combine_groups` |
| Visualization & objects | `rmv_view`, `rmv_hide`, `rmv_create_object`, `rmv_bg_color`, `rmv_toggle` |
| Color | `rmv_set_color`, `rmv_color`, `rmv_colors` |
| Analysis & export | `rmv_super`, `rmv_align`, `rmv_save`, `rmv_pair`, `rmv_pair_batch` |
| External pipelines | `rmv_fr3d`, `rmv_rmsx`, `rmv_rmsx_doctor` |
| Diagnostics & session | `rmv_chains`, `rmv_loaded`, `rmv_debug`, `rmv_help`, `rmv_reset`, `rmv_reset cache`, `rmv_reset session` |

For full syntax and examples see [../README.md](../README.md) and
[TUTORIAL.md](TUTORIAL.md); for the in-PyMOL reference, run `rmv_help`.
