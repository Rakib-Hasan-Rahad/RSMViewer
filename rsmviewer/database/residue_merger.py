"""
Residue-based motif merging - order-aware, residue-set redundancy filtering.

Merges motif datasets from multiple sources by residue-set overlap, following
the manuscript's redundancy-filtering rule: two residue sets are treated as the
same structural fragment when their Jaccard index is at least the threshold
(0.60 by default), or when one set is largely contained in the other; otherwise
they are kept as distinct fragments. Sources are processed right to left so the
leftmost source in the user's input takes precedence when a tie must be broken.

Residue-set based, name-agnostic merging:
  Whether two instances are the "same" physical motif is decided **only**
  by residue/chain overlap (exact match, subset/superset, or Jaccard >=
  threshold) - never by motif type name.  An instance labelled ``IL`` by
  one source and ``K-TURN`` by another still merges/dedups correctly if
  their residues overlap enough; no family name is inferred as "more
  specific" than another.  Every source's own label for a matched
  instance is preserved separately in the hierarchy cache
  (``motif_hierarchy_cache.py``) so the display table can show each
  source's raw label side by side.

  Motif type names are still used for exactly one thing: choosing which
  dict key a *kept, non-overlapping* updater instance is filed under in
  the merged dict returned to the legacy per-type loader (e.g. so
  ``rmv_show HL`` still groups same-category instances together). This
  naming step never affects the overlap/merge decision itself.

Algorithm:
  Given ordered sources [Src1, Src2, Src3]:
    1. result = pairwise_merge(Src2, Src3)   # Src2 takes precedence
    2. final  = pairwise_merge(Src1, result)  # Src1 takes precedence

  In each pairwise merge:
    - For each motif in the lower-precedence set (updater):
      - Compare against every ref instance regardless of its reported
        motif label (no category/name gate):
        1. If updater ⊂ ref (strict subset): DISCARD updater (annotate ref).
        2. If ref ⊂ updater (strict superset): REPLACE ref with updater
           (file under ref's key — precedence-source naming convention).
        3. Else if Jaccard >= threshold: DISCARD updater (annotate ref).
        4. Else: KEEP updater (file under the precedence source's key when
           there is exactly one ref key for the base category, otherwise
           under the updater's own key).
    - Final = ref motifs + non-discarded updater motifs

  Subset/superset checks always take precedence over Jaccard-based
  discard so that the larger, more informative residue set is preserved
  regardless of which source ranked higher.

Uses Jaccard similarity: J(A,B) = |A ∩ B| / |A ∪ B|
Threshold: 0.60 by default

Author: CBB Lab
Version: 1.0.0
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Set, Tuple

from .base_provider import MotifInstance, ResidueSpec

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Name harmonization – canonical mapping for motif type variants
# ---------------------------------------------------------------------------
# Different annotation tools use different names for the same structural motif.
# This map normalises every known variant to a single canonical form so that
# the residue merger compares motifs **within the same biological category**.
#
# All keys MUST be uppercase.  Values are the canonical display form.
MOTIF_CANONICAL_MAP: Dict[str, str] = {
    # K-turns
    'KTURN':            'K-TURN',
    'K_TURN':           'K-TURN',
    'KINK-TURN':        'K-TURN',
    'KINK_TURN':        'K-TURN',
    'KINKTURN':         'K-TURN',
    # C-loops
    'CLOOP':            'C-LOOP',
    'C_LOOP':           'C-LOOP',
    # Sarcin-ricin
    'SARCIN':           'SARCIN-RICIN',
    'SARCINRICIN':      'SARCIN-RICIN',
    'SARCIN_RICIN':     'SARCIN-RICIN',
    # Reverse K-turn
    'REVERSE-KTURN':    'REVERSE-K-TURN',
    'REVERSE_KTURN':    'REVERSE-K-TURN',
    'REVERSEKTURN':     'REVERSE-K-TURN',
    'REVERSE_K_TURN':   'REVERSE-K-TURN',
    'REVERSE_KINK_TURN':'REVERSE-K-TURN',
    'REVERSE-KINK-TURN':'REVERSE-K-TURN',
    # E-loops
    'ELOOP':            'E-LOOP',
    'E_LOOP':           'E-LOOP',
    # T-loops
    'TLOOP':            'T-LOOP',
    'T_LOOP':           'T-LOOP',
}


def _harmonize_name(name: str) -> str:
    """Return the canonical motif type for *name*, or UPPER(name) if unknown."""
    key = name.strip().upper()
    return MOTIF_CANONICAL_MAP.get(key, key)


# Regex that matches a trailing numeric suffix: _1, _2, -1, -2, etc.
_TRAILING_NUM_RE = re.compile(r'[-_]\d+$')


def _get_base_category(name: str) -> str:
    """
    Extract the base motif category for comparison grouping.

    Strips a trailing numeric suffix (``_1``, ``-2``, …) then maps the
    remainder through :data:`MOTIF_CANONICAL_MAP`.

    Examples::

        SARCIN_RICIN_1  → SARCIN-RICIN
        SARCIN_RICIN_2  → SARCIN-RICIN
        K_TURN_1        → K-TURN
        K-TURN          → K-TURN
        SARCIN-RICIN    → SARCIN-RICIN
        GNRA            → GNRA
        J3              → J3   (no separator before digit → no strip)
    """
    key = name.strip().upper()
    stripped = _TRAILING_NUM_RE.sub('', key)
    return MOTIF_CANONICAL_MAP.get(stripped, stripped)


def harmonize_motif_dict(
    motifs: Dict[str, List[MotifInstance]],
) -> Dict[str, List[MotifInstance]]:
    """
    Re-key a motif dict so that variant names are merged under their canonical form.

    Example: {'KINK-TURN': [...], 'K-TURN': [...]} → {'K-TURN': [...combined...]}

    Note: this **replaces** original keys.  The residue merger uses
    :func:`_get_base_category` internally and preserves original keys instead.
    """
    harmonized: Dict[str, List[MotifInstance]] = {}
    for mtype, instances in motifs.items():
        canonical = _harmonize_name(mtype)
        if canonical in harmonized:
            harmonized[canonical].extend(instances)
        else:
            harmonized[canonical] = list(instances)
    return harmonized


def _get_residue_set(instance: MotifInstance) -> Set[Tuple[str, int]]:
    """
    Extract a set of (chain, residue_number) tuples from a MotifInstance.
    
    Using (chain, residue_number) pairs ensures chain-aware comparison.
    """
    return {(r.chain, r.residue_number) for r in instance.residues}


def _jaccard(set_a: Set, set_b: Set) -> float:
    """
    Calculate Jaccard similarity between two sets.
    
    J(A, B) = |A intersect B| / |A union B|
    """
    if not set_a or not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


def _annotate_overlap(ref_inst: MotifInstance, discarded_inst: MotifInstance) -> None:
    """
    Annotate a surviving ref instance with source info from a discarded overlapping instance.
    
    Adds/extends the '_also_found_in' metadata list on the ref instance so that
    the display layer can show e.g. "NoBIAS + RMSX" for instances with cross-source support.
    Also propagates any existing '_also_found_in' from the discarded instance (for 3+ source merges).
    """
    if not hasattr(ref_inst, 'metadata') or ref_inst.metadata is None:
        return
    if not hasattr(discarded_inst, 'metadata') or discarded_inst.metadata is None:
        return
    
    also = ref_inst.metadata.get('_also_found_in', [])
    ref_label = ref_inst.metadata.get('_source_label', '')
    
    # Add the discarded instance's own source label
    # Skip if it matches the ref's own primary source label (same-source subset)
    d_label = discarded_inst.metadata.get('_source_label', '')
    if d_label and d_label not in also and d_label != ref_label:
        also.append(d_label)
    
    # Propagate any _also_found_in from the discarded instance (3+ source merge)
    for label in discarded_inst.metadata.get('_also_found_in', []):
        if label and label not in also and label != ref_label:
            also.append(label)
    
    ref_inst.metadata['_also_found_in'] = also


class ResidueMerger:
    """
    Merges motif datasets from multiple sources by residue-set overlap.
    
    Key principles:
    - Merge proceeds right to left
    - Overlap/merge decisions are residue-set based ONLY (exact match,
      subset/superset, or Jaccard >= threshold) 
      gates whether two instances are treated as the same physical motif
    - Motif type names are used only to choose a display key for the
      merged dict; each source's own raw label is preserved separately in
      the hierarchy cache for per-source display
    """
    
    def __init__(self, jaccard_threshold: float = 0.60):
        """
        Initialize the ResidueMerger.
        
        Args:
            jaccard_threshold: Minimum Jaccard similarity to consider
                              two motifs as overlapping (default 0.60)
        """
        self.threshold = jaccard_threshold
    
    def merge_sources(
        self,
        ordered_sources: List[Dict[str, List[MotifInstance]]],
        source_labels: Optional[List[str]] = None,
        pdb_id: Optional[str] = None,
        source_keys: Optional[List[str]] = None,
    ) -> Dict[str, List[MotifInstance]]:
        """
        Merge multiple motif datasets right to left.
        
        Args:
            ordered_sources: List of motif dicts, ordered by precedence (index 0 = highest precedence).
                            Each dict maps motif_type -> [MotifInstance, ...]
            source_labels: Optional labels for logging (e.g., ['Source 1', 'Source 3'])
            pdb_id: PDB ID, used only to persist the source label hierarchy to
                    the SQLite cache (rsmviewer/database/motif_hierarchy_cache.py)
            source_keys: Canonical source names matching source_labels order.
        
        Returns:
            Merged motif dict: motif_type -> [MotifInstance, ...]
        """
        if not ordered_sources:
            return {}
        
        if len(ordered_sources) == 1:
            result = ordered_sources[0]
            if pdb_id and source_keys and source_labels:
                self._persist_hierarchy(result, pdb_id, source_keys, source_labels)
            return result
        
        labels = source_labels or [f"Source_{i}" for i in range(len(ordered_sources))]
        
        logger.info(f"[MERGE] Merging {len(ordered_sources)} sources: {labels}")
        
        # Right-to-left merge: start from the end
        # Start with the last (lowest precedence) source
        result = ordered_sources[-1]
        logger.info(
            f"[MERGE] Starting with {labels[-1]}: "
            f"{self._count_motifs(result)} motifs in {len(result)} categories"
        )
        
        # Merge each source from right to left (second-to-last back to first)
        for i in range(len(ordered_sources) - 2, -1, -1):
            ref = ordered_sources[i]
            updater = result
            
            ref_label = labels[i]
            updater_label = f"intermediate" if i < len(ordered_sources) - 2 else labels[i + 1]
            
            result = self._pairwise_merge(ref, updater, ref_label, updater_label)
            
            logger.info(
                f"[MERGE] After merging {ref_label} (precedence) + {updater_label}: "
                f"{self._count_motifs(result)} motifs in {len(result)} categories"
            )
        
        logger.info(
            f"[MERGE] Final result: "
            f"{self._count_motifs(result)} motifs in {len(result)} categories"
        )
        
        # Final pass: remove subset instances within the same category.
        # The pairwise merge only checks updater-vs-ref.  Two instances from
        # the SAME source (both ref) are never compared during the merge,
        # so a small instance can survive even though it is a complete subset
        # of a larger instance from the same source.
        result = self._deduplicate_subsets(result)
        
        if pdb_id and source_keys and source_labels:
            self._persist_hierarchy(result, pdb_id, source_keys, source_labels)
        
        return result

    def _persist_hierarchy(
        self,
        merged: Dict[str, List[MotifInstance]],
        pdb_id: str,
        source_keys: List[str],
        source_labels: List[str],
    ) -> None:
        """Record each final instance's per-source label chain in the hierarchy cache."""
        from .motif_hierarchy_cache import get_hierarchy_cache, residue_key

        label_to_rank = {label: i for i, label in enumerate(source_labels)}
        label_to_key = {label: key for label, key in zip(source_labels, source_keys)}
        cache = get_hierarchy_cache()

        for motif_label, instances in merged.items():
            for inst in instances:
                rset = _get_residue_set(inst)
                if not rset:
                    continue
                r_key = residue_key(list(rset))
                meta = inst.metadata or {}
                primary_label = meta.get('_source_label', '')
                contributing = [primary_label] + list(meta.get('_also_found_in', []))
                # The surviving merged instance's own generic parent type
                # (e.g. BGSU's IL before it was refined to K-TURN) is only
                # known for its primary/original source - other contributing
                # sources' generic type was discarded when their instance
                # was dropped during the merge, so they get one level only.
                generic = meta.get('_generic_type') or meta.get('loop_type')
                for label in contributing:
                    if not label or label not in label_to_key:
                        continue
                    if label == primary_label and generic and generic.upper() != motif_label.upper():
                        cache.record(
                            pdb_id, r_key, label_to_key[label], label,
                            generic.upper(), label_to_rank.get(label, 99), hierarchy_level=1,
                            jaccard_threshold=self.threshold,
                        )
                        cache.record(
                            pdb_id, r_key, label_to_key[label], label,
                            motif_label, label_to_rank.get(label, 99), hierarchy_level=2,
                            jaccard_threshold=self.threshold,
                        )
                    else:
                        cache.record(
                            pdb_id, r_key, label_to_key[label], label,
                            motif_label, label_to_rank.get(label, 99),
                            jaccard_threshold=self.threshold,
                        )
    
    def _pairwise_merge(
        self,
        ref: Dict[str, List[MotifInstance]],
        updater: Dict[str, List[MotifInstance]],
        ref_label: str = "ref",
        updater_label: str = "updater",
    ) -> Dict[str, List[MotifInstance]]:
        """
        Merge two motif datasets with ref having precedence over updater.

        **Residue/chain-based, name-agnostic**: overlap is checked between
        every updater instance and every ref instance regardless of their
        reported motif label (no category/name gate).  Two instances that
        occupy the same residues are treated as the same physical motif
        even if one source calls it ``IL`` and another calls it ``K-TURN``;
        both labels are preserved via ``_persist_hierarchy`` for the
        hierarchy table, while the merged dict below keeps one display key
        per the existing precedence-naming convention.

        Original motif-type keys from each source are **preserved** in the
        output so that the precedence source's naming convention is respected.
        """
        # --- Start with all ref motifs under their original (uppercased) keys ---
        merged: Dict[str, List[MotifInstance]] = {}
        for mtype, instances in ref.items():
            key = mtype.upper()
            if key in merged:
                merged[key].extend(instances)
            else:
                merged[key] = list(instances)

        # Pre-compute ref residue sets grouped by **base category** (kept
        # only for the precedence-naming resolution of unique updater motifs;
        # the overlap comparison itself uses the flattened `ref_all` below).
        # Each entry stores (original_key, instance, residue_set)
        ref_by_cat: Dict[str, List[Tuple[str, MotifInstance, Set]]] = {}
        for mtype, instances in ref.items():
            cat = _get_base_category(mtype)
            bucket = ref_by_cat.setdefault(cat, [])
            for inst in instances:
                rset = _get_residue_set(inst)
                if rset:
                    bucket.append((mtype.upper(), inst, rset))

        # Flat candidate pool for overlap comparison: every ref instance,
        # regardless of its base category.  This is what makes the merge
        # residue/chain-based instead of name-first: an updater instance is
        # compared against ALL ref instances, not just same-named ones.
        ref_all: List[Tuple[str, MotifInstance, Set]] = [
            entry for bucket in ref_by_cat.values() for entry in bucket
        ]

        # Map: base_category -> set of distinct ref keys that share that
        # category.  Used to decide whether a surviving (unique) updater
        # motif should adopt the ref's key (precedence naming) or keep its own.
        # Rule: if exactly ONE ref key represents the base category, the
        # updater inherits that key so all motifs of the same biological
        # category appear together (e.g. KINK-TURN survivors filed under
        # K-TURN when source 7 takes precedence).  When the ref source uses
        # multiple suffixed keys for the same category (e.g. Rfam
        # SARCIN_RICIN_1 + SARCIN_RICIN_2), the updater keeps its own key
        # to avoid arbitrary bucketing.
        ref_keys_by_cat: Dict[str, List[str]] = {}
        for mtype in ref.keys():
            cat = _get_base_category(mtype)
            keys = ref_keys_by_cat.setdefault(cat, [])
            k_upper = mtype.upper()
            if k_upper not in keys:
                keys.append(k_upper)

        def _resolve_target_key(u_cat: str, u_key: str) -> str:
            """Choose the dict key under which a kept updater motif should
            be filed.  Prefers a single matching ref key for the same
            category (precedence naming convention)."""
            ref_keys = ref_keys_by_cat.get(u_cat, [])
            if len(ref_keys) == 1:
                return ref_keys[0]
            return u_key

        # Check each updater motif against every ref motif (residue/chain based)
        stats = {'kept': 0, 'discarded': 0, 'subsumed': 0, 'replaced': 0}

        for u_type, u_instances in updater.items():
            u_key = u_type.upper()
            u_cat = _get_base_category(u_type)
            # Candidate pool is every ref instance, not just same-named ones:
            # overlap is decided by residues/chain, never by motif label.
            same_cat_refs = ref_all

            for u_inst in u_instances:
                u_rset = _get_residue_set(u_inst)
                if not u_rset:
                    # No residues – keep it anyway, filed under canonical key
                    target_key = _resolve_target_key(u_cat, u_key)
                    merged.setdefault(target_key, []).append(u_inst)
                    stats['kept'] += 1
                    continue

                # --- Phase 1: Asymmetric containment (subset / superset) ---
                # Subset relationships ALWAYS take precedence over Jaccard
                # precedence-based discard: the larger residue set is the
                # more informative annotation regardless of which source
                # has higher precedence.  
                #
                # • u ⊂ r (strict) : updater adds nothing  → discard updater
                # • r ⊂ u (strict) : updater is a superset → replace ref
                #                    with updater (file under ref's key)
                containment_action = None

                for _r_key, r_inst, r_rset in same_cat_refs:
                    if not r_rset or u_rset == r_rset:
                        # Equal sets fall through to Jaccard branch (which
                        # will discard the updater as duplicate, J=1.0).
                        continue
                    if u_rset.issubset(r_rset):
                        containment_action = ('discard_updater', r_inst, _r_key, r_rset)
                        break
                    if r_rset.issubset(u_rset):
                        containment_action = ('replace_ref', r_inst, _r_key, r_rset)
                        break

                if containment_action:
                    action, c_ref_inst, c_ref_key, c_ref_rset = containment_action
                    if action == 'discard_updater':
                        stats['subsumed'] += 1
                        _annotate_overlap(c_ref_inst, u_inst)
                        logger.debug(
                            f"[MERGE] Subsumed {updater_label} motif "
                            f"'{u_key}/{u_inst.instance_id}' - subset of "
                            f"{ref_label} '{c_ref_key}' ({len(u_rset)}/{len(c_ref_rset)} residues)"
                        )
                    else:
                        # Replace smaller ref with larger updater's residues,
                        # but keep ref's KEY (precedence naming convention).
                        stats['replaced'] += 1
                        # Remove ref instance from merged
                        if c_ref_key in merged:
                            merged[c_ref_key] = [
                                inst for inst in merged[c_ref_key]
                                if inst is not c_ref_inst
                            ]
                            if not merged[c_ref_key]:
                                del merged[c_ref_key]
                        # Use updater's residues but file under REF's key
                        # (preserves precedence source naming convention)
                        merged.setdefault(c_ref_key, []).append(u_inst)
                        # Build combined source label (e.g. "NoBIAS + RMSX")
                        upd_label_str = u_inst.metadata.get('_source_label', '')
                        ref_label_str = c_ref_inst.metadata.get('_source_label', '')
                        if ref_label_str:
                            u_inst.metadata['_source_label'] = ref_label_str
                        also = list(u_inst.metadata.get('_also_found_in', []))
                        if upd_label_str and upd_label_str != ref_label_str and upd_label_str not in also:
                            also.append(upd_label_str)
                        for lbl in c_ref_inst.metadata.get('_also_found_in', []):
                            if lbl and lbl != ref_label_str and lbl not in also:
                                also.append(lbl)
                        u_inst.metadata['_also_found_in'] = also
                        # Update ref_all so later comparisons use the larger set
                        for idx, (rk, ri, rs) in enumerate(ref_all):
                            if ri is c_ref_inst:
                                ref_all[idx] = (c_ref_key, u_inst, u_rset)
                                break
                        logger.debug(
                            f"[MERGE] Replaced {ref_label} motif "
                            f"'{c_ref_key}' with larger {updater_label} "
                            f"'{u_key}/{u_inst.instance_id}' "
                            f"({len(c_ref_rset)} → {len(u_rset)} residues)"
                        )
                    continue

                # --- Phase 2: Standard Jaccard overlap ---
                # Only reached when neither set is a strict subset of the
                # other.  If overlap is high enough, prefer the precedence
                # (ref) version and discard the updater.
                is_overlap = False
                best_j = 0.0
                best_ref_inst = None

                for _r_key, r_inst, r_rset in same_cat_refs:
                    j = _jaccard(u_rset, r_rset)
                    if j >= self.threshold:
                        is_overlap = True
                        if j > best_j:
                            best_j = j
                            best_ref_inst = r_inst
                        break

                if is_overlap:
                    stats['discarded'] += 1
                    if best_ref_inst is not None:
                        _annotate_overlap(best_ref_inst, u_inst)
                    logger.debug(
                        f"[MERGE] Discarded {updater_label} motif "
                        f"'{u_key}/{u_inst.instance_id}' - overlaps with "
                        f"{ref_label} instance by residues (Jaccard={best_j:.3f})"
                    )
                    continue

                # --- Truly unique ---
                # File under the precedence (ref) key when the base category
                # has exactly one ref key, otherwise keep updater's key.
                target_key = _resolve_target_key(u_cat, u_key)
                merged.setdefault(target_key, []).append(u_inst)
                stats['kept'] += 1

        logger.info(
            f"[MERGE] Pairwise merge {ref_label} + {updater_label}: "
            f"kept {stats['kept']}, discarded {stats['discarded']} (Jaccard), "
            f"subsumed {stats['subsumed']}, replaced {stats['replaced']}"
        )

        return merged
    
    def _count_motifs(self, motifs: Dict[str, List[MotifInstance]]) -> int:
        """Count total motif instances across all types."""
        return sum(len(instances) for instances in motifs.values())

    # ------------------------------------------------------------------
    # Post-merge subset deduplication
    # ------------------------------------------------------------------
    def _deduplicate_subsets(
        self,
        motifs: Dict[str, List[MotifInstance]],
    ) -> Dict[str, List[MotifInstance]]:
        """
        Remove instances whose residue set is a complete subset of another
        instance, regardless of motif label (residue/chain based, not name-gated).

        This catches subset pairs that the pairwise merge cannot see — e.g.
        two instances from the same source where one has 2 residues and the
        other has 11 residues that fully contain the first two.

        The larger (superset) instance is kept and annotated with the
        smaller instance's source info.  If both sets are identical this
        is a no-op (within-source exact dedup already ran).
        """
        # Flatten every instance across all keys into one comparison pool.
        # Each entry: (original_key, instance, residue_set)
        all_entries: List[Tuple[str, MotifInstance, Set]] = []
        for mtype, instances in motifs.items():
            for inst in instances:
                rset = _get_residue_set(inst)
                if rset:
                    all_entries.append((mtype, inst, rset))

        # Find strict subset pairs across the whole pool
        to_remove: Set[int] = set()          # id(inst) of instances to drop
        total_subsumed = 0

        n = len(all_entries)
        if n >= 2:
            # Sort by residue-set size descending so we compare smaller
            # against larger first.  This is O(n^2) but n is typically small.
            entries_sorted = sorted(all_entries, key=lambda e: len(e[2]), reverse=True)
            for i in range(n):
                if id(entries_sorted[i][1]) in to_remove:
                    continue
                i_key, i_inst, i_rset = entries_sorted[i]
                for j in range(i + 1, n):
                    if id(entries_sorted[j][1]) in to_remove:
                        continue
                    j_key, j_inst, j_rset = entries_sorted[j]
                    # j is same size or smaller (sorted desc)
                    if j_rset.issubset(i_rset) and j_rset != i_rset:
                        # j is a strict subset of i → remove j, annotate i
                        to_remove.add(id(j_inst))
                        _annotate_overlap(i_inst, j_inst)
                        total_subsumed += 1
                        logger.debug(
                            f"[MERGE] Post-merge: subsumed '{j_key}' "
                            f"({len(j_rset)} residues) into '{i_key}' "
                            f"({len(i_rset)} residues)"
                        )

        if total_subsumed == 0:
            return motifs

        # Rebuild motif dict without subsumed instances
        cleaned: Dict[str, List[MotifInstance]] = {}
        for mtype, instances in motifs.items():
            surviving = [inst for inst in instances if id(inst) not in to_remove]
            if surviving:
                cleaned[mtype] = surviving

        logger.info(
            f"[MERGE] Post-merge subset dedup: removed {total_subsumed} "
            f"subset instance(s)"
        )
        return cleaned
