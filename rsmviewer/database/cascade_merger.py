"""
Cascade Merger - Priority-ordered, category-aware motif merging

Merges motif datasets from multiple sources using a right-to-left cascade
strategy with strict priority ordering.  The leftmost source in the
user's input is the highest priority.

Category-aware merging:
  Before overlap comparison, every motif type name is reduced to a
  **base category** by stripping trailing numeric suffixes (_1, _2, …)
  and mapping variant spellings to a canonical form (e.g. KINK-TURN →
  K-TURN, SARCIN_RICIN_2 → SARCIN-RICIN).  Jaccard overlap is then
  checked **only within the same base category**, preserving biologically
  distinct motifs that happen to share spatial residues.  Original
  motif-type keys from each source are kept in the output.

Algorithm:
  Given ordered sources [Src1, Src2, Src3]:
    1. result = pairwise_merge(Src2, Src3)   # Src2 has priority
    2. final  = pairwise_merge(Src1, result)  # Src1 has priority
  
  In each pairwise merge:
    - For each motif in the lower-priority set (updater):
      - Compute its base category
      - Check Jaccard overlap against ref motifs of the SAME base category
      - If Jaccard >= threshold: DISCARD updater motif
      - Otherwise: KEEP updater motif (unique within that category)
    - Final = ref motifs + non-overlapping updater motifs

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
# the cascade merger compares motifs **within the same biological category**.
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

    Note: this **replaces** original keys.  The cascade merger uses
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
    Also propagates any existing '_also_found_in' from the discarded instance (for 3+ source cascades).
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
    
    # Propagate any _also_found_in from the discarded instance (3+ source cascade)
    for label in discarded_inst.metadata.get('_also_found_in', []):
        if label and label not in also and label != ref_label:
            also.append(label)
    
    ref_inst.metadata['_also_found_in'] = also


class CascadeMerger:
    """
    Merges motif datasets from multiple sources using cascade strategy.
    
    Key principles:
    - Sources are priority-ordered (first = highest priority)
    - Merge proceeds right-to-left in cascade fashion
    - Motif type names are harmonized to canonical forms before comparison
    - Jaccard overlap is checked within the same canonical category only
    - Overlapping motifs are resolved by keeping the higher-priority version
    - Non-overlapping motifs from lower-priority sources are added
    """
    
    def __init__(self, jaccard_threshold: float = 0.60):
        """
        Initialize the CascadeMerger.
        
        Args:
            jaccard_threshold: Minimum Jaccard similarity to consider
                              two motifs as overlapping (default 0.60)
        """
        self.threshold = jaccard_threshold
    
    def merge_sources(
        self,
        ordered_sources: List[Dict[str, List[MotifInstance]]],
        source_labels: Optional[List[str]] = None,
    ) -> Dict[str, List[MotifInstance]]:
        """
        Merge multiple motif datasets using right-to-left cascade.
        
        Args:
            ordered_sources: List of motif dicts, ordered by priority
                            (index 0 = highest priority).
                            Each dict maps motif_type -> [MotifInstance, ...]
            source_labels: Optional labels for logging (e.g., ['Source 1', 'Source 3'])
        
        Returns:
            Merged motif dict: motif_type -> [MotifInstance, ...]
        """
        if not ordered_sources:
            return {}
        
        if len(ordered_sources) == 1:
            return ordered_sources[0]
        
        labels = source_labels or [f"Source_{i}" for i in range(len(ordered_sources))]
        
        logger.info(f"[CASCADE] Merging {len(ordered_sources)} sources: {labels}")
        
        # Right-to-left cascade: merge from the end
        # Start with the last (lowest priority) source
        result = ordered_sources[-1]
        logger.info(
            f"[CASCADE] Starting with {labels[-1]}: "
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
                f"[CASCADE] After merging {ref_label} (priority) + {updater_label}: "
                f"{self._count_motifs(result)} motifs in {len(result)} categories"
            )
        
        logger.info(
            f"[CASCADE] Final result: "
            f"{self._count_motifs(result)} motifs in {len(result)} categories"
        )
        
        # Final pass: remove subset instances within the same category.
        # The pairwise merge only checks updater-vs-ref.  Two instances from
        # the SAME source (both ref) are never compared during the cascade,
        # so a small instance can survive even though it is a complete subset
        # of a larger instance from the same source.
        result = self._deduplicate_subsets(result)
        
        return result
    
    def _pairwise_merge(
        self,
        ref: Dict[str, List[MotifInstance]],
        updater: Dict[str, List[MotifInstance]],
        ref_label: str = "ref",
        updater_label: str = "updater",
    ) -> Dict[str, List[MotifInstance]]:
        """
        Merge two motif datasets with ref having priority over updater.

        **Category-aware**: Jaccard overlap is only checked between motifs
        whose *base category* matches (trailing numeric suffixes like ``_1``,
        ``_2`` are stripped, then variant names are harmonized).  This means
        ``SARCIN_RICIN_1`` and ``SARCIN-RICIN`` are compared (both base =
        ``SARCIN-RICIN``), but ``K-TURN`` and ``C-LOOP`` are never compared.

        Original motif-type keys from each source are **preserved** in the
        output so that the priority source's naming convention is respected.
        """
        # --- Start with all ref motifs under their original (uppercased) keys ---
        merged: Dict[str, List[MotifInstance]] = {}
        for mtype, instances in ref.items():
            key = mtype.upper()
            if key in merged:
                merged[key].extend(instances)
            else:
                merged[key] = list(instances)

        # Pre-compute ref residue sets grouped by **base category**
        # Each entry stores (original_key, instance, residue_set)
        ref_by_cat: Dict[str, List[Tuple[str, MotifInstance, Set]]] = {}
        for mtype, instances in ref.items():
            cat = _get_base_category(mtype)
            bucket = ref_by_cat.setdefault(cat, [])
            for inst in instances:
                rset = _get_residue_set(inst)
                if rset:
                    bucket.append((mtype.upper(), inst, rset))

        # Check each updater motif against same-base-category ref motifs
        stats = {'kept': 0, 'discarded': 0, 'subsumed': 0, 'replaced': 0}

        for u_type, u_instances in updater.items():
            u_key = u_type.upper()
            u_cat = _get_base_category(u_type)
            same_cat_refs = ref_by_cat.get(u_cat, [])

            for u_inst in u_instances:
                u_rset = _get_residue_set(u_inst)
                if not u_rset:
                    # No residues – keep it anyway
                    merged.setdefault(u_key, []).append(u_inst)
                    stats['kept'] += 1
                    continue

                # --- Phase 1: Standard Jaccard overlap ---
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
                        break  # Found overlap, skip remaining

                if is_overlap:
                    stats['discarded'] += 1
                    if best_ref_inst is not None:
                        _annotate_overlap(best_ref_inst, u_inst)
                    logger.debug(
                        f"[CASCADE] Discarded {updater_label} motif "
                        f"'{u_key}/{u_inst.instance_id}' - overlaps with "
                        f"{ref_label} category '{u_cat}' (Jaccard={best_j:.3f})"
                    )
                    continue

                # --- Phase 2: Asymmetric containment (subset check) ---
                # If standard Jaccard < threshold, one set may still be a
                # complete subset of the other (e.g. 2 residues inside 11).
                # When |A ∩ B| / |smaller| = 1.0 the smaller set adds no new
                # information — keep the larger, discard the smaller.
                containment_action = None

                for _r_key, r_inst, r_rset in same_cat_refs:
                    if not r_rset:
                        continue
                    if u_rset.issubset(r_rset):
                        # Updater fully contained in ref → discard updater
                        containment_action = ('discard_updater', r_inst, _r_key, r_rset)
                        break
                    elif r_rset.issubset(u_rset):
                        # Ref fully contained in updater → replace ref with updater
                        containment_action = ('replace_ref', r_inst, _r_key, r_rset)
                        break

                if containment_action:
                    action, c_ref_inst, c_ref_key, c_ref_rset = containment_action
                    if action == 'discard_updater':
                        stats['subsumed'] += 1
                        _annotate_overlap(c_ref_inst, u_inst)
                        logger.debug(
                            f"[CASCADE] Subsumed {updater_label} motif "
                            f"'{u_key}/{u_inst.instance_id}' - subset of "
                            f"{ref_label} '{u_cat}' ({len(u_rset)}/{len(c_ref_rset)} residues)"
                        )
                    else:
                        # Replace smaller ref with larger updater's residues,
                        # but keep ref's KEY (priority naming convention).
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
                        # (preserves priority source naming convention)
                        merged.setdefault(c_ref_key, []).append(u_inst)
                        # Build combined source label (e.g. "NoBIAS + RMSX")
                        # Cannot use _annotate_overlap here because we need
                        # to relabel _source_label to the ref's (priority)
                        # label.  _annotate_overlap would add the ref's label
                        # to _also_found_in, and the subsequent relabel would
                        # make _source_label == _also_found_in entry, which
                        # the display layer deduplicates → losing the updater
                        # attribution entirely.  Instead, manually build
                        # _also_found_in with the updater's original label.
                        upd_label_str = u_inst.metadata.get('_source_label', '')
                        ref_label_str = c_ref_inst.metadata.get('_source_label', '')
                        # Set primary label to ref's (priority) label
                        if ref_label_str:
                            u_inst.metadata['_source_label'] = ref_label_str
                        # _also_found_in: updater's original label + propagated
                        also = list(u_inst.metadata.get('_also_found_in', []))
                        if upd_label_str and upd_label_str != ref_label_str and upd_label_str not in also:
                            also.append(upd_label_str)
                        for lbl in c_ref_inst.metadata.get('_also_found_in', []):
                            if lbl and lbl != ref_label_str and lbl not in also:
                                also.append(lbl)
                        u_inst.metadata['_also_found_in'] = also
                        # Update ref_by_cat so future comparisons use the larger set
                        cat_list = ref_by_cat.get(u_cat, [])
                        for idx, (rk, ri, rs) in enumerate(cat_list):
                            if ri is c_ref_inst:
                                cat_list[idx] = (c_ref_key, u_inst, u_rset)
                                break
                        logger.debug(
                            f"[CASCADE] Replaced {ref_label} motif "
                            f"'{c_ref_key}' with larger {updater_label} "
                            f"'{u_key}/{u_inst.instance_id}' "
                            f"({len(c_ref_rset)} → {len(u_rset)} residues)"
                        )
                else:
                    # Truly unique – keep under updater's original key
                    merged.setdefault(u_key, []).append(u_inst)
                    stats['kept'] += 1

        logger.info(
            f"[CASCADE] Pairwise merge {ref_label} + {updater_label}: "
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
        instance in the **same base category**.

        This catches subset pairs that the pairwise merge cannot see — e.g.
        two instances from the same source where one has 2 residues and the
        other has 11 residues that fully contain the first two.

        The larger (superset) instance is kept and annotated with the
        smaller instance's source info.  If both sets are identical this
        is a no-op (within-source exact dedup already ran).
        """
        # Group ALL instances across keys by base category
        # Each entry: (original_key, instance, residue_set, list_ref)
        by_cat: Dict[str, List[Tuple[str, MotifInstance, Set]]] = {}
        for mtype, instances in motifs.items():
            cat = _get_base_category(mtype)
            bucket = by_cat.setdefault(cat, [])
            for inst in instances:
                rset = _get_residue_set(inst)
                if rset:
                    bucket.append((mtype, inst, rset))

        # For each category, find strict subset pairs
        to_remove: Set[int] = set()          # id(inst) of instances to drop
        total_subsumed = 0

        for cat, entries in by_cat.items():
            n = len(entries)
            if n < 2:
                continue
            # Sort by residue-set size descending so we compare smaller
            # against larger first.  This is O(n^2) but n is typically small.
            entries_sorted = sorted(entries, key=lambda e: len(e[2]), reverse=True)
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
                            f"[CASCADE] Post-merge: subsumed '{j_key}' "
                            f"({len(j_rset)} residues) into '{i_key}' "
                            f"({len(i_rset)} residues) in category '{cat}'"
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
            f"[CASCADE] Post-merge subset dedup: removed {total_subsumed} "
            f"subset instance(s)"
        )
        return cleaned
