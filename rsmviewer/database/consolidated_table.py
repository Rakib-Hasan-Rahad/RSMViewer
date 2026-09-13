"""Deterministic, residue-set indexed annotation table.

The table deliberately keeps source labels side by side. It does not attempt
to decide whether two sources use the same biological definition of a motif.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from .base_provider import MotifInstance, ResidueSpec


ResidueKey = Tuple[str, int, str, int]
DEFAULT_JACCARD_THRESHOLD = 0.60
# Overlap-coefficient threshold for containment-aware merging. Two annotations
# that describe the same motif at different granularity (e.g. a tight Atlas
# core nested inside an extended Rfam region) can have a low Jaccard yet a high
# containment. When the smaller residue set is at least this fraction inside
# the larger, the rows are merged even if their Jaccard is below threshold.
DEFAULT_CONTAINMENT_THRESHOLD = 0.80


def normalize_residue_set(residues: Iterable[ResidueSpec]) -> Tuple[ResidueKey, ...]:
    """Return a sorted, duplicate-free, structure-local residue identity."""
    normalized = {
        (
            str(residue.chain),
            int(residue.residue_number),
            str(residue.insertion_code or ""),
            int(residue.model),
        )
        for residue in residues
    }
    return tuple(sorted(normalized))


def residue_jaccard(
    first: Sequence[ResidueKey], second: Sequence[ResidueKey]
) -> float:
    """Calculate residue-set Jaccard similarity."""
    first_set = set(first)
    second_set = set(second)
    union = first_set | second_set
    return len(first_set & second_set) / len(union) if union else 0.0


def residue_containment(
    first: Sequence[ResidueKey], second: Sequence[ResidueKey]
) -> float:
    """Return the overlap coefficient |A∩B| / min(|A|, |B|).

    This is 1.0 when the smaller residue set is fully contained in the larger
    one, so it captures nested annotations (a motif core inside an extended
    region) that Jaccard would score low.
    """
    first_set = set(first)
    second_set = set(second)
    smaller = min(len(first_set), len(second_set))
    return len(first_set & second_set) / smaller if smaller else 0.0


def stable_motif_id(structure_id: str, residue_set: Sequence[ResidueKey]) -> str:
    """Return the public five-digit motif-ID format placeholder.

    The table assigns the sequence number deterministically when rows are
    materialized; this helper remains as the public ID factory API.
    """
    prefix = "".join(character for character in structure_id.strip().upper() if character.isalnum()) or "LOCAL"
    return f"{prefix}_00001"


@dataclass
class AnnotationRow:
    """One consolidated physical motif fragment."""

    motif_id: str
    structure_id: str
    residue_set: Tuple[ResidueKey, ...]
    source_annotations: Dict[str, Tuple[str, ...]] = field(default_factory=dict)
    source_hierarchy: Dict[str, Tuple[str, ...]] = field(default_factory=dict)
    provenance: Dict[str, Dict[str, str]] = field(default_factory=dict)


class ConsolidatedAnnotationTable:
    """Deterministic table of structure-local residue-set clusters."""

    def __init__(
        self,
        jaccard_threshold: float = DEFAULT_JACCARD_THRESHOLD,
        containment_threshold: float = DEFAULT_CONTAINMENT_THRESHOLD,
        merge_enabled: bool = True,
    ) -> None:
        if not 0.0 < jaccard_threshold <= 1.0:
            raise ValueError("jaccard_threshold must be between 0 and 1")
        if not 0.0 < containment_threshold <= 1.0:
            raise ValueError("containment_threshold must be between 0 and 1")
        self.jaccard_threshold = jaccard_threshold
        self.containment_threshold = containment_threshold
        # When False the table keeps every annotation as its own row and only
        # co-locates byte-for-byte identical residue sets (used at rmv_db load
        # time so overlapping/contained annotations survive until rmv_select).
        self.merge_enabled = merge_enabled
        self._rows: Dict[str, AnnotationRow] = {}

    @property
    def rows(self) -> Tuple[AnnotationRow, ...]:
        """Return rows in deterministic motif-ID order."""
        return tuple(sorted(self._rows.values(), key=self._row_sort_key))

    def clear(self) -> None:
        self._rows.clear()

    def add_annotations(
        self,
        structure_id: str,
        source_name: str,
        motifs: Mapping[str, Sequence[MotifInstance]],
        provenance: Optional[Mapping[str, str]] = None,
    ) -> List[AnnotationRow]:
        """Merge source records into the table without losing source labels."""
        structure = structure_id.strip().upper()
        source = source_name.strip()
        if not structure:
            raise ValueError("structure_id cannot be empty")
        if not source:
            raise ValueError("source_name cannot be empty")

        for motif_type in sorted(motifs):
            for instance in sorted(motifs[motif_type], key=self._instance_sort_key):
                residue_set = normalize_residue_set(instance.residues)
                if not residue_set:
                    continue
                row = self._find_matching_row(structure, residue_set)
                if row is None:
                    row = self._create_row(structure, residue_set)
                    self._rows[row.motif_id] = row
                    self._reindex_structure(structure)

                labels = list(row.source_annotations.get(source, ()))
                if motif_type not in labels:
                    labels.append(motif_type)
                row.source_annotations[source] = tuple(sorted(labels))

                hierarchy = self._hierarchy_for(instance, motif_type)
                existing_hierarchy = list(row.source_hierarchy.get(source, ()))
                for level in hierarchy:
                    if level not in existing_hierarchy:
                        existing_hierarchy.append(level)
                row.source_hierarchy[source] = tuple(existing_hierarchy)

                source_provenance = dict(provenance or {})
                source_provenance.setdefault("retrieved_at", datetime.now(timezone.utc).isoformat())
                row.provenance[source] = source_provenance

        return list(self.rows)

    def get(self, motif_id: str) -> Optional[AnnotationRow]:
        return self._rows.get(motif_id)

    def for_structure(self, structure_id: str) -> Tuple[AnnotationRow, ...]:
        structure = structure_id.strip().upper()
        return tuple(row for row in self.rows if row.structure_id == structure)

    def _create_row(
        self, structure_id: str, residue_set: Tuple[ResidueKey, ...]
    ) -> AnnotationRow:
        return AnnotationRow(
            motif_id=f"{structure_id}__pending_{len(self._rows) + 1}",
            structure_id=structure_id,
            residue_set=residue_set,
        )

    def _reindex_structure(self, structure_id: str) -> None:
        rows = sorted(
            (row for row in self._rows.values() if row.structure_id == structure_id),
            key=lambda row: row.residue_set,
        )
        for row in rows:
            self._rows.pop(row.motif_id, None)
        for index, row in enumerate(rows, start=1):
            new_id = f"{structure_id}_{index:05d}"
            row.motif_id = new_id
            self._rows[new_id] = row

    @staticmethod
    def _row_sort_key(row: AnnotationRow):
        return row.structure_id, row.residue_set

    def _find_matching_row(
        self, structure_id: str, residue_set: Tuple[ResidueKey, ...]
    ) -> Optional[AnnotationRow]:
        # No-merge mode: only fold byte-for-byte identical residue sets together
        # so overlapping/contained annotations remain separate rows.
        if not self.merge_enabled:
            for row in self.rows:
                if row.structure_id == structure_id and row.residue_set == residue_set:
                    return row
            return None
        candidates = []
        for row in self.rows:
            if row.structure_id != structure_id:
                continue
            jaccard = residue_jaccard(row.residue_set, residue_set)
            containment = residue_containment(row.residue_set, residue_set)
            if jaccard >= self.jaccard_threshold or containment >= self.containment_threshold:
                candidates.append((jaccard, containment, row))
        if not candidates:
            return None
        # Prefer the strongest match; break ties deterministically by motif_id.
        best = max(candidates, key=lambda item: (item[0], item[1], ))
        best_score = (best[0], best[1])
        tied = [row for jac, cont, row in candidates if (jac, cont) == best_score]
        return min(tied, key=lambda row: row.motif_id)

    @staticmethod
    def _hierarchy_for(instance: MotifInstance, motif_type: str) -> Tuple[str, ...]:
        raw_hierarchy = instance.metadata.get("hierarchy", ()) if instance.metadata else ()
        if isinstance(raw_hierarchy, str):
            raw_hierarchy = (raw_hierarchy,)
        levels = [str(level) for level in raw_hierarchy if str(level).strip()]
        if motif_type not in levels:
            levels.append(motif_type)
        return tuple(levels)

    @staticmethod
    def _instance_sort_key(instance: MotifInstance) -> Tuple[str, Tuple[ResidueKey, ...]]:
        return instance.instance_id, normalize_residue_set(instance.residues)