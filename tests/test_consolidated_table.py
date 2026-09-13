import unittest

from rsmviewer.database.base_provider import MotifInstance, ResidueSpec
from rsmviewer.database.consolidated_table import (
    ConsolidatedAnnotationTable,
    normalize_residue_set,
    residue_jaccard,
    stable_motif_id,
)


def motif(instance_id, motif_type, residues, metadata=None):
    return MotifInstance(
        instance_id=instance_id,
        motif_id=motif_type,
        pdb_id="1S72",
        residues=[ResidueSpec(chain, number) for chain, number in residues],
        metadata=metadata or {},
    )


class ConsolidatedTableTests(unittest.TestCase):
    def test_normalization_and_id_are_deterministic(self):
        residues = [ResidueSpec("A", 2), ResidueSpec("A", 1), ResidueSpec("A", 1)]
        normalized = normalize_residue_set(residues)
        self.assertEqual(normalized, (("A", 1, "", 1), ("A", 2, "", 1)))
        self.assertEqual(stable_motif_id("1s72", normalized), stable_motif_id("1S72", normalized))

    def test_jaccard_boundary(self):
        first = (("A", 1, "", 1), ("A", 2, "", 1), ("A", 3, "", 1))
        second = (("A", 1, "", 1), ("A", 2, "", 1), ("A", 3, "", 1), ("A", 4, "", 1), ("A", 5, "", 1))
        self.assertAlmostEqual(residue_jaccard(first, second), 0.6)

    def test_sources_merge_and_hierarchy_is_preserved(self):
        table = ConsolidatedAnnotationTable()
        table.add_annotations("1S72", "RNA3DMotifAtlas", {"SR": [motif("a", "SR", [("A", 1), ("A", 2)])]})
        table.add_annotations(
            "1S72",
            "Rfam",
            {"GNRA": [motif("b", "GNRA", [("A", 1), ("A", 2)], {"hierarchy": ["family", "GNRA"]})]},
        )
        self.assertEqual(len(table.rows), 1)
        row = table.rows[0]
        self.assertEqual(row.source_annotations, {"RNA3DMotifAtlas": ("SR",), "Rfam": ("GNRA",)})
        self.assertEqual(row.source_hierarchy["Rfam"], ("family", "GNRA"))

    def test_different_structures_never_merge(self):
        table = ConsolidatedAnnotationTable()
        source_motifs = {"SR": [motif("a", "SR", [("A", 1), ("A", 2)])]}
        table.add_annotations("1S72", "Rfam", source_motifs)
        table.add_annotations("1FFK", "Rfam", source_motifs)
        self.assertEqual(len(table.rows), 2)

    def test_below_threshold_stays_distinct(self):
        # Partial overlap that is below BOTH the Jaccard threshold (0.60) and
        # the containment threshold (0.80): intersection {A4,A5} -> Jaccard
        # 2/8 = 0.25, containment 2/5 = 0.40. These must remain distinct rows.
        table = ConsolidatedAnnotationTable()
        table.add_annotations("1S72", "FR3D", {"SR": [motif("a", "SR", [("A", 1), ("A", 2), ("A", 3), ("A", 4), ("A", 5)])]})
        table.add_annotations("1S72", "FR3D", {"KT": [motif("b", "KT", [("A", 4), ("A", 5), ("A", 6), ("A", 7), ("A", 8)])]})
        self.assertEqual(len(table.rows), 2)

    def test_nested_annotations_merge_by_containment(self):
        # A tight core fully contained in an extended region (Jaccard 3/6 = 0.5,
        # below threshold) merges because containment is 3/3 = 1.0. This mirrors
        # an Atlas motif core nested inside an extended Rfam annotation.
        table = ConsolidatedAnnotationTable()
        table.add_annotations("1S72", "RNA3DMotifAtlas", {"SR": [motif("a", "SR", [("A", 1), ("A", 2), ("A", 3)])]})
        table.add_annotations(
            "1S72",
            "Rfam",
            {"SR": [motif("b", "SR", [("A", 1), ("A", 2), ("A", 3), ("A", 4), ("A", 5), ("A", 6)])]},
        )
        self.assertEqual(len(table.rows), 1)
        self.assertEqual(
            set(table.rows[0].source_annotations), {"RNA3DMotifAtlas", "Rfam"}
        )

    def test_public_ids_use_five_digit_per_structure_sequence(self):
        table = ConsolidatedAnnotationTable()
        table.add_annotations("1S72", "Rfam", {
            "A": [motif("a", "A", [("A", 10)])],
            "B": [motif("b", "B", [("A", 20)])],
        })
        self.assertEqual([row.motif_id for row in table.rows], ["1S72_00001", "1S72_00002"])


if __name__ == "__main__":
    unittest.main()