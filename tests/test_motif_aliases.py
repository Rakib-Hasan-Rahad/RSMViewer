import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rsmviewer.database.motif_aliases import (
    canonical_motif,
    family_from_text,
    is_known_family,
    labels_match_motif,
)


class MotifAliasTests(unittest.TestCase):
    def test_abbreviations_resolve(self):
        self.assertEqual(canonical_motif("SR"), "SARCIN-RICIN")
        self.assertEqual(canonical_motif("KT"), "K-TURN")
        self.assertEqual(canonical_motif("CL"), "C-LOOP")
        self.assertEqual(canonical_motif("EL"), "E-LOOP")

    def test_spelling_and_separator_variants(self):
        for variant in ("sarcin-ricin", "sarcin_ricin", "SARCINRICIN", "Sarcin Ricin"):
            self.assertEqual(canonical_motif(variant), "SARCIN-RICIN")
        for variant in ("k-turn", "kink-turn", "kink turn", "KTURN"):
            self.assertEqual(canonical_motif(variant), "K-TURN")

    def test_family_index_suffix_stripped(self):
        self.assertEqual(canonical_motif("sarcin-ricin-1"), "SARCIN-RICIN")
        self.assertEqual(canonical_motif("sarcin-ricin-2"), "SARCIN-RICIN")
        self.assertEqual(canonical_motif("k-turn-2"), "K-TURN")

    def test_kturn_and_reverse_kturn_are_distinct(self):
        self.assertEqual(canonical_motif("reverse-kturn"), "REVERSE-K-TURN")
        self.assertNotEqual(canonical_motif("reverse-kturn"), canonical_motif("k-turn"))
        # A K-TURN query must NOT match a REVERSE-K-TURN label.
        self.assertFalse(labels_match_motif("KT", ["REVERSE-K-TURN"]))
        self.assertTrue(labels_match_motif("KT", ["K-TURN"]))

    def test_labels_match_across_sources(self):
        # Atlas, Rfam, and RMSX spellings of the same family all match "SR".
        self.assertTrue(labels_match_motif("SR", ["Sarcin-Ricin"]))
        self.assertTrue(labels_match_motif("SR", ["sarcin-ricin-1"]))
        self.assertTrue(labels_match_motif("SR", ["sarcin-ricin"]))

    def test_unknown_compound_label_substring_fallback(self):
        # A compound/novel label still matches via the guarded substring path.
        self.assertTrue(labels_match_motif("SR", ["sarcin-ricin core region"]))

    def test_short_codes_do_not_substring_bleed(self):
        # Short canonical codes must not match unrelated words containing them.
        self.assertFalse(labels_match_motif("HL", ["THLXYZ"]))
        self.assertFalse(labels_match_motif("IL", ["SILENCER"]))
        self.assertFalse(labels_match_motif("J3", ["MAJ3OR"]))
        # But the real family labels still match.
        self.assertTrue(labels_match_motif("HL", ["Hairpin Loop (HL)"]))
        self.assertTrue(labels_match_motif("IL", ["Internal Loop (IL)"]))
        self.assertTrue(labels_match_motif("J3", ["3-way Junction (J3)"]))

    def test_no_alias_collisions(self):
        # Importing the module registers every alias; a collision raises at
        # import time, so a successful import already proves uniqueness. This
        # asserts the registry is populated and internally consistent.
        from rsmviewer.database.motif_aliases import _ALIAS_TO_CANONICAL, _alnum
        self.assertEqual(_ALIAS_TO_CANONICAL[_alnum("SR")], "SARCIN-RICIN")
        self.assertEqual(_ALIAS_TO_CANONICAL[_alnum("PK")], "PSEUDOKNOT")

    def test_is_known_family(self):
        self.assertTrue(is_known_family("SR"))
        self.assertTrue(is_known_family("reverse-kturn"))
        self.assertFalse(is_known_family("totally-made-up-motif"))

    def test_enriched_atlas_families_are_known(self):
        # Atlas descriptive families are now recognized, so short abbreviations
        # do not falsely match them by substring (e.g. "AG tHS outside loop").
        self.assertTrue(is_known_family("AG THS OUTSIDE LOOP"))
        self.assertFalse(labels_match_motif("EL", ["AG tHS outside loop"]))
        self.assertTrue(is_known_family("Anticodon loop related"))
        self.assertTrue(is_known_family("Major groove platform"))

    def test_fr3d_free_text_family(self):
        # FR3D query names embed the family loosely.
        self.assertEqual(family_from_text("geometric_5_sarcin_ricin"), "SARCIN-RICIN")
        self.assertEqual(family_from_text("geometric_3_sarcin3geometric"), "SARCIN-RICIN")
        self.assertEqual(family_from_text("reverse_kturn_query"), "REVERSE-K-TURN")
        self.assertEqual(family_from_text("kink_turn_5nt"), "K-TURN")
        self.assertIsNone(family_from_text("geometric_2_cWW_in_DNA"))

    def test_fr3d_free_text_matching_keeps_reverse_distinct(self):
        # FR3D names are passed as free_text_labels; reverse-kturn must not
        # match a plain K-TURN query.
        self.assertTrue(labels_match_motif("SR", [], ["geometric_5_sarcin_ricin"]))
        self.assertTrue(labels_match_motif("KT", [], ["kink_turn_5nt"]))
        self.assertFalse(labels_match_motif("KT", [], ["reverse_kturn_query"]))


if __name__ == "__main__":
    unittest.main()
