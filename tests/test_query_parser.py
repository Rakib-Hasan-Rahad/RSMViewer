import unittest

from rsmviewer.database.query_parser import QuerySyntaxError, parse_selection_query


class QueryParserTests(unittest.TestCase):
    def test_parses_alias_and_source_precedence(self):
        query = parse_selection_query(
            "SR, 1S72 and 1FFK, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP"
        )
        self.assertEqual(query.motif, "SARCIN-RICIN")
        self.assertEqual(query.structures, ("1S72", "1FFK"))
        self.assertEqual(query.group, "group_FP")
        self.assertFalse(query.sources.matches(["RNA3DMotifAtlas"]))
        self.assertTrue(query.sources.matches(["RNAMotifScanX"]))
        self.assertTrue(query.matches_structure("1ffk"))

    def test_and_binds_tighter_than_or(self):
        query = parse_selection_query(
            "KT, all, RNA3DMotifAtlas or Rfam and FR3D, as group_KT"
        )
        self.assertTrue(query.sources.matches(["RNA3DMotifAtlas"]))
        self.assertTrue(query.sources.matches(["Rfam", "FR3D"]))
        self.assertFalse(query.sources.matches(["Rfam"]))

    def test_rejects_unknown_source_with_location(self):
        with self.assertRaisesRegex(QuerySyntaxError, "Unknown source 'BGSU'"):
            parse_selection_query("SR, all, BGSU, as group_SR")

    def test_rejects_missing_group(self):
        with self.assertRaisesRegex(QuerySyntaxError, "Expected a group name"):
            parse_selection_query("SR, all, Rfam, group_SR")

    def test_all_matches_any_loaded_structure(self):
        query = parse_selection_query("CL, all, Rfam, as group_CL")
        self.assertTrue(query.matches_structure("1S72"))


if __name__ == "__main__":
    unittest.main()