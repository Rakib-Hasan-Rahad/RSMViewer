import unittest

from rsmviewer.database.source_registry import SourceRegistry


class SourceRegistryTests(unittest.TestCase):
    def setUp(self):
        self.registry = SourceRegistry()

    def test_parses_names_case_insensitively_in_input_order(self):
        self.assertEqual(
            self.registry.parse_names("rna3dmotifatlas,RFAM,fr3d,RnaMotifScanX"),
            ["RNA3DMotifAtlas", "Rfam", "FR3D", "RNAMotifScanX"],
        )

    def test_maps_public_names_to_internal_provider_ids(self):
        self.assertEqual(
            self.registry.get_provider_ids(["RNA3DMotifAtlas", "Rfam"]),
            ["bgsu_api", "rfam_api"],
        )

    def test_rejects_numeric_source_ids(self):
        with self.assertRaisesRegex(ValueError, "Unknown source '3'"):
            self.registry.parse_names("3")

    def test_rejects_removed_source_names(self):
        for name in ("BGSU", "RNAMotifScan", "NoBIAS", "atlas"):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, "Unknown source"):
                    self.registry.parse_names(name)

    def test_rejects_case_insensitive_duplicates(self):
        with self.assertRaisesRegex(ValueError, "Duplicate source 'Rfam'"):
            self.registry.parse_names("Rfam,rfam")

    def test_rejects_empty_items(self):
        for value in ("", "Rfam,", ",Rfam", "Rfam,,FR3D"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "cannot be empty"):
                    self.registry.parse_names(value)


if __name__ == "__main__":
    unittest.main()