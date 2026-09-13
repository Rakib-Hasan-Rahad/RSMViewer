import tempfile
import unittest
from pathlib import Path

from rsmviewer.database.motif_hierarchy_cache import MotifHierarchyCache


class HierarchyCacheSourceKeyTests(unittest.TestCase):
    def test_records_canonical_textual_source_key(self):
        with tempfile.TemporaryDirectory() as directory:
            cache = MotifHierarchyCache(str(Path(directory) / "hierarchy.sqlite3"))
            cache.record("1S72", "A:1", "Rfam", "Rfam", "GNRA", 0)
            rows = cache.get_hierarchy_for_pdb("1S72")["A:1"]
            self.assertEqual(rows[0]["source_key"], "Rfam")

    def test_schema_contains_no_integer_source_column(self):
        with tempfile.TemporaryDirectory() as directory:
            cache = MotifHierarchyCache(str(Path(directory) / "hierarchy.sqlite3"))
            columns = {
                row[1]
                for row in cache._conn.execute("PRAGMA table_info(motif_hierarchy)")
            }
            self.assertIn("source_key", columns)
            self.assertNotIn("source_id", columns)


if __name__ == "__main__":
    unittest.main()