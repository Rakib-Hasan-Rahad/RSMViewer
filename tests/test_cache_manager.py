import tempfile
import unittest
from pathlib import Path

from rsmviewer.database.base_provider import MotifInstance, ResidueSpec
from rsmviewer.database.cache_manager import CacheManager


class CacheManagerTests(unittest.TestCase):
    def test_provenance_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            manager = CacheManager(Path(directory))
            motifs = {
                "SR": [
                    MotifInstance(
                        instance_id="one",
                        motif_id="SR",
                        pdb_id="1S72",
                        residues=[ResidueSpec("A", 1)],
                    )
                ]
            }
            self.assertTrue(
                manager.cache_motifs(
                    "1S72",
                    "rfam_api",
                    motifs,
                    provenance={"endpoint": "https://example.invalid", "release": "test"},
                )
            )
            metadata_path = Path(directory) / "1S72_rfam_api.meta.json"
            metadata = metadata_path.read_text(encoding="utf-8")
            self.assertIn('"endpoint": "https://example.invalid"', metadata)
            restored = manager.get_cached_motifs("1S72", "rfam_api")
            self.assertEqual(restored["SR"][0].instance_id, "one")

    def test_old_metadata_without_provenance_is_readable(self):
        with tempfile.TemporaryDirectory() as directory:
            manager = CacheManager(Path(directory))
            metadata = {
                "pdb_id": "1S72",
                "source": "rfam_api",
                "fetched_at": "2026-09-11T00:00:00",
                "expires_at": "2999-09-11T00:00:00",
                "version": "2.1",
            }
            path = Path(directory) / "1S72_rfam_api.meta.json"
            path.write_text(__import__("json").dumps(metadata), encoding="utf-8")
            self.assertEqual(manager.get_cached_motifs("1S72", "rfam_api"), None)


if __name__ == "__main__":
    unittest.main()