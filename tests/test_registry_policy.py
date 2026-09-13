import tempfile
import unittest
from pathlib import Path

from rsmviewer.database.registry import initialize_registry


class RegistryPolicyTests(unittest.TestCase):
    def test_bundled_sources_are_not_registered_by_default(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory)
            (database / "RNA 3D motif atlas").mkdir()
            (database / "Rfam motif database").mkdir()
            registry = initialize_registry(str(database), enable_api=False)
            self.assertNotIn("atlas", registry.get_provider_ids())
            self.assertNotIn("rfam", registry.get_provider_ids())


if __name__ == "__main__":
    unittest.main()