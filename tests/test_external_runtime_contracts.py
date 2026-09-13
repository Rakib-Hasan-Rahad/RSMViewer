import os
import tempfile
import unittest
from pathlib import Path

from rsmviewer.tools.runtime_layout import (
    discover_rmsx_runtime_dir,
    get_runtime_platform_dir,
    validate_fr3d_checkout,
)


class ExternalRuntimeContractTests(unittest.TestCase):
    def test_missing_fr3d_checkout_is_actionable(self):
        result = validate_fr3d_checkout("/path/that/does/not/exist")
        self.assertFalse(result["ok"])
        self.assertIn("checkout-does-not-exist", result["missing"])

    def test_missing_rmsx_runtime_returns_empty(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = str(Path(directory) / "missing")
            self.assertEqual(discover_rmsx_runtime_dir(missing, [missing]), "")

    def test_explicit_rmsx_runtime_wins(self):
        with tempfile.TemporaryDirectory() as directory:
            old_value = os.environ.get("RSMVIEWER_RMSX_RUNTIME_DIR")
            try:
                os.environ["RSMVIEWER_RMSX_RUNTIME_DIR"] = directory
                self.assertEqual(discover_rmsx_runtime_dir("", []), str(Path(directory).resolve()))
            finally:
                if old_value is None:
                    os.environ.pop("RSMVIEWER_RMSX_RUNTIME_DIR", None)
                else:
                    os.environ["RSMVIEWER_RMSX_RUNTIME_DIR"] = old_value

    def test_platform_directory_is_known(self):
        self.assertIn(get_runtime_platform_dir(), {"macos-arm64", "macos-x86_64", "linux-x86_64", "windows-x86_64"})


if __name__ == "__main__":
    unittest.main()
