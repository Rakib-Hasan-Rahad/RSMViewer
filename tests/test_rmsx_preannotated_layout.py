import tempfile
import unittest
from pathlib import Path

from rsmviewer.tools.rmsx_runner import _copy_prebuilt_targets_from_directory


class RmsxPreannotatedLayoutTests(unittest.TestCase):
    def test_expanded_lab_layout_is_discovered(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "rmsx_work_default" / "1s72" / "0"
            root.mkdir(parents=True)
            source = root / "1s72_0.rmsx.in"
            source.write_text("fixture", encoding="utf-8")
            output = Path(directory) / "output"
            result = _copy_prebuilt_targets_from_directory(
                str(Path(directory) / "rmsx_work_default"),
                "1S72",
                str(output),
                ["0"],
            )
            self.assertIn("0", result)
            self.assertEqual(Path(result["0"]).read_text(encoding="utf-8"), "fixture")


if __name__ == "__main__":
    unittest.main()
