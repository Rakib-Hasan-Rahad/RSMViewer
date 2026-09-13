import tempfile
import unittest
from pathlib import Path

from rsmviewer.database.user_annotations.converters import FR3DConverter, RNAMotifScanXConverter


class ExternalConverterTests(unittest.TestCase):
    def test_fr3d_csv_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fr3d.csv"
            path.write_text(
                "Motif order,Motif type,Positions,Sequence,Description\n"
                "1,Hairpin,1S72|1|A|10-12,AGU,fixture\n",
                encoding="utf-8",
            )
            result = FR3DConverter.convert_file(str(path))
            self.assertIn("Hairpin", result)
            self.assertEqual(len(result["Hairpin"][0].residues), 3)

    def test_rmsx_replay_and_pvalue_filter(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.log"
            path.write_text(
                "#fragment_ID\taligned_regions\talignment_score\tP-value\n"
                "1S72_A:10-12\t\t12.0\t0.001\n"
                "1S72_A:20-22\t\t2.0\t0.9\n",
                encoding="utf-8",
            )
            result = RNAMotifScanXConverter.convert_file(
                str(path), "k-turn", apply_filters=True
            )
            self.assertEqual(len(result["K-TURN"]), 1)
            self.assertEqual(result["K-TURN"][0].metadata["p_value"], 0.001)

    def test_rmsx_alignment_report_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "k-turn_consensus.log"
            path.write_text(
                "#  Aligning k-turn_consensus and 1a51_A:30-37_3-10:\n"
                "#  Alignment score:  23.8\n"
                "#  P-value:  1.000e+00\n",
                encoding="utf-8",
            )
            result = RNAMotifScanXConverter.convert_file(
                str(path), "k-turn", apply_filters=False
            )
            self.assertEqual(len(result["K-TURN"]), 1)
            self.assertEqual(result["K-TURN"][0].residues[0][1], 30)


if __name__ == "__main__":
    unittest.main()
