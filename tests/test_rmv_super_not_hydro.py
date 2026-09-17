"""Unit tests: rmv_super fits over non-hydrogen atoms; rmv_align is unchanged.

These tests load ``rsmviewer/alignment.py`` in isolation with a fake ``pymol``
module so the exact selection strings passed to ``cmd.super`` / ``cmd.align``
can be captured without a running PyMOL.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path

ALIGNMENT_PATH = Path(__file__).resolve().parents[1] / "rsmviewer" / "alignment.py"


class FakeCmd:
    """Minimal PyMOL ``cmd`` stand-in that records fit calls."""

    def __init__(self):
        self.super_calls = []
        self.align_calls = []
        self.empty_selections = set()

    # object bookkeeping -------------------------------------------------
    def create(self, name, selection=None, *a, **k):
        return None

    def copy(self, dst, src, *a, **k):
        return None

    def delete(self, name, *a, **k):
        return None

    def count_atoms(self, selection, *a, **k):
        return 0 if selection in self.empty_selections else 12

    # fitting ------------------------------------------------------------
    def super(self, mobile, target, *a, **k):
        self.super_calls.append((mobile, target))
        return (0.5, 12, 1, 0.6, 14, 100.0, 8)

    def align(self, mobile, target, *a, **k):
        self.align_calls.append((mobile, target))
        return (0.5, 12, 1, 0.6, 14, 100.0, 8)


def _load_alignment(fake_cmd):
    fake_pymol = types.ModuleType("pymol")
    fake_pymol.cmd = fake_cmd
    sys.modules["pymol"] = fake_pymol
    spec = importlib.util.spec_from_file_location(
        "rsmviewer_alignment_under_test", ALIGNMENT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RmvSuperNotHydroTests(unittest.TestCase):
    def setUp(self):
        self.cmd = FakeCmd()
        self.mod = _load_alignment(self.cmd)
        self.objects = ["m1", "m2", "m3"]

    # 1. pairwise super calls are non-hydrogen ---------------------------
    def test_pairwise_super_uses_not_hydro(self):
        self.mod.compute_pairwise_rmsd(
            self.objects, method="super", atom_selection=self.mod.ALIGNMENT_ATOM_SELECTION
        )
        self.assertTrue(self.cmd.super_calls)
        for mobile, target in self.cmd.super_calls:
            self.assertTrue(mobile.endswith("and not hydro"), mobile)
            self.assertTrue(target.endswith("and not hydro"), target)
        self.assertEqual(self.cmd.align_calls, [])

    # 2. final medoid super calls are non-hydrogen -----------------------
    def test_final_super_uses_not_hydro(self):
        results = self.mod.superimpose_onto_medoid(
            self.objects, 0, method="super",
            atom_selection=self.mod.ALIGNMENT_ATOM_SELECTION,
        )
        self.assertEqual(len(results), 3)
        self.assertTrue(self.cmd.super_calls)
        for mobile, target in self.cmd.super_calls:
            self.assertTrue(mobile.endswith("and not hydro"), mobile)
            self.assertTrue(target.endswith("and not hydro"), target)

    # 3. medoid uses the RMSDs returned from the filtered calls ----------
    def test_medoid_uses_filtered_rmsds(self):
        matrix, skipped = self.mod.compute_pairwise_rmsd(
            self.objects, method="super", atom_selection=self.mod.ALIGNMENT_ATOM_SELECTION
        )
        self.assertEqual(skipped, [])
        medoid_idx, avg = self.mod.find_medoid(matrix)
        self.assertIn(medoid_idx, range(len(self.objects)))
        self.assertTrue(all(v != float("inf") for v in avg))

    # 4. rmv_align uses the shared non-hydrogen selection -----------------
    def test_align_path_unchanged(self):
        self.mod.compute_pairwise_rmsd(
            self.objects, method="align", atom_selection=self.mod.ALIGNMENT_ATOM_SELECTION
        )
        self.mod.superimpose_onto_medoid(
            self.objects, 0, method="align",
            atom_selection=self.mod.ALIGNMENT_ATOM_SELECTION,
        )
        self.assertTrue(self.cmd.align_calls)
        self.assertEqual(self.cmd.super_calls, [])
        for mobile, target in self.cmd.align_calls:
            self.assertTrue(mobile.endswith("and not hydro"), mobile)
            self.assertTrue(target.endswith("and not hydro"), target)

    # 5. rmv_super and rmv_align use their respective PyMOL calls --------
    def test_super_and_align_use_distinct_pymol_calls(self):
        self.mod.compute_pairwise_rmsd(
            self.objects, method="super",
            atom_selection=self.mod.ALIGNMENT_ATOM_SELECTION,
        )
        self.mod.compute_pairwise_rmsd(
            self.objects, method="align",
            atom_selection=self.mod.ALIGNMENT_ATOM_SELECTION,
        )
        self.assertTrue(self.cmd.super_calls)
        self.assertTrue(self.cmd.align_calls)

    # Helper still supports unfiltered calls for unrelated callers. ---------
    def test_no_selection_passes_plain_names(self):
        self.mod.superimpose_onto_medoid(
            self.objects, 0, method="super", atom_selection=None
        )
        self.assertTrue(self.cmd.super_calls)
        for mobile, target in self.cmd.super_calls:
            self.assertNotIn("not hydro", mobile)
            self.assertEqual(mobile, mobile.strip())

    # 7. empty non-hydrogen selections are reported and skipped ----------
    def test_empty_selection_is_skipped(self):
        # Every mobile/target non-hydro selection is empty.
        self.cmd.empty_selections = {
            f"(_medoid_ref_{i}) and not hydro" for i in range(len(self.objects))
        }
        self.cmd.empty_selections.add("(_medoid_work) and not hydro")
        matrix, skipped = self.mod.compute_pairwise_rmsd(
            self.objects, method="super", atom_selection=self.mod.ALIGNMENT_ATOM_SELECTION
        )
        pairs = len(self.objects) * (len(self.objects) - 1) // 2
        self.assertEqual(len(skipped), pairs)
        self.assertEqual(self.cmd.super_calls, [])
        # Skipped pairs are recorded as infinite, never as 0.0 fits.
        for i in range(len(self.objects)):
            for j in range(i + 1, len(self.objects)):
                self.assertEqual(matrix[i][j], float("inf"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
