"""Live PyMOL regression for rmv_super's explicit non-hydrogen fitting."""

from __future__ import annotations

import math
import json
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("RSMVIEWER_ROOT", Path(__file__).resolve().parents[1])).resolve()
sys.path.insert(0, str(ROOT))

from pymol import cmd
from rsmviewer.plugin import __init_plugin__
from rsmviewer import alignment


EXPECTED_MEDOID = "motif_group_SR_007"
EXPECTED_AVERAGE = 0.362
EXPECTED_FINAL_RMSDS = {
    "motif_group_SR_001": 0.242,
    "motif_group_SR_002": 0.428,
    "motif_group_SR_003": 0.158,
    "motif_group_SR_004": 0.772,
    "motif_group_SR_005": 0.575,
    "motif_group_SR_006": 0.146,
    "motif_group_SR_008": 0.215,
}
EXPECTED_ALIGN_MEDOID = "motif_group_SR_008"
EXPECTED_ALIGN_AVERAGE = 2.5269308346
EXPECTED_ALIGN_FINAL_RMSDS = {
    "motif_group_SR_001": 6.961824,
    "motif_group_SR_002": 0.397266,
    "motif_group_SR_003": 0.646945,
    "motif_group_SR_004": 5.697439,
    "motif_group_SR_005": 2.789911,
    "motif_group_SR_006": 0.619200,
    "motif_group_SR_007": 0.575930,
}
TOLERANCE = 0.02


def assert_close(actual, expected, label):
    if not math.isclose(actual, expected, abs_tol=TOLERANCE):
        raise AssertionError(f"{label}: expected {expected}, got {actual}")


def main():
    __init_plugin__(None)
    cif = ROOT / "1s72.cif"
    cmd.do(f"rmv_fetch {cif}")
    cmd.do("rmv_db RNA3DMotifAtlas")
    cmd.do("rmv_select Sarcin-Ricin, 1S72, RNA3DMotifAtlas, as group_SR")
    cmd.do("rmv_create_object group_SR")

    objects = sorted(cmd.get_object_list("group_SR"))
    if len(objects) != 8:
        raise AssertionError(f"Expected 8 group_SR objects, got {objects}")

    align_objects = []
    for index, source in enumerate(objects, 1):
        name = f"_align_baseline_{index}"
        cmd.create(name, source)
        align_objects.append(name)

    align_matrix, align_skipped = alignment.compute_pairwise_rmsd(
        align_objects, method="align",
        atom_selection=alignment.ALIGNMENT_ATOM_SELECTION,
    )
    if align_skipped:
        raise AssertionError(f"Unexpected skipped zero-hydrogen align pairs: {align_skipped}")
    align_medoid_idx, align_averages = alignment.find_medoid(align_matrix)
    if align_objects[align_medoid_idx] != "_align_baseline_8":
        raise AssertionError(f"Expected align medoid _align_baseline_8, got {align_objects[align_medoid_idx]}")
    assert_close(align_averages[align_medoid_idx], EXPECTED_ALIGN_AVERAGE, "align medoid average RMSD")
    align_results = alignment.superimpose_onto_medoid(
        align_objects, align_medoid_idx, method="align",
        atom_selection=alignment.ALIGNMENT_ATOM_SELECTION,
    )
    for obj, rmsd, ok in align_results:
        if obj == "_align_baseline_8":
            continue
        if not ok:
            raise AssertionError(f"Final align failed for {obj}")
        suffix = int(obj.rsplit("_", 1)[-1])
        expected = EXPECTED_ALIGN_FINAL_RMSDS[f"motif_group_SR_{suffix:03d}"]
        assert_close(rmsd, expected, f"final align RMSD for {obj}")

    # The current implementation's pairwise matrix and medoid path.
    matrix, skipped = alignment.compute_pairwise_rmsd(
        objects, method="super", atom_selection=alignment.ALIGNMENT_ATOM_SELECTION
    )
    if skipped:
        raise AssertionError(f"Unexpected skipped zero-hydrogen pairs: {skipped}")
    medoid_idx, averages = alignment.find_medoid(matrix)
    if objects[medoid_idx] != EXPECTED_MEDOID:
        raise AssertionError(f"Expected medoid {EXPECTED_MEDOID}, got {objects[medoid_idx]}")
    assert_close(averages[medoid_idx], EXPECTED_AVERAGE, "medoid average RMSD")

    final_results = alignment.superimpose_onto_medoid(
        objects, medoid_idx, method="super",
        atom_selection=alignment.ALIGNMENT_ATOM_SELECTION,
    )
    for obj, rmsd, ok in final_results:
        if obj == EXPECTED_MEDOID:
            continue
        if not ok:
            raise AssertionError(f"Final superimposition failed for {obj}")
        assert_close(rmsd, EXPECTED_FINAL_RMSDS[obj], f"final RMSD for {obj}")

    # Add real hydrogen atoms to copied objects and capture every fit selection.
    hydrogen_objects = []
    for index, source in enumerate(objects[:3], 1):
        name = f"_hydro_test_{index}"
        cmd.create(name, source)
        cmd.pseudoatom(name, pos=[0.0, 0.0, float(index)], name="H1", elem="H")
        if cmd.count_atoms(f"{name} and hydro") < 1:
            raise AssertionError(f"Failed to add a hydrogen atom to {name}")
        hydrogen_objects.append(name)

    captured_super = []
    captured_align = []
    original_super = cmd.super
    original_align = cmd.align

    def recording_super(mobile, target, *args, **kwargs):
        captured_super.append((mobile, target))
        if cmd.count_atoms(f"({mobile}) and hydro") or cmd.count_atoms(f"({target}) and hydro"):
            raise AssertionError("A cmd.super selection included hydrogen atoms")
        return original_super(mobile, target, *args, **kwargs)

    def recording_align(mobile, target, *args, **kwargs):
        captured_align.append((mobile, target))
        if cmd.count_atoms(f"({mobile}) and hydro") or cmd.count_atoms(f"({target}) and hydro"):
            raise AssertionError("A cmd.align selection included hydrogen atoms")
        return original_align(mobile, target, *args, **kwargs)

    cmd.super = recording_super
    cmd.align = recording_align
    try:
        _, hydrogen_skipped = alignment.compute_pairwise_rmsd(
            hydrogen_objects, method="super",
            atom_selection=alignment.ALIGNMENT_ATOM_SELECTION,
        )
        hydrogen_final = alignment.superimpose_onto_medoid(
            hydrogen_objects, 0, method="super",
            atom_selection=alignment.ALIGNMENT_ATOM_SELECTION,
        )
        _, hydrogen_align_skipped = alignment.compute_pairwise_rmsd(
            hydrogen_objects, method="align",
            atom_selection=alignment.ALIGNMENT_ATOM_SELECTION,
        )
        hydrogen_align_final = alignment.superimpose_onto_medoid(
            hydrogen_objects, 0, method="align",
            atom_selection=alignment.ALIGNMENT_ATOM_SELECTION,
        )
    finally:
        cmd.super = original_super
        cmd.align = original_align

    if hydrogen_skipped:
        raise AssertionError(f"Unexpected skipped hydrogen test pairs: {hydrogen_skipped}")
    expected_calls = len(hydrogen_objects) * (len(hydrogen_objects) - 1) // 2
    expected_calls += len(hydrogen_objects) - 1
    if len(captured_super) != expected_calls or len(captured_align) != expected_calls:
        raise AssertionError(
            f"Expected {expected_calls} super and align calls, got "
            f"{len(captured_super)} and {len(captured_align)}"
        )
    for calls, label in ((captured_super, "cmd.super"), (captured_align, "cmd.align")):
        if any("not hydro" not in mobile or "not hydro" not in target for mobile, target in calls):
            raise AssertionError(f"Unfiltered {label} selection captured: {calls}")
    if not all(ok for _, _, ok in hydrogen_final):
        raise AssertionError(f"Hydrogen final superimposition failed: {hydrogen_final}")
    if hydrogen_align_skipped or not all(ok for _, _, ok in hydrogen_align_final):
        raise AssertionError("Hydrogen final alignment failed")

    cmd.group("direct_align_test", "_hydro_test_*")
    direct_align_calls = []
    original_align = cmd.align

    def direct_recording_align(mobile, target, *args, **kwargs):
        direct_align_calls.append((mobile, target))
        if cmd.count_atoms(f"({mobile}) and hydro") or cmd.count_atoms(f"({target}) and hydro"):
            raise AssertionError("Direct rmv_align passed hydrogen atoms")
        return original_align(mobile, target, *args, **kwargs)

    cmd.align = direct_recording_align
    try:
        cmd.do("rmv_align direct_align_test")
    finally:
        cmd.align = original_align
    if len(direct_align_calls) != expected_calls:
        raise AssertionError(f"Expected {expected_calls} direct align calls, got {len(direct_align_calls)}")

    cmd.group("direct_super_test", "_hydro_test_*")
    direct_super_calls = []
    original_super = cmd.super

    def direct_recording_super(mobile, target, *args, **kwargs):
        direct_super_calls.append((mobile, target))
        if cmd.count_atoms(f"({mobile}) and hydro") or cmd.count_atoms(f"({target}) and hydro"):
            raise AssertionError("Direct rmv_super passed hydrogen atoms")
        return original_super(mobile, target, *args, **kwargs)

    cmd.super = direct_recording_super
    try:
        cmd.do("rmv_super direct_super_test")
    finally:
        cmd.super = original_super
    if len(direct_super_calls) != expected_calls:
        raise AssertionError(f"Expected {expected_calls} direct super calls, got {len(direct_super_calls)}")

    report = {
        "status": "ok",
        "group_SR_objects": len(objects),
        "medoid": objects[medoid_idx],
        "medoid_average_rmsd": round(averages[medoid_idx], 6),
        "align_medoid": align_objects[align_medoid_idx],
        "align_medoid_average_rmsd": round(align_averages[align_medoid_idx], 6),
        "zero_hydrogen_skipped_pairs": len(skipped),
        "zero_hydrogen_align_skipped_pairs": len(align_skipped),
        "hydrogen_super_calls": len(captured_super),
        "hydrogen_align_calls": len(captured_align),
        "direct_align_calls": len(direct_align_calls),
        "direct_super_calls": len(direct_super_calls),
        "hydrogen_test_skipped_pairs": len(hydrogen_skipped),
        "hydrogen_align_test_skipped_pairs": len(hydrogen_align_skipped),
    }
    report_path = Path(os.environ.get("RSM_SUPER_LIVE_REPORT", "/tmp/rmv_super_not_hydro.json"))
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("LIVE_RMV_SUPER_NOT_HYDRO_OK", json.dumps(report), flush=True)


main()
