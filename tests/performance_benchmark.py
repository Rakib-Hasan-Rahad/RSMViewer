"""Reproducible cached consolidated-table benchmark.

Run with:
    python3 tests/performance_benchmark.py

The benchmark isolates annotation consolidation from one-time network work.
PyMOL medoid timings require a PyMOL session and are reported separately by
`tests/pymol_multistructure_smoke.py`.
"""

from __future__ import annotations

import platform
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rsmviewer.database.base_provider import MotifInstance, ResidueSpec
from rsmviewer.database.consolidated_table import ConsolidatedAnnotationTable


def build_instances(count: int):
    return [
        MotifInstance(
            instance_id=f"benchmark_{index}",
            motif_id="SR",
            pdb_id="BENCH",
            residues=[
                ResidueSpec("A", index * 3 + offset)
                for offset in range(3)
            ],
        )
        for index in range(count)
    ]


def main() -> None:
    count = 1000
    instances = build_instances(count)
    table = ConsolidatedAnnotationTable()
    started = time.perf_counter()
    table.add_annotations("BENCH", "RNA3DMotifAtlas", {"SR": instances})
    elapsed = time.perf_counter() - started
    print(f"hardware={platform.platform()}")
    print(f"python={platform.python_version()}")
    print("workload=cached consolidated annotation insertion")
    print(f"motifs={count}")
    print(f"rows={len(table.rows)}")
    print(f"seconds={elapsed:.6f}")
    if len(table.rows) != count:
        raise RuntimeError(f"Expected {count} deterministic rows, got {len(table.rows)}")


if __name__ == "__main__":
    main()
