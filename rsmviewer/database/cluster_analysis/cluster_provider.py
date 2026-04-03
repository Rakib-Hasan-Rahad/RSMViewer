"""
RSMViewer - Cluster Provider
Parses cluster definition CSV files for motif-based cluster analysis.

CSV format (one cluster per line):
    cluster_name,PDBID_CHAIN:range1_range2_...,PDBID_CHAIN:range1_range2_...,...

Each underscore-separated range (e.g. 2552-2555) is a separate motif region.
"""

from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, field


@dataclass
class ClusterEntry:
    """A single PDB entry within a cluster."""
    pdb_id: str
    chain: str
    regions: List[tuple]  # list of (start, end) residue ranges
    index: int = 0        # disambiguator when same pdb+chain appears multiple times

    @property
    def label(self) -> str:
        """Short label for display, e.g. '1S72_0'."""
        return f"{self.pdb_id}_{self.chain}"

    def all_residue_numbers(self) -> List[int]:
        """Expand all regions into individual residue numbers."""
        nums = []
        for start, end in self.regions:
            nums.extend(range(start, end + 1))
        return sorted(set(nums))

    def region_strings(self) -> List[str]:
        """Return list of 'start-end' strings for display."""
        return [f"{s}-{e}" for s, e in self.regions]


class ClusterProvider:
    """Load and query cluster definitions from CSV files."""

    def __init__(self):
        self._clusters: Dict[str, List[ClusterEntry]] = {}
        self._source_path: Optional[Path] = None

    # ------------------------------------------------------------------ #
    # Loading
    # ------------------------------------------------------------------ #

    def load_csv(self, filepath: str) -> bool:
        """
        Parse a cluster CSV file.

        Returns True on success, False on error.
        """
        fpath = Path(filepath)
        if not fpath.exists():
            print(f"  [cluster] Error: file not found: {filepath}")
            return False

        self._clusters.clear()
        self._source_path = fpath

        for line_no, raw_line in enumerate(fpath.read_text().splitlines(), start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = [p.strip() for p in line.split(",")]
            # Remove empty trailing parts (trailing comma)
            while parts and not parts[-1]:
                parts.pop()

            if len(parts) < 2:
                print(f"  [cluster] Warning: skipping line {line_no} "
                      f"(need at least cluster name + 1 entry)")
                continue

            cluster_name = parts[0]
            entries: List[ClusterEntry] = []

            # Track duplicates for index assignment
            seen_labels: Dict[str, int] = {}

            for token in parts[1:]:
                entry = self._parse_entry(token, line_no)
                if entry is None:
                    continue

                label = entry.label
                if label in seen_labels:
                    seen_labels[label] += 1
                    entry.index = seen_labels[label]
                else:
                    seen_labels[label] = 0

                entries.append(entry)

            if entries:
                self._clusters[cluster_name] = entries

        return True

    def load_bundled(self) -> bool:
        """Load the bundled ML_train_Motif_input.csv from this package."""
        bundled = Path(__file__).parent / "ML_train_Motif_input.csv"
        return self.load_csv(str(bundled))

    # ------------------------------------------------------------------ #
    # Queries
    # ------------------------------------------------------------------ #

    def get_cluster_names(self) -> List[str]:
        """Return list of loaded cluster names."""
        return list(self._clusters.keys())

    def get_cluster(self, name: str) -> Optional[List[ClusterEntry]]:
        """Return entries for a given cluster, or None if not found."""
        return self._clusters.get(name)

    def get_cluster_summary(self) -> Dict[str, int]:
        """Return {cluster_name: entry_count} for all clusters."""
        return {name: len(entries) for name, entries in self._clusters.items()}

    @property
    def source_path(self) -> Optional[Path]:
        return self._source_path

    @property
    def is_loaded(self) -> bool:
        return len(self._clusters) > 0

    # ------------------------------------------------------------------ #
    # Internal parsing
    # ------------------------------------------------------------------ #

    @staticmethod
    def _parse_entry(token: str, line_no: int) -> Optional[ClusterEntry]:
        """
        Parse a single entry like '1S72_0:2552-2555_2580-2582_2596-2602'.

        Format: PDBID_CHAIN:range1_range2_...
        """
        if ":" not in token:
            print(f"  [cluster] Warning: skipping malformed entry on "
                  f"line {line_no}: '{token}' (missing ':')")
            return None

        pdb_chain, range_part = token.split(":", 1)

        # Split PDB ID and chain: PDB is always 4 chars
        if len(pdb_chain) < 5 or "_" not in pdb_chain:
            print(f"  [cluster] Warning: skipping entry on line {line_no}: "
                  f"'{token}' (cannot parse PDB_CHAIN)")
            return None

        # PDB ID is first 4 characters, chain is everything after the first underscore
        underscore_pos = pdb_chain.index("_")
        pdb_id = pdb_chain[:underscore_pos].upper()
        chain = pdb_chain[underscore_pos + 1:]

        # Parse residue ranges
        regions = []
        for rng in range_part.split("_"):
            rng = rng.strip()
            if not rng:
                continue
            if "-" in rng:
                parts = rng.split("-", 1)
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    regions.append((start, end))
                except ValueError:
                    print(f"  [cluster] Warning: bad range '{rng}' in "
                          f"entry '{token}' on line {line_no}")
            else:
                try:
                    num = int(rng)
                    regions.append((num, num))
                except ValueError:
                    print(f"  [cluster] Warning: bad value '{rng}' in "
                          f"entry '{token}' on line {line_no}")

        if not regions:
            print(f"  [cluster] Warning: no valid ranges in entry "
                  f"'{token}' on line {line_no}")
            return None

        return ClusterEntry(pdb_id=pdb_id, chain=chain, regions=regions)
