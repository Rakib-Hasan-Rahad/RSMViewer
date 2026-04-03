"""
Structure Exporter Module for RSMViewer.

Extracts motif instances as mmCIF files with ORIGINAL coordinates
directly from the on-disk CIF file (not PyMOL's internal coordinates,
which may be slightly modified during loading).

Folder structure mirrors the image saver:
    motif_structures/pdb_id/MOTIF_TYPE/TYPE-NO-CHAIN-RESIDUES.cif
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


class MotifStructureExporter:
    """Export motif instances as mmCIF files from the original CIF source."""

    def __init__(self, cmd):
        """
        Initialize the exporter.

        Args:
            cmd: PyMOL cmd module for querying fetch_path
        """
        self.cmd = cmd
        self.logger = logging.getLogger('rsmviewer.structure_exporter')

    # ------------------------------------------------------------------
    # Folder / filename helpers  (mirrors image_saver.py)
    # ------------------------------------------------------------------

    def create_folder_hierarchy(self, pdb_id: str,
                                output_base_dir: Optional[str] = None) -> Path:
        """Create and return the PDB output folder.

        Default base directory: ``<plugin_dir>/../motif_structures/``
        """
        if output_base_dir:
            base_dir = Path(output_base_dir)
        else:
            plugin_dir = Path(__file__).parent
            base_dir = plugin_dir.parent / 'motif_structures'

        pdb_folder = base_dir / pdb_id.lower()
        pdb_folder.mkdir(parents=True, exist_ok=True)
        return pdb_folder

    def create_motif_type_folder(self, pdb_folder: Path,
                                 motif_type: str) -> Path:
        motif_folder = pdb_folder / motif_type.upper()
        motif_folder.mkdir(parents=True, exist_ok=True)
        return motif_folder

    def generate_instance_filename(self, instance_no: int,
                                   motif_details: Dict,
                                   motif_type: str = '') -> str:
        """Build filename: TYPE-NO-CHAIN_RESIDUES.cif"""
        residues = motif_details.get('residues', [])
        chain_residues: Dict[str, List[str]] = {}

        for res in residues:
            if isinstance(res, tuple) and len(res) >= 3:
                _nuc, resi, chain = res[0], res[1], res[2]
                chain_residues.setdefault(chain, []).append(str(resi))

        chain_parts = []
        for chain in sorted(chain_residues.keys()):
            condensed = self._condense_residues(chain_residues[chain])
            chain_parts.append(f"{chain}_{condensed}")

        chain_residue_str = "_".join(chain_parts) if chain_parts else "unknown"

        if motif_type:
            return f"{motif_type}-{instance_no}-{chain_residue_str}.cif"
        return f"instance_{instance_no}_{chain_residue_str}.cif"

    # ------------------------------------------------------------------
    # CIF parsing helpers
    # ------------------------------------------------------------------

    def _find_cif_file(self, pdb_id: str) -> Optional[str]:
        """Locate the original CIF file on disk via PyMOL's fetch_path."""
        fetch_path = "."
        try:
            fp = self.cmd.get("fetch_path")
            if fp:
                fetch_path = str(fp).strip()
        except Exception:
            pass

        for name in [f"{pdb_id.lower()}.cif", f"{pdb_id.upper()}.cif"]:
            candidate = os.path.join(fetch_path, name)
            if os.path.exists(candidate):
                return candidate
        return None

    def _parse_atom_site_header(self, lines: List[str],
                                start_idx: int) -> Tuple[List[str], int, int, int, int]:
        """Parse ``_atom_site.*`` column headers starting at *start_idx*.

        Returns
        -------
        column_names : list[str]
            Full list of ``_atom_site.*`` column names in order.
        auth_asym_col : int
            Column index for auth_asym_id (chain).
        auth_seq_col : int
            Column index for auth_seq_id (residue number).
        data_start : int
            Line index of the first data row.
        num_columns : int
            Total number of columns.
        """
        column_names: List[str] = []
        auth_asym_col = -1
        auth_seq_col = -1
        idx = start_idx

        while idx < len(lines):
            stripped = lines[idx].strip()
            if stripped.startswith('_atom_site.'):
                col_name = stripped.split()[0]
                column_names.append(col_name)
                if col_name == '_atom_site.auth_asym_id':
                    auth_asym_col = len(column_names) - 1
                elif col_name == '_atom_site.auth_seq_id':
                    auth_seq_col = len(column_names) - 1
                idx += 1
            else:
                break

        return column_names, auth_asym_col, auth_seq_col, idx, len(column_names)

    def _build_residue_set(self, motif_details: Dict) -> Set[Tuple[str, str]]:
        """Return a set of (chain, residue_number) tuples for a motif instance."""
        residues = motif_details.get('residues', [])
        pairs: Set[Tuple[str, str]] = set()
        for res in residues:
            if isinstance(res, tuple) and len(res) >= 3:
                _nuc, resi, chain = res[0], res[1], res[2]
                pairs.add((str(chain), str(resi)))
        return pairs

    # ------------------------------------------------------------------
    # Core extraction
    # ------------------------------------------------------------------

    def extract_instance_cif(self, cif_path: str,
                             motif_details: Dict,
                             output_path: Path,
                             pdb_id: str) -> bool:
        """Extract atoms for one motif instance and write a minimal mmCIF.

        Reads the original CIF, filters ``_atom_site`` rows to those
        matching the instance's (chain, residue) pairs, and writes
        a standalone mmCIF.

        Returns True on success.
        """
        residue_set = self._build_residue_set(motif_details)
        if not residue_set:
            self.logger.warning("No residues in motif details — skipping")
            return False

        try:
            with open(cif_path, 'r') as fh:
                all_lines = fh.readlines()
        except Exception as exc:
            self.logger.error(f"Cannot read CIF file: {exc}")
            return False

        # ---- Locate _atom_site loop_ ----
        atom_site_loop_start: Optional[int] = None
        atom_site_header_start: Optional[int] = None

        for i, line in enumerate(all_lines):
            stripped = line.strip()
            if stripped == 'loop_':
                # Check if the next non-blank line starts with _atom_site.
                j = i + 1
                while j < len(all_lines) and not all_lines[j].strip():
                    j += 1
                if j < len(all_lines) and all_lines[j].strip().startswith('_atom_site.'):
                    atom_site_loop_start = i
                    atom_site_header_start = j
                    break

        if atom_site_loop_start is None or atom_site_header_start is None:
            self.logger.error("Could not find _atom_site loop in CIF file")
            return False

        # ---- Parse column headers ----
        (column_names, auth_asym_col, auth_seq_col,
         data_start, num_cols) = self._parse_atom_site_header(
            all_lines, atom_site_header_start)

        if auth_asym_col < 0 or auth_seq_col < 0:
            self.logger.error(
                "CIF missing auth_asym_id or auth_seq_id columns "
                f"(asym={auth_asym_col}, seq={auth_seq_col})")
            return False

        # ---- Filter data rows ----
        matched_rows: List[str] = []
        atom_site_end = data_start

        for i in range(data_start, len(all_lines)):
            stripped = all_lines[i].strip()
            if not stripped or stripped.startswith('#') or stripped.startswith('loop_') or stripped.startswith('_'):
                atom_site_end = i
                break

            tokens = stripped.split()
            if len(tokens) < num_cols:
                # Possibly a continuation / semicolon value — skip
                continue

            chain_val = tokens[auth_asym_col]
            resi_val = tokens[auth_seq_col]
            if (chain_val, resi_val) in residue_set:
                matched_rows.append(all_lines[i])
        else:
            # Reached EOF while still in data
            atom_site_end = len(all_lines)

        if not matched_rows:
            self.logger.warning("No matching atoms found in CIF for this instance")
            return False

        # ---- Write output mmCIF ----
        try:
            with open(output_path, 'w') as out:
                # Data block header
                out.write(f"data_{pdb_id.upper()}_motif\n")
                out.write("#\n")

                # Write key metadata blocks from original
                # (copy _cell, _symmetry, _entity etc. if present)
                self._write_metadata_blocks(all_lines, atom_site_loop_start, out)

                # Write filtered _atom_site block
                out.write("loop_\n")
                for col_name in column_names:
                    out.write(f"{col_name}\n")
                for row in matched_rows:
                    out.write(row)
                out.write("#\n")

            return True

        except Exception as exc:
            self.logger.error(f"Failed to write CIF: {exc}")
            return False

    def _write_metadata_blocks(self, lines: List[str],
                               atom_site_start: int,
                               out) -> None:
        """Copy useful metadata blocks that precede _atom_site.

        Includes _cell, _symmetry, _entity, _struct, _pdbx_struct_assembly,
        _audit_author and similar blocks — everything before the _atom_site
        loop (except the original data_ header which we already wrote).
        """
        i = 0
        # Skip the original data_ line
        while i < atom_site_start:
            stripped = lines[i].strip()
            if stripped.startswith('data_'):
                i += 1
                continue
            # Copy everything else up to the _atom_site loop_ line
            out.write(lines[i])
            i += 1

    # ------------------------------------------------------------------
    # High-level export methods  (parallel to image_saver methods)
    # ------------------------------------------------------------------

    def export_instance(self, motif_folder: Path, instance_no: int,
                        motif_type: str, motif_details: Dict,
                        cif_path: str, pdb_id: str) -> bool:
        """Export a single motif instance to mmCIF.

        Returns True on success.
        """
        filename = self.generate_instance_filename(instance_no, motif_details,
                                                   motif_type)
        filepath = motif_folder / filename

        success = self.extract_instance_cif(cif_path, motif_details,
                                            filepath, pdb_id)
        if success:
            self.logger.info(f"Exported: {filename}")
        return success

    def export_all_motifs(self, loaded_motifs: Dict, pdb_id: str,
                          output_base_dir: Optional[str] = None) -> Dict:
        """Export all loaded motif instances as mmCIF files.

        Returns statistics dict.
        """
        cif_path = self._find_cif_file(pdb_id)
        if not cif_path:
            print(f"\n Original CIF file not found for {pdb_id}.")
            print(f"     Make sure the structure was fetched via rmv_fetch or rmv_load.")
            return {'total_saved': 0, 'total_failed': 0}

        pdb_folder = self.create_folder_hierarchy(pdb_id, output_base_dir)

        stats: Dict = {
            'pdb_id': pdb_id,
            'output_dir': str(pdb_folder),
            'motif_types': {},
            'total_saved': 0,
            'total_failed': 0,
        }

        if not loaded_motifs:
            print("  No motifs loaded to export")
            return stats

        print("=" * 60)
        print(f"  EXPORTING MOTIF STRUCTURES (mmCIF) — {pdb_id}")
        print("=" * 60)
        print(f"  Source CIF: {cif_path}")
        print(f"  Using ORIGINAL coordinates from disk")
        print("-" * 60)

        for motif_type, info in sorted(loaded_motifs.items()):
            motif_folder = self.create_motif_type_folder(pdb_folder, motif_type)
            motif_details = info.get('motif_details', [])
            saved = 0
            failed = 0

            print(f"\n  {motif_type}: Exporting {len(motif_details)} instances...")

            for instance_no, detail in enumerate(motif_details, 1):
                ok = self.export_instance(motif_folder, instance_no,
                                          motif_type, detail,
                                          cif_path, pdb_id)
                if ok:
                    saved += 1
                    stats['total_saved'] += 1
                else:
                    failed += 1
                    stats['total_failed'] += 1

            stats['motif_types'][motif_type] = {
                'saved': saved,
                'failed': failed,
                'total': len(motif_details),
                'folder': str(motif_folder),
            }
            print(f"Exported {saved}/{len(motif_details)} instances")
            if failed:
                print(f"    ⚠️  Failed: {failed}")

        print("\n" + "=" * 60)
        print(f"  SUMMARY")
        print("=" * 60)
        print(f"  Total exported: {stats['total_saved']}")
        if stats['total_failed']:
            print(f"  Total failed:   {stats['total_failed']}")
        print(f"  Output folder:  {pdb_folder}")
        print("=" * 60 + "\n")
        return stats

    def export_motif_type(self, loaded_motifs: Dict, motif_type: str,
                          pdb_id: str,
                          output_base_dir: Optional[str] = None) -> Dict:
        """Export all instances of a specific motif type as mmCIF."""
        cif_path = self._find_cif_file(pdb_id)
        if not cif_path:
            print(f"\n  ❌ Original CIF file not found for {pdb_id}.")
            print(f"     Make sure the structure was fetched via rmv_fetch or rmv_load.")
            return {'saved': 0, 'failed': 0}

        pdb_folder = self.create_folder_hierarchy(pdb_id, output_base_dir)
        motif_type_upper = motif_type.upper()

        stats: Dict = {
            'pdb_id': pdb_id,
            'motif_type': motif_type_upper,
            'output_dir': str(pdb_folder),
            'saved': 0,
            'failed': 0,
        }

        if motif_type_upper not in loaded_motifs:
            print(f"Motif type '{motif_type}' not found")
            avail = ', '.join(sorted(loaded_motifs.keys()))
            print(f"     Available: {avail}")
            return stats

        info = loaded_motifs[motif_type_upper]
        motif_details = info.get('motif_details', [])
        if not motif_details:
            print(f"  No instances for {motif_type_upper}")
            return stats

        motif_folder = self.create_motif_type_folder(pdb_folder, motif_type_upper)

        print("=" * 60)
        print(f"  EXPORTING MOTIF STRUCTURES (mmCIF) — {pdb_id}")
        print("=" * 60)
        print(f"  Source CIF:  {cif_path}")
        print(f"  Motif Type:  {motif_type_upper}")
        print(f"  Instances:   {len(motif_details)}")
        print(f"  Using ORIGINAL coordinates from disk")
        print("-" * 60)

        for instance_no, detail in enumerate(motif_details, 1):
            ok = self.export_instance(motif_folder, instance_no,
                                      motif_type_upper, detail,
                                      cif_path, pdb_id)
            if ok:
                stats['saved'] += 1
            else:
                stats['failed'] += 1

        print("\n" + "=" * 60)
        print(f"  SUMMARY")
        print("=" * 60)
        print(f"  Total exported: {stats['saved']}/{len(motif_details)}")
        if stats['failed']:
            print(f"  Total failed:   {stats['failed']}")
        print(f"  Output folder:  {motif_folder}")
        print("=" * 60 + "\n")
        return stats

    def export_motif_instance(self, loaded_motifs: Dict, motif_type: str,
                              instance_id: int, pdb_id: str,
                              output_base_dir: Optional[str] = None) -> bool:
        """Export a single specific motif instance as mmCIF.

        Returns True on success.
        """
        cif_path = self._find_cif_file(pdb_id)
        if not cif_path:
            print(f"\n  ❌ Original CIF file not found for {pdb_id}.")
            print(f"     Make sure the structure was fetched via rmv_fetch or rmv_load.")
            return False

        pdb_folder = self.create_folder_hierarchy(pdb_id, output_base_dir)
        motif_type_upper = motif_type.upper()

        if motif_type_upper not in loaded_motifs:
            print(f"Motif type '{motif_type}' not found")
            avail = ', '.join(sorted(loaded_motifs.keys()))
            print(f"     Available: {avail}")
            return False

        info = loaded_motifs[motif_type_upper]
        motif_details = info.get('motif_details', [])

        if not motif_details:
            print(f"No instances for {motif_type_upper}")
            return False

        if instance_id < 1 or instance_id > len(motif_details):
            print(f"Instance ID {instance_id} out of range (1-{len(motif_details)})")
            return False

        motif_folder = self.create_motif_type_folder(pdb_folder, motif_type_upper)
        detail = motif_details[instance_id - 1]

        print("=" * 60)
        print(f"  EXPORTING MOTIF STRUCTURE (mmCIF) — {pdb_id}")
        print("=" * 60)
        print(f"  Source CIF:  {cif_path}")
        print(f"  Motif Type:  {motif_type_upper}")
        print(f"  Instance:    {instance_id}/{len(motif_details)}")
        print(f"  Using ORIGINAL coordinates from disk")
        print("-" * 60)

        ok = self.export_instance(motif_folder, instance_id,
                                  motif_type_upper, detail,
                                  cif_path, pdb_id)

        if ok:
            filename = self.generate_instance_filename(instance_id, detail,
                                                       motif_type_upper)
            filepath = motif_folder / filename
            print(f"\n Exported successfully")
            print(f"  Output file:  {filepath}")
        else:
            print(f"\n  ❌ Failed to export instance")
        print("=" * 60 + "\n")
        return ok

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _condense_residues(residue_list: List[str]) -> str:
        """Condense ['1','2','3','5','6'] → '1-3_5-6'."""
        try:
            residues = sorted(int(r) for r in residue_list)
        except (ValueError, TypeError):
            return "_".join(residue_list[:3])

        if not residues:
            return ""

        ranges: List[str] = []
        start = end = residues[0]

        for resi in residues[1:]:
            if resi == end + 1:
                end = resi
            else:
                ranges.append(str(start) if start == end else f"{start}-{end}")
                start = end = resi

        ranges.append(str(start) if start == end else f"{start}-{end}")
        return "_".join(ranges)
