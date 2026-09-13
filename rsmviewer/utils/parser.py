"""
RSMViewer - Parser Module
Handles parsing of PDB/mmCIF filenames (PDB ID extraction) and selection formatting.
"""

import os
import logging

logger = logging.getLogger('RMV')


class PDBParser:
    """Simple parser for PDB metadata (minimal - mostly handled by PyMOL)."""
    
    @staticmethod
    def extract_pdb_id(filepath_or_id):
        """
        Extract PDB ID from file path or return if already a PDB ID.
        
        Args:
            filepath_or_id (str): Either a PDB ID or file path
        
        Returns:
            str: PDB ID (4 characters) or None if invalid
        """
        # Check if it's already a PDB ID (4 characters, alphanumeric)
        if len(filepath_or_id) == 4 and filepath_or_id.isalnum():
            return filepath_or_id.upper()
        
        # Try to extract from filename
        filename = os.path.basename(filepath_or_id)
        # PDB files often named like "1s72.pdb" or "1S72.cif"
        if filename:
            name_without_ext = os.path.splitext(filename)[0]
            if len(name_without_ext) >= 4:
                potential_id = name_without_ext[:4]
                if potential_id.isalnum():
                    return potential_id.upper()
        
        return None
    
    @staticmethod
    def is_valid_pdb_id(pdb_id):
        """
        Check if a string is a valid PDB ID.
        
        Args:
            pdb_id (str): Potential PDB ID
        
        Returns:
            bool: True if valid format
        """
        if not isinstance(pdb_id, str):
            return False
        return len(pdb_id) == 4 and pdb_id.isalnum()


class SelectionParser:
    """Parser for creating PyMOL selection strings from residue data."""
    
    @staticmethod
    def create_selection_string(chain, residues, structure_name=None, use_segi=False):
        """
        Create a PyMOL selection string from chain and residue numbers.
        
        Layer 2 Chain Validation: When structure_name is provided, validates
        that the chain exists in the loaded structure before building the
        selection. Logs a warning if the chain is not found.
        
        Args:
            chain (str): Chain identifier
            residues (list): List of residue numbers
            structure_name (str, optional): PyMOL object name for chain validation
            use_segi (bool): If True, use 'segi' field instead of 'chain' (Layer 3 fallback)
        
        Returns:
            str: PyMOL selection string (e.g., "chain A and resi 77-82")
        """
        if not residues:
            return None
        
        # Layer 2: Validate chain exists in structure
        if structure_name:
            try:
                from pymol import cmd as pymol_cmd
                actual_chains = pymol_cmd.get_chains(structure_name)
                if actual_chains and chain not in actual_chains:
                    logger.warning(
                        f"Chain '{chain}' not found in '{structure_name}'. "
                        f"Available chains: {actual_chains}"
                    )
            except Exception:
                pass  # PyMOL not available or structure not loaded yet
        
        residues = sorted(residues)
        chain_keyword = "segi" if use_segi else "chain"
        # Build compact range notation: group contiguous residues into
        # ranges (e.g. 1805-1806+2874-2875) instead of a single span
        # from first to last which would include unrelated residues.
        ranges = []
        i = 0
        while i < len(residues):
            start = residues[i]
            end = start
            while i + 1 < len(residues) and residues[i + 1] == end + 1:
                i += 1
                end = residues[i]
            ranges.append(f"{start}-{end}" if end != start else str(start))
            i += 1
        resi_expr = "+".join(ranges)
        selection = f"{chain_keyword} {chain} and resi {resi_expr}"
        return selection
    
    @staticmethod
    def create_detailed_selection(chain, residues):
        """
        Create a detailed PyMOL selection string listing all residues.
        
        Args:
            chain (str): Chain identifier
            residues (list): List of residue numbers
        
        Returns:
            str: Detailed PyMOL selection string
        """
        if not residues:
            return None
        
        residue_list = "+".join([f"resi {r}" for r in sorted(residues)])
        selection = f"chain {chain} and ({residue_list})"
        return selection


def validate_motif_data(motif_entry):
    """
    Validate a motif entry has required fields.
    
    Args:
        motif_entry (dict): Motif dictionary to validate
    
    Returns:
        bool: True if valid
    """
    required_fields = ['chain', 'residues', 'motif_id']
    return all(field in motif_entry for field in required_fields) and \
           isinstance(motif_entry.get('residues'), list) and \
           len(motif_entry.get('residues', [])) > 0
