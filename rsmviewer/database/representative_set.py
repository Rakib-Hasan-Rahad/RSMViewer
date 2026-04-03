"""
Representative Set Loader - NR List CSV Parser

Parses the BGSU Non-Redundant (NR) list CSV to provide chain-specific
PDB-to-Representative mappings. Used by the homolog enricher to find
the representative structure for a given PDB chain.

NR List CSV Format:
  Column 1: Equivalence class ID (e.g., "NR_all_26150.5")
  Column 2: Representative IFE (e.g., "4V9F|1|0")
  Column 3: All members (comma-separated IFEs, e.g., "4V9F|1|0,1S72|1|0,...")

IFE Format: PDB_ID|Model|Chain
  - Chain can be single char ('A'), digit ('0'), multi-char ('hB', '1x')
  - Some entries have '+' for multi-chain IFEs (e.g., "4V9F|1|A+4V9F|1|B")

Author: CBB Lab
Version: 1.0.0
"""

from __future__ import annotations

import csv
import os
import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# Default path to the NR list CSV (bundled with the plugin)
_DEFAULT_CSV_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'nrlist_4.24_all.csv'
)


class RepresentativeSetLoader:
    """
    Parses the BGSU NR list CSV and provides chain-specific
    PDB-to-Representative lookups.
    
    Creates mapping:
      (member_pdb, member_chain) -> (rep_pdb, rep_chain, equivalence_class)
    
    Example:
      ('1S72', '0') -> ('4V9F', '0', 'NR_all_26150.5')
      ('1S72', '9') -> ('4V9F', '9', 'NR_all_25303.5')
    """
    
    def __init__(self, csv_path: str = None):
        """
        Initialize and parse the NR list CSV.
        
        Args:
            csv_path: Path to nrlist CSV file. Uses bundled file if None.
        """
        self.csv_path = csv_path or _DEFAULT_CSV_PATH
        # {(pdb_upper, chain): (rep_pdb_upper, rep_chain, eq_class_id)}
        self.mapping: Dict[Tuple[str, str], Tuple[str, str, str]] = {}
        self._parse_csv()
    
    def _parse_ife(self, ife_str: str) -> List[Tuple[str, str]]:
        """
        Parse an IFE string into (PDB, Chain) pairs.
        
        Handles:
          - Simple: "4V9F|1|0"  -> [("4V9F", "0")]
          - Multi-chain: "4V9F|1|A+4V9F|1|B" -> [("4V9F", "A"), ("4V9F", "B")]
        
        Returns:
            List of (pdb_id, chain) tuples
        """
        results = []
        # Split on '+' for multi-chain IFEs
        sub_ifes = ife_str.strip().split('+')
        for sub in sub_ifes:
            parts = sub.strip().split('|')
            if len(parts) >= 3:
                pdb_id = parts[0].strip().upper()
                chain = parts[2].strip()
                if pdb_id and chain:
                    results.append((pdb_id, chain))
        return results
    
    def _parse_csv(self):
        """Parse the NR list CSV file and build the mapping."""
        if not os.path.exists(self.csv_path):
            logger.warning(f"NR list CSV not found: {self.csv_path}")
            return
        
        count = 0
        try:
            with open(self.csv_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                for row in reader:
                    if len(row) < 3:
                        continue
                    
                    eq_class_id = row[0].strip().strip('"')
                    rep_ife_str = row[1].strip().strip('"')
                    members_str = row[2].strip().strip('"')
                    
                    # Parse representative IFE(s)
                    rep_pairs = self._parse_ife(rep_ife_str)
                    if not rep_pairs:
                        continue
                    
                    # For multi-chain representatives, use first pair as primary
                    rep_pdb, rep_chain = rep_pairs[0]
                    
                    # Parse all members
                    member_ifes = members_str.split(',')
                    for member_ife in member_ifes:
                        member_pairs = self._parse_ife(member_ife.strip())
                        for i, (mem_pdb, mem_chain) in enumerate(member_pairs):
                            # For multi-chain IFEs, map each member chain
                            # to corresponding representative chain
                            if i < len(rep_pairs):
                                target_rep_pdb, target_rep_chain = rep_pairs[i]
                            else:
                                target_rep_pdb, target_rep_chain = rep_pdb, rep_chain
                            
                            key = (mem_pdb, mem_chain)
                            self.mapping[key] = (
                                target_rep_pdb,
                                target_rep_chain,
                                eq_class_id
                            )
                            count += 1
            
            logger.info(f"Loaded {count} NR list mappings from {len(set(v[2] for v in self.mapping.values()))} equivalence classes")
            
        except Exception as e:
            logger.error(f"Error parsing NR list CSV: {e}")
    
    def get_representative(
        self, pdb_id: str, chain: str
    ) -> Optional[Tuple[str, str]]:
        """
        Get the representative PDB and chain for a given (PDB, chain).
        
        Args:
            pdb_id: PDB ID (case-insensitive)
            chain: Chain identifier
            
        Returns:
            (rep_pdb_id, rep_chain) or None if not found
        """
        key = (pdb_id.upper(), chain)
        entry = self.mapping.get(key)
        if entry:
            return (entry[0], entry[1])
        return None
    
    def get_all_representatives(
        self, pdb_id: str
    ) -> Dict[str, Tuple[str, str]]:
        """
        Get representative mappings for ALL chains of a given PDB.
        
        Args:
            pdb_id: PDB ID (case-insensitive)
            
        Returns:
            Dict mapping member_chain -> (rep_pdb, rep_chain)
            Example: {'0': ('4V9F', '0'), '9': ('4V9F', '9')}
        """
        pdb_upper = pdb_id.upper()
        result = {}
        for (mem_pdb, mem_chain), (rep_pdb, rep_chain, _) in self.mapping.items():
            if mem_pdb == pdb_upper:
                result[mem_chain] = (rep_pdb, rep_chain)
        return result
    
    def is_self_representative(self, pdb_id: str, chain: str) -> bool:
        """
        Check if a PDB chain is its own representative.
        
        Args:
            pdb_id: PDB ID
            chain: Chain identifier
            
        Returns:
            True if (pdb, chain) maps to itself as representative
        """
        rep = self.get_representative(pdb_id, chain)
        if rep is None:
            return False
        return rep[0] == pdb_id.upper() and rep[1] == chain
    
    def has_pdb(self, pdb_id: str) -> bool:
        """Check if a PDB is in the NR list at all."""
        pdb_upper = pdb_id.upper()
        return any(mem_pdb == pdb_upper for mem_pdb, _ in self.mapping)
    
    def get_equivalence_class(
        self, pdb_id: str, chain: str
    ) -> Optional[str]:
        """Get the equivalence class ID for a (PDB, chain) pair."""
        key = (pdb_id.upper(), chain)
        entry = self.mapping.get(key)
        return entry[2] if entry else None


# Singleton instance
_loader_instance: Optional[RepresentativeSetLoader] = None


def get_representative_loader(csv_path: str = None) -> RepresentativeSetLoader:
    """Get or create the singleton RepresentativeSetLoader instance."""
    global _loader_instance
    if _loader_instance is None:
        _loader_instance = RepresentativeSetLoader(csv_path)
    return _loader_instance


def reset_representative_loader():
    """Reset the singleton (useful for testing)."""
    global _loader_instance
    _loader_instance = None
