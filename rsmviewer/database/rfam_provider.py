"""
RSMViewer - Rfam Database Provider
Implements BaseProvider for Rfam motif database (Stockholm SEED files).

This provider loads motif data from Rfam SEED files in Stockholm format.
Each motif type (GNRA, T-loop, etc.) is in its own subdirectory.

Author: Structural Biology Lab
Version: 2.0.0
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from .base_provider import (
    BaseProvider,
    DatabaseInfo,
    DatabaseSourceType,
    MotifInstance,
    MotifType,
    ResidueSpec,
)
from .converters import StockholmConverter


class RfamProvider(BaseProvider):
    """
    Database provider for Rfam motif database.
    
    Directory structure expected:
    motif_database/Rfam motif database/
    ├── GNRA/
    │   ├── SEED
    │   ├── CM
    │   └── SEED.png
    ├── T-loop/
    │   ├── SEED
    │   └── ...
    └── ...
    
    Each subdirectory contains a SEED file in Stockholm format.
    """
    
    def __init__(self, database_path: str):
        """
        Initialize the Rfam provider.
        
        Args:
            database_path: Path to Rfam motif database directory
                          (e.g., 'motif_database/Rfam motif database/')
        """
        super().__init__(database_path)
        self._converter = StockholmConverter()
        self._motif_dirs: Dict[str, Path] = {}
        
    @property
    def info(self) -> DatabaseInfo:
        """Get database metadata."""
        return DatabaseInfo(
            id='rfam',
            name='Rfam Motif Database',
            description='Curated RNA structural motifs from Rfam database',
            version='1.0.0',
            source_type=DatabaseSourceType.LOCAL_DIRECTORY,
            motif_types=list(self._motif_types.keys()),
            pdb_count=len(self._pdb_index),
            last_updated=''
        )
    
    def initialize(self) -> bool:
        """
        Initialize the provider by discovering and loading motif directories.
        
        Returns:
            True if initialization successful
        """
        try:
            # Discover motif directories
            self._motif_dirs = self._discover_motif_directories()
            
            if not self._motif_dirs:
                print(f"Warning: No motif directories found in {self.database_path}")
                return False
            
            # Load all motif types
            for motif_name, dir_path in self._motif_dirs.items():
                self._load_motif_directory(motif_name, dir_path)
            
            # Build PDB index
            self._build_pdb_index()
            
            self._initialized = True
            return True
            
        except Exception as e:
            print(f"Error initializing Rfam provider: {e}")
            return False
    
    def _discover_motif_directories(self) -> Dict[str, Path]:
        """
        Discover subdirectories containing SEED files.
        
        Returns:
            Dict mapping motif name to directory path
        """
        discovered: Dict[str, Path] = {}
        
        if not self.database_path.exists():
            return discovered
        
        for entry in self.database_path.iterdir():
            if not entry.is_dir():
                continue
            
            # Skip hidden directories and common non-motif folders
            if entry.name.startswith('.') or entry.name.startswith('__'):
                continue
            
            # Check for SEED file
            seed_file = entry / 'SEED'
            if seed_file.exists():
                discovered[entry.name] = entry
        
        return discovered
    
    def _load_motif_directory(self, motif_name: str, dir_path: Path) -> None:
        """
        Load motif data from a directory's SEED file.
        
        Args:
            motif_name: Name of the motif (directory name)
            dir_path: Path to the motif directory
        """
        seed_file = dir_path / 'SEED'
        
        if not seed_file.exists():
            return
        
        try:
            motif_types = self._converter.convert_file(seed_file)
            
            for mt in motif_types:
                # Use original directory name as display name
                mt.name = motif_name
                # Normalize type_id
                type_id = self._normalize_type_id(motif_name)
                mt.type_id = type_id
                
                self._motif_types[type_id] = mt
                
        except Exception as e:
            print(f"Error loading motif {motif_name}: {e}")
    
    def _normalize_type_id(self, name: str) -> str:
        """
        Normalize motif name to consistent type ID.
        
        Examples:
            'GNRA' -> 'GNRA'
            'T-loop' -> 'T_LOOP'
            'k-turn-1' -> 'K_TURN_1'
        """
        # Replace dashes and spaces with underscores, uppercase
        normalized = name.replace('-', '_').replace(' ', '_')
        return normalized.upper()
    
    def _build_pdb_index(self) -> None:
        """Build index of PDB -> motif instances."""
        self._pdb_index.clear()
        
        for motif_type in self._motif_types.values():
            for instance in motif_type.instances:
                pdb_id = instance.pdb_id.upper()
                if pdb_id not in self._pdb_index:
                    self._pdb_index[pdb_id] = []
                self._pdb_index[pdb_id].append(instance)
    
    def get_available_motif_types(self) -> List[str]:
        """Get list of available motif type IDs."""
        return sorted(self._motif_types.keys())
    
    def get_motif_type(self, type_id: str) -> Optional[MotifType]:
        """Get a specific motif type."""
        # Try exact match first
        type_id_upper = type_id.upper()
        if type_id_upper in self._motif_types:
            return self._motif_types[type_id_upper]
        
        # Try normalized lookup
        normalized = self._normalize_type_id(type_id)
        return self._motif_types.get(normalized)
    
    def get_motifs_for_pdb(self, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """Get all motifs for a PDB structure, grouped by type."""
        pdb_id = pdb_id.upper()
        instances = self._pdb_index.get(pdb_id, [])
        
        grouped: Dict[str, List[MotifInstance]] = {}
        for inst in instances:
            type_id = self._normalize_type_id(inst.motif_id)
            if type_id not in grouped:
                grouped[type_id] = []
            grouped[type_id].append(inst)
        
        return grouped
    
    def get_available_pdb_ids(self) -> List[str]:
        """Get list of all PDB IDs with motifs."""
        return sorted(self._pdb_index.keys())
    
    def get_motif_residues(self, pdb_id: str, motif_type: str,
                          instance_id: str) -> List[ResidueSpec]:
        """Get residues for a specific motif instance."""
        pdb_id = pdb_id.upper()
        
        instances = self._pdb_index.get(pdb_id, [])
        
        for inst in instances:
            if inst.instance_id == instance_id:
                return inst.residues
        
        return []
    
    def get_motif_names(self) -> Dict[str, str]:
        """
        Get mapping of type IDs to original names.
        
        Useful for display purposes.
        
        Returns:
            Dict mapping type_id -> original_name
        """
        return {
            type_id: mt.name 
            for type_id, mt in self._motif_types.items()
        }
    
    def get_motif_description(self, type_id: str) -> str:
        """Get description for a motif type."""
        mt = self.get_motif_type(type_id)
        if mt:
            return mt.description
        return ''
    
    def get_instances_for_pdb(self, pdb_id: str, motif_type: str) -> List[MotifInstance]:
        """Get all instances of a motif type in a PDB."""
        motifs = self.get_motifs_for_pdb(pdb_id)
        return motifs.get(self._normalize_type_id(motif_type), [])
