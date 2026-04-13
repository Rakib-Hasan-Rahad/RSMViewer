"""
RSMViewer - RNA 3D Atlas Database Provider
Implements BaseProvider for RNA 3D Motif Atlas JSON files.

This provider loads motif data from the RNA 3D Motif Atlas format,
which contains hairpin loops (HL), internal loops (IL), and junctions (J3-J7).

Author: CBB Lab
Version: 1.0.0
"""

from __future__ import annotations

import json
import os
import re
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
from .converters import AtlasJSONConverter


class RNA3DAtlasProvider(BaseProvider):
    """
    Database provider for RNA 3D Motif Atlas.
    
    Supports loading motifs from JSON files in the Atlas format.
    Automatically discovers and loads the latest version of each motif type.
    """
    
    # Default motif types to look for
    DEFAULT_MOTIF_TYPES = ['HL', 'IL', 'J3', 'J4', 'J5', 'J6', 'J7']
    
    def __init__(self, database_path: str):
        """
        Initialize the RNA 3D Atlas provider.
        
        Args:
            database_path: Path to directory containing Atlas JSON files
                          (e.g., 'motif_database/RNA 3D motif atlas/')
        """
        super().__init__(database_path)
        self._converter = AtlasJSONConverter()
        self._registry: Dict[str, Any] = {}
        self._resolved_files: Dict[str, Path] = {}
        self._version = '4.5'
        
    @property
    def info(self) -> DatabaseInfo:
        """Get database metadata."""
        return DatabaseInfo(
            id='rna_3d_atlas',
            name='RNA 3D Motif Atlas',
            description='Structural motifs from the RNA 3D Motif Atlas (BGSU)',
            version=self._version,
            source_type=DatabaseSourceType.LOCAL_DIRECTORY,
            motif_types=list(self._motif_types.keys()),
            pdb_count=len(self._pdb_index),
            last_updated=self._registry.get('last_updated', '')
        )
    
    def initialize(self) -> bool:
        """
        Initialize the provider by loading registry and building index.
        
        Returns:
            True if initialization successful
        """
        try:
            # Load registry if exists
            self._load_registry()
            
            # Resolve motif files (find latest versions)
            self._resolved_files = self._resolve_motif_files()
            
            if not self._resolved_files:
                print(f"Warning: No motif files found in {self.database_path}")
                return False
            
            # Load all motif types
            for motif_type, file_path in self._resolved_files.items():
                self._load_motif_type(motif_type, file_path)
            
            # Build PDB index
            self._build_pdb_index()
            
            self._initialized = True
            return True
            
        except Exception as e:
            print(f"Error initializing RNA 3D Atlas provider: {e}")
            return False
    
    def _load_registry(self) -> None:
        """Load motif registry file if it exists."""
        registry_path = self.database_path / 'motif_registry.json'
        if registry_path.exists():
            try:
                with open(registry_path, 'r') as f:
                    self._registry = json.load(f)
            except Exception:
                self._registry = {}
    
    def _resolve_motif_files(self) -> Dict[str, Path]:
        """
        Resolve motif type to JSON file paths.
        
        Supports:
        - Explicit registry entries
        - Auto-discovery of versioned files (e.g., hl_4.5.json)
        
        Returns:
            Dict mapping motif type to file path
        """
        override_version = os.environ.get('RNA_MOTIF_ATLAS_VERSION', '').strip() or None
        
        # Get motif types from registry or use defaults
        registry_files = self._registry.get('motif_files', {})
        motif_types = list(registry_files.keys()) if registry_files else self.DEFAULT_MOTIF_TYPES
        
        resolved: Dict[str, Path] = {}
        
        for motif_type in motif_types:
            # Try explicit registry path first
            cfg = registry_files.get(motif_type, {})
            if isinstance(cfg, dict) and cfg.get('file'):
                explicit_path = (self.database_path / cfg['file']).resolve()
                if explicit_path.exists():
                    resolved[motif_type] = explicit_path
                    continue
            
            # Auto-discover versioned file
            discovered = self._discover_latest_version(motif_type, override_version)
            if discovered:
                resolved[motif_type] = discovered
        
        return resolved
    
    def _discover_latest_version(self, motif_type: str, 
                                  override_version: Optional[str]) -> Optional[Path]:
        """
        Find the latest version of a motif type file.
        
        Pattern: <prefix>_<version>.json (e.g., hl_4.5.json)
        """
        prefix = motif_type.lower()
        candidates = list(self.database_path.glob(f'{prefix}_*.json'))
        
        if not candidates:
            return None
        
        versioned: List[tuple] = []
        override_match: Optional[Path] = None
        
        for path in candidates:
            stem = path.stem
            match = re.match(rf'^{re.escape(prefix)}_(.+)$', stem)
            if not match:
                continue
            
            version_str = match.group(1)
            
            if override_version and version_str == override_version:
                override_match = path
            
            parsed = self._parse_version(version_str)
            if parsed:
                versioned.append((parsed, path))
        
        if override_match:
            return override_match.resolve()
        
        if not versioned:
            # Fall back to newest by filename
            return sorted(candidates, key=lambda p: p.name)[-1].resolve()
        
        versioned.sort(key=lambda vp: vp[0])
        latest_path = versioned[-1][1]
        
        # Extract version for metadata
        match = re.match(rf'^{re.escape(prefix)}_(.+)$', latest_path.stem)
        if match:
            self._version = match.group(1)
        
        return latest_path.resolve()
    
    @staticmethod
    def _parse_version(version: str) -> Optional[tuple]:
        """Parse version string to tuple of ints."""
        parts = [p for p in version.split('.') if p]
        if not parts:
            return None
        
        try:
            return tuple(int(p) for p in parts)
        except ValueError:
            return None
    
    def _load_motif_type(self, motif_type: str, file_path: Path) -> None:
        """Load a motif type from JSON file."""
        motif_types = self._converter.convert_file(file_path)
        
        for mt in motif_types:
            # Ensure type ID is uppercase and consistent
            mt.type_id = motif_type.upper()
            self._motif_types[mt.type_id] = mt
    
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
        return self._motif_types.get(type_id.upper())
    
    def get_motifs_for_pdb(self, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """Get all motifs for a PDB structure, grouped by type."""
        pdb_id = pdb_id.upper()
        instances = self._pdb_index.get(pdb_id, [])
        
        grouped: Dict[str, List[MotifInstance]] = {}
        for inst in instances:
            motif_type = self._get_instance_type(inst)
            if motif_type not in grouped:
                grouped[motif_type] = []
            grouped[motif_type].append(inst)
        
        return grouped
    
    def _get_instance_type(self, instance: MotifInstance) -> str:
        """Get motif type for an instance."""
        # Extract type from motif_id (e.g., HL_00317.1 -> HL)
        motif_id = instance.motif_id
        if '_' in motif_id:
            return motif_id.split('_')[0].upper()
        return 'UNKNOWN'
    
    def get_available_pdb_ids(self) -> List[str]:
        """Get list of all PDB IDs with motifs."""
        return sorted(self._pdb_index.keys())
    
    def get_motif_residues(self, pdb_id: str, motif_type: str,
                          instance_id: str) -> List[ResidueSpec]:
        """Get residues for a specific motif instance."""
        pdb_id = pdb_id.upper()
        motif_type = motif_type.upper()
        
        instances = self._pdb_index.get(pdb_id, [])
        
        for inst in instances:
            if inst.instance_id == instance_id:
                return inst.residues
        
        return []
    
    def get_instances_for_pdb(self, pdb_id: str, motif_type: str) -> List[MotifInstance]:
        """Get all instances of a motif type in a PDB."""
        motifs = self.get_motifs_for_pdb(pdb_id)
        return motifs.get(motif_type.upper(), [])
