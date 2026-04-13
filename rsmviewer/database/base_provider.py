"""
RSMViewer - Base Provider Classes
Defines abstract base classes and data structures for database providers.

This module provides:
- ResidueSpec: Standard representation for nucleotide residues
- MotifInstance: Individual motif occurrence in a PDB structure
- MotifType: Collection of motif instances of the same type
- DatabaseInfo: Metadata about a database provider
- BaseProvider: Abstract base class that all database providers must implement

All database providers (RNA 3D Atlas, Rfam, etc.) must inherit from BaseProvider
and implement the required methods to ensure consistent behavior across the plugin.

Author: Structural Biology Lab
Version: 1.0.0
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


class DatabaseSourceType(Enum):
    """Database source type."""
    LOCAL_DIRECTORY = 'local'
    API = 'api'


@dataclass
class ResidueSpec:
    """
    Standard representation for a nucleotide residue.
    
    Attributes:
        chain: Chain identifier (e.g., 'A', '1A', 'AA')
        residue_number: Residue sequence number
        nucleotide: Nucleotide type (e.g., 'A', 'U', 'G', 'C', or modified nucleotides)
        insertion_code: Optional insertion code for special residues
        model: Model number (default 1 for single-model structures)
    """
    chain: str
    residue_number: int
    nucleotide: str = ''
    insertion_code: str = ''
    model: int = 1
    
    def to_tuple(self) -> Tuple[str, int, str]:
        """Convert to legacy tuple format (nucleotide, residue_number, chain)."""
        return (self.nucleotide, self.residue_number, self.chain)
    
    @classmethod
    def from_tuple(cls, t: Tuple[str, int, str]) -> 'ResidueSpec':
        """Create from legacy tuple format."""
        return cls(nucleotide=t[0], residue_number=t[1], chain=t[2])
    
    @classmethod
    def from_atlas_spec(cls, spec: str) -> Optional['ResidueSpec']:
        """
        Parse RNA 3D Atlas residue specification.
        
        Format: PDB|Model|Chain|Nucleotide|ResNum[|insertion_code][|extra]
        Example: 4V9F|1|A|G|303
        """
        parts = str(spec).split('|')
        if len(parts) < 5:
            return None
        
        try:
            model = int(parts[1]) if parts[1].isdigit() else 1
            chain = parts[2]
            nucleotide = parts[3]
            res_num = int(parts[4])
            insertion_code = parts[5] if len(parts) > 5 and parts[5] else ''
            
            return cls(
                chain=chain,
                residue_number=res_num,
                nucleotide=nucleotide,
                insertion_code=insertion_code,
                model=model
            )
        except (ValueError, IndexError):
            return None
    
    def __hash__(self) -> int:
        return hash((self.chain, self.residue_number, self.insertion_code, self.model))
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ResidueSpec):
            return False
        return (self.chain == other.chain and 
                self.residue_number == other.residue_number and
                self.insertion_code == other.insertion_code and
                self.model == other.model)


@dataclass
class MotifInstance:
    """
    Represents a single occurrence of a motif in a PDB structure.
    
    Attributes:
        instance_id: Unique identifier for this instance
        motif_id: ID of the parent motif class/family
        pdb_id: PDB structure ID where this instance occurs
        residues: List of residues that make up this motif instance
        annotation: Optional annotation/description
        metadata: Additional metadata specific to the database
    """
    instance_id: str
    motif_id: str
    pdb_id: str
    residues: List[ResidueSpec] = field(default_factory=list)
    annotation: str = ''
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Ensure metadata is always a unique dict instance."""
        if self.metadata is None:
            object.__setattr__(self, 'metadata', {})
        elif not isinstance(self.metadata, dict):
            object.__setattr__(self, 'metadata', {})
    
    def get_chains(self) -> Set[str]:
        """Get all unique chains involved in this motif."""
        return set(r.chain for r in self.residues)
    
    def get_residue_numbers(self, chain: Optional[str] = None) -> List[int]:
        """Get residue numbers, optionally filtered by chain."""
        if chain:
            return sorted(r.residue_number for r in self.residues if r.chain == chain)
        return sorted(r.residue_number for r in self.residues)
    
    def to_legacy_format(self) -> List[Dict]:
        """
        Convert to legacy format expected by MotifSelector.
        
        Returns list of dicts: [{motif_id, chain, residues}, ...]
        """
        by_chain: Dict[str, List[int]] = {}
        for r in self.residues:
            by_chain.setdefault(r.chain, []).append(r.residue_number)
        
        result = []
        for chain, res_nums in by_chain.items():
            result.append({
                'motif_id': str(self.motif_id),
                'chain': str(chain),
                'residues': sorted(set(res_nums)),
            })
        return result


@dataclass
class MotifType:
    """
    Represents a class/family of motifs.
    
    Attributes:
        type_id: Identifier for this motif type (e.g., 'HL', 'IL', 'GNRA')
        name: Human-readable name
        description: Detailed description
        instances: List of all instances of this motif type
        source: Database source (e.g., 'RNA 3D Atlas', 'Rfam')
        metadata: Additional type-level metadata
    """
    type_id: str
    name: str = ''
    description: str = ''
    instances: List[MotifInstance] = field(default_factory=list)
    source: str = ''
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def get_instances_for_pdb(self, pdb_id: str) -> List[MotifInstance]:
        """Get all instances in a specific PDB structure."""
        pdb_id = pdb_id.upper()
        return [inst for inst in self.instances if inst.pdb_id.upper() == pdb_id]
    
    def get_all_pdb_ids(self) -> Set[str]:
        """Get all PDB IDs that have instances of this motif type."""
        return set(inst.pdb_id.upper() for inst in self.instances)


@dataclass
class DatabaseInfo:
    """
    Metadata about a database provider.
    
    Attributes:
        id: Unique identifier for this database
        name: Human-readable name
        description: Detailed description
        version: Database version
        source_type: Type of data source
        motif_types: List of available motif type IDs
        pdb_count: Number of PDB structures indexed
        last_updated: Last update timestamp
    """
    id: str
    name: str
    description: str = ''
    version: str = '1.0.0'
    source_type: DatabaseSourceType = DatabaseSourceType.LOCAL_DIRECTORY
    motif_types: List[str] = field(default_factory=list)
    pdb_count: int = 0
    last_updated: str = ''


class BaseProvider(ABC):
    """
    Abstract base class for all motif database providers.
    
    All database implementations (RNA 3D Atlas, Rfam, custom, API-based)
    must inherit from this class and implement all abstract methods.
    
    This ensures consistent behavior across the plugin regardless of
    the underlying data source or format.
    """
    
    def __init__(self, database_path: str):
        """
        Initialize the provider.
        
        Args:
            database_path: Path to the database directory or configuration file
        """
        self.database_path = Path(database_path)
        self._initialized = False
        self._pdb_index: Dict[str, List[MotifInstance]] = {}
        self._motif_types: Dict[str, MotifType] = {}
    
    @property
    @abstractmethod
    def info(self) -> DatabaseInfo:
        """
        Get database metadata.
        
        Returns:
            DatabaseInfo object with database details
        """
        pass
    
    @abstractmethod
    def initialize(self) -> bool:
        """
        Initialize the database provider.
        
        This should:
        - Load configuration
        - Build indices
        - Validate data availability
        
        Returns:
            True if initialization successful, False otherwise
        """
        pass
    
    @abstractmethod
    def get_available_motif_types(self) -> List[str]:
        """
        Get list of all available motif type IDs.
        
        Returns:
            List of motif type identifiers (e.g., ['HL', 'IL', 'J3'])
        """
        pass
    
    @abstractmethod
    def get_motif_type(self, type_id: str) -> Optional[MotifType]:
        """
        Get a specific motif type with all its instances.
        
        Args:
            type_id: Motif type identifier
            
        Returns:
            MotifType object or None if not found
        """
        pass
    
    @abstractmethod
    def get_motifs_for_pdb(self, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """
        Get all motifs for a specific PDB structure.
        
        Args:
            pdb_id: PDB structure identifier
            
        Returns:
            Dictionary mapping motif_type -> list of instances
        """
        pass
    
    @abstractmethod
    def get_available_pdb_ids(self) -> List[str]:
        """
        Get list of all PDB IDs that have motifs in this database.
        
        Returns:
            Sorted list of PDB IDs
        """
        pass
    
    @abstractmethod
    def get_motif_residues(self, pdb_id: str, motif_type: str, 
                          instance_id: str) -> List[ResidueSpec]:
        """
        Get residues for a specific motif instance.
        
        Args:
            pdb_id: PDB structure identifier
            motif_type: Motif type identifier
            instance_id: Instance identifier
            
        Returns:
            List of ResidueSpec objects
        """
        pass
    
    # Optional methods with default implementations
    
    def has_pdb(self, pdb_id: str) -> bool:
        """Check if database has motifs for a PDB structure."""
        return pdb_id.upper() in [p.upper() for p in self.get_available_pdb_ids()]
    
    def count_motifs_for_pdb(self, pdb_id: str) -> int:
        """Count total motifs for a PDB structure."""
        motifs = self.get_motifs_for_pdb(pdb_id)
        return sum(len(instances) for instances in motifs.values())
    
    def get_summary(self, pdb_id: str) -> str:
        """Get human-readable summary of motifs for a PDB."""
        motifs = self.get_motifs_for_pdb(pdb_id)
        if not motifs:
            return f"No motifs found in {pdb_id.upper()} (database: {self.info.name})"
        
        lines = [f"Motifs in {pdb_id.upper()} ({self.info.name}):"]
        for motif_type in sorted(motifs.keys()):
            lines.append(f"  {motif_type}: {len(motifs[motif_type])} instances")
        return "\n".join(lines)
    
    def refresh(self) -> bool:
        """
        Refresh database (reload files, rebuild indices).
        
        Returns:
            True if refresh successful
        """
        self._initialized = False
        self._pdb_index.clear()
        self._motif_types.clear()
        return self.initialize()
    
    @property
    def is_initialized(self) -> bool:
        """Check if provider has been initialized."""
        return self._initialized
