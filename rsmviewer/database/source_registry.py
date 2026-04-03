"""
Source Registry - Manages available data sources and combining

Provides a central registry of all motif data sources (local, online, user)
with IDs for combining multiple sources in a single query.
"""

from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class SourceInfo:
    """Information about a data source"""
    source_id: int
    name: str
    provider_id: str
    source_type: str  # 'local', 'online', 'user'
    description: str
    coverage: str
    confidence: float  # 0.0 - 1.0


class SourceRegistry:
    """Central registry of all available sources"""
    
    SOURCE_MAP = {
        # Local sources (offline)
        1: SourceInfo(
            source_id=1,
            name='atlas',
            provider_id='atlas',
            source_type='local',
            description='RNA 3D Motif Atlas',
            coverage='759 PDBs',
            confidence=0.95
        ),
        2: SourceInfo(
            source_id=2,
            name='rfam',
            provider_id='rfam',
            source_type='local',
            description='Rfam Database',
            coverage='173 PDBs',
            confidence=0.90
        ),
        
        # Online sources (requires internet)
        3: SourceInfo(
            source_id=3,
            name='bgsu_api',
            provider_id='bgsu_api',
            source_type='online',
            description='BGSU RNA 3D Hub',
            coverage='~3000+ PDBs',
            confidence=0.85
        ),
        4: SourceInfo(
            source_id=4,
            name='rfam_api',
            provider_id='rfam_api',
            source_type='online',
            description='Rfam API',
            coverage='All Rfam',
            confidence=0.85
        ),
        
        # User annotations
        5: SourceInfo(
            source_id=5,
            name='fr3d',
            provider_id='user',
            source_type='user',
            description='FR3D Annotations',
            coverage='Custom',
            confidence=0.80
        ),
        6: SourceInfo(
            source_id=6,
            name='rnamotifscan',
            provider_id='user',
            source_type='user',
            description='RNAMotifScan',
            coverage='Custom',
            confidence=0.75
        ),
    }
    
    def __init__(self):
        self.sources = self.SOURCE_MAP.copy()
    
    def get_source(self, source_id: int) -> Optional[SourceInfo]:
        """Get source by ID"""
        return self.sources.get(source_id)
    
    def get_all_sources(self) -> Dict[int, SourceInfo]:
        """Get all sources"""
        return self.sources
    
    def get_source_by_name(self, name: str) -> Optional[SourceInfo]:
        """Get source by name"""
        for source in self.sources.values():
            if source.name == name.lower():
                return source
        return None
    
    def validate_source_ids(self, source_ids: List[int]) -> Tuple[bool, str]:
        """Validate that source IDs are valid"""
        invalid = [sid for sid in source_ids if sid not in self.sources]
        if invalid:
            return False, f"Invalid source IDs: {invalid}"
        if len(source_ids) != len(set(source_ids)):
            return False, "Duplicate source IDs not allowed"
        return True, "OK"
    
    def get_source_names(self, source_ids: List[int]) -> List[str]:
        """Get source names for IDs"""
        names = []
        for sid in source_ids:
            source = self.get_source(sid)
            if source:
                names.append(source.name)
        return names
    
    def get_source_descriptions(self, source_ids: List[int]) -> List[str]:
        """Get source descriptions for IDs"""
        descs = []
        for sid in source_ids:
            source = self.get_source(sid)
            if source:
                descs.append(f"{source.description} ({source.coverage})")
        return descs


# Global registry instance
_source_registry = SourceRegistry()


def get_source_registry() -> SourceRegistry:
    """Get global source registry"""
    return _source_registry
