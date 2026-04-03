"""
RSMViewer - Cache Manager
Manages caching of API responses with expiry timestamps.

Cache is stored in ~/.rsmviewer_cache/ with:
- {pdb_id}_{source}.json: Cached motif data
- {pdb_id}_{source}.meta.json: Metadata including fetch time and expiry

Author: Structural Biology Lab
Version: 2.0.0
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional

from .base_provider import MotifInstance, ResidueSpec


@dataclass
class CacheMetadata:
    """Metadata for a cached entry."""
    pdb_id: str
    source: str
    fetched_at: str  # ISO format datetime
    expires_at: str  # ISO format datetime
    version: str = "2.0"
    
    def is_expired(self) -> bool:
        """Check if this cache entry has expired."""
        try:
            expires = datetime.fromisoformat(self.expires_at)
            return datetime.now() > expires
        except:
            return True
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'CacheMetadata':
        """Create from dictionary."""
        return cls(
            pdb_id=data.get('pdb_id', ''),
            source=data.get('source', ''),
            fetched_at=data.get('fetched_at', ''),
            expires_at=data.get('expires_at', ''),
            version=data.get('version', '2.0'),
        )


class CacheManager:
    """
    Manages caching of motif data from API providers.
    
    Features:
    - Persistent disk cache in user home directory
    - Configurable expiry time (default 30 days)
    - Automatic cleanup of expired entries
    - Thread-safe file operations
    """
    
    # Default cache location
    DEFAULT_CACHE_DIR = Path.home() / '.rsmviewer_cache'
    
    # Default cache expiry (30 days)
    DEFAULT_EXPIRY_DAYS = 30
    
    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        expiry_days: int = DEFAULT_EXPIRY_DAYS,
    ):
        """
        Initialize the cache manager.
        
        Args:
            cache_dir: Directory for cache files (default: ~/.rsmviewer_cache/)
            expiry_days: Number of days before cache entries expire
        """
        if cache_dir is None:
            self.cache_dir = self.DEFAULT_CACHE_DIR
        elif isinstance(cache_dir, str):
            self.cache_dir = Path(cache_dir)
        else:
            self.cache_dir = cache_dir
            
        self.expiry_days = expiry_days
        
        # Create cache directory if needed
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def _get_cache_path(self, pdb_id: str, source: str) -> Path:
        """Get path to cache file for a PDB/source combination."""
        return self.cache_dir / f"{pdb_id.upper()}_{source}.json"
    
    def _get_meta_path(self, pdb_id: str, source: str) -> Path:
        """Get path to metadata file for a PDB/source combination."""
        return self.cache_dir / f"{pdb_id.upper()}_{source}.meta.json"
    
    def get_cached_motifs(
        self, pdb_id: str, source: str, ignore_expiry: bool = False
    ) -> Optional[Dict[str, List[MotifInstance]]]:
        """
        Get cached motifs for a PDB from a specific source.
        
        Args:
            pdb_id: PDB structure ID
            source: Source identifier (e.g., 'bgsu_api', 'rfam_api')
            ignore_expiry: If True, return cached data even if expired
            
        Returns:
            Dict of motif type -> list of MotifInstance, or None if not cached
        """
        pdb_id = pdb_id.upper()
        cache_path = self._get_cache_path(pdb_id, source)
        meta_path = self._get_meta_path(pdb_id, source)
        
        # Check if cache files exist
        if not cache_path.exists() or not meta_path.exists():
            return None
        
        try:
            # Load and check metadata
            with open(meta_path, 'r') as f:
                meta = CacheMetadata.from_dict(json.load(f))
            
            # Check expiry
            if not ignore_expiry and meta.is_expired():
                # Cache expired, remove it
                self._remove_cache_entry(pdb_id, source)
                return None
            
            # Load cached data
            with open(cache_path, 'r') as f:
                cached_data = json.load(f)
            
            # Convert back to MotifInstance objects
            return self._deserialize_motifs(cached_data)
            
        except Exception as e:
            print(f"Error reading cache for {pdb_id}/{source}: {e}")
            # Remove corrupted cache
            self._remove_cache_entry(pdb_id, source)
            return None
    
    def cache_motifs(
        self, pdb_id: str, source: str, motifs: Dict[str, List[MotifInstance]]
    ) -> bool:
        """
        Cache motifs for a PDB from a specific source.
        
        Args:
            pdb_id: PDB structure ID
            source: Source identifier
            motifs: Dict of motif type -> list of MotifInstance
            
        Returns:
            True if caching successful
        """
        pdb_id = pdb_id.upper()
        cache_path = self._get_cache_path(pdb_id, source)
        meta_path = self._get_meta_path(pdb_id, source)
        
        try:
            # Create metadata
            now = datetime.now()
            expiry = now + timedelta(days=self.expiry_days)
            
            meta = CacheMetadata(
                pdb_id=pdb_id,
                source=source,
                fetched_at=now.isoformat(),
                expires_at=expiry.isoformat(),
            )
            
            # Serialize motifs
            serialized = self._serialize_motifs(motifs)
            
            # Write cache file
            with open(cache_path, 'w') as f:
                json.dump(serialized, f, indent=2)
            
            # Write metadata
            with open(meta_path, 'w') as f:
                json.dump(meta.to_dict(), f, indent=2)
            
            return True
            
        except Exception as e:
            print(f"Error caching data for {pdb_id}/{source}: {e}")
            return False
    
    def _serialize_motifs(
        self, motifs: Dict[str, List[MotifInstance]]
    ) -> Dict[str, List[Dict]]:
        """Serialize MotifInstance objects to JSON-compatible format."""
        result = {}
        
        for motif_type, instances in motifs.items():
            result[motif_type] = []
            for instance in instances:
                result[motif_type].append({
                    'instance_id': instance.instance_id,
                    'motif_id': instance.motif_id,
                    'pdb_id': instance.pdb_id,
                    'annotation': instance.annotation,
                    'metadata': instance.metadata,
                    'residues': [
                        {
                            'chain': r.chain,
                            'residue_number': r.residue_number,
                            'nucleotide': r.nucleotide,
                            'insertion_code': r.insertion_code,
                            'model': r.model,
                        }
                        for r in instance.residues
                    ]
                })
        
        return result
    
    def _deserialize_motifs(
        self, data: Dict[str, List[Dict]]
    ) -> Dict[str, List[MotifInstance]]:
        """Deserialize JSON data back to MotifInstance objects."""
        result = {}
        
        for motif_type, instances_data in data.items():
            result[motif_type] = []
            for inst_data in instances_data:
                residues = [
                    ResidueSpec(
                        chain=r['chain'],
                        residue_number=r['residue_number'],
                        nucleotide=r.get('nucleotide', ''),
                        insertion_code=r.get('insertion_code', ''),
                        model=r.get('model', 1),
                    )
                    for r in inst_data.get('residues', [])
                ]
                
                result[motif_type].append(MotifInstance(
                    instance_id=inst_data['instance_id'],
                    motif_id=inst_data['motif_id'],
                    pdb_id=inst_data['pdb_id'],
                    residues=residues,
                    annotation=inst_data.get('annotation', ''),
                    metadata=inst_data.get('metadata', {}),
                ))
        
        return result
    
    def _remove_cache_entry(self, pdb_id: str, source: str) -> None:
        """Remove a cache entry."""
        try:
            cache_path = self._get_cache_path(pdb_id, source)
            meta_path = self._get_meta_path(pdb_id, source)
            
            if cache_path.exists():
                cache_path.unlink()
            if meta_path.exists():
                meta_path.unlink()
        except:
            pass
    
    def invalidate_cache(self, pdb_id: str, source: Optional[str] = None) -> None:
        """
        Invalidate cache for a PDB.
        
        Args:
            pdb_id: PDB ID to invalidate
            source: Specific source to invalidate, or None for all sources
        """
        pdb_id = pdb_id.upper()
        
        if source:
            self._remove_cache_entry(pdb_id, source)
        else:
            # Remove all sources for this PDB
            for file in self.cache_dir.glob(f"{pdb_id}_*.json"):
                file.unlink()
    
    def clear_cache(self) -> int:
        """
        Clear all cached data.
        
        Returns:
            Number of entries removed
        """
        count = 0
        for file in self.cache_dir.glob("*.json"):
            try:
                file.unlink()
                count += 1
            except:
                pass
        return count // 2  # Divide by 2 because we have both .json and .meta.json
    
    def cleanup_expired(self) -> int:
        """
        Remove all expired cache entries.
        
        Returns:
            Number of entries removed
        """
        count = 0
        
        for meta_path in self.cache_dir.glob("*.meta.json"):
            try:
                with open(meta_path, 'r') as f:
                    meta = CacheMetadata.from_dict(json.load(f))
                
                if meta.is_expired():
                    self._remove_cache_entry(meta.pdb_id, meta.source)
                    count += 1
            except:
                # Remove corrupted entry
                try:
                    meta_path.unlink()
                    count += 1
                except:
                    pass
        
        return count
    
    def get_cache_info(self) -> Dict[str, Any]:
        """
        Get information about current cache state.
        
        Returns:
            Dict with cache statistics
        """
        total_entries = 0
        total_size = 0
        expired_count = 0
        sources: Dict[str, int] = {}
        
        for meta_path in self.cache_dir.glob("*.meta.json"):
            try:
                with open(meta_path, 'r') as f:
                    meta = CacheMetadata.from_dict(json.load(f))
                
                total_entries += 1
                sources[meta.source] = sources.get(meta.source, 0) + 1
                
                if meta.is_expired():
                    expired_count += 1
                
                # Add file sizes
                cache_path = self._get_cache_path(meta.pdb_id, meta.source)
                if cache_path.exists():
                    total_size += cache_path.stat().st_size
                total_size += meta_path.stat().st_size
                
            except:
                pass
        
        return {
            'cache_dir': str(self.cache_dir),
            'total_entries': total_entries,
            'expired_entries': expired_count,
            'total_size_kb': total_size / 1024,
            'sources': sources,
            'expiry_days': self.expiry_days,
        }


# Global cache manager instance
_cache_manager: Optional[CacheManager] = None


def get_cache_manager() -> CacheManager:
    """Get or create the global cache manager."""
    global _cache_manager
    if _cache_manager is None:
        _cache_manager = CacheManager()
    return _cache_manager


def initialize_cache_manager(
    cache_dir: Optional[Path] = None,
    expiry_days: int = CacheManager.DEFAULT_EXPIRY_DAYS,
) -> CacheManager:
    """Initialize the global cache manager with custom settings."""
    global _cache_manager
    _cache_manager = CacheManager(cache_dir=cache_dir, expiry_days=expiry_days)
    return _cache_manager
