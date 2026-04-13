"""
RSMViewer - Source Selector
Intelligently selects and combines motif data from multiple sources.

The source selector implements a fallback chain:
1. Local database (bundled, fast)
2. BGSU RNA 3D Hub API (comprehensive, ~3000+ PDBs)
3. Rfam API (named motifs)

It handles:
- Automatic source selection based on availability
- Combining results from multiple sources
- Cache management for API responses
- Graceful fallback when sources fail

Author: Structural Biology Lab
Version: 1.0.0
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

from .base_provider import BaseProvider, MotifInstance
from .cache_manager import CacheManager, get_cache_manager
from .config import get_config, FreshnessPolicy, SourceMode


class SourceSelector:
    """
    Selects and combines motif data from multiple sources.
    
    Implements smart fallback logic:
    - Tries local sources first (faster, offline)
    - Falls back to API sources if local doesn't have data
    - Combines results when configured to use all sources
    - Caches API responses for future use
    """
    
    def __init__(
        self,
        providers: Dict[str, BaseProvider],
        cache_manager: Optional[CacheManager] = None,
    ):
        """
        Initialize the source selector.
        
        Args:
            providers: Dict mapping source IDs to provider instances
            cache_manager: Optional cache manager for API caching
        """
        self.providers = providers
        self.cache_manager = cache_manager or get_cache_manager()
        self._last_source_used: Optional[str] = None
    
    def get_motifs_for_pdb(
        self,
        pdb_id: str,
        source_override: Optional[str] = None,
        force_refresh: bool = False,
    ) -> Tuple[Dict[str, List[MotifInstance]], str]:
        """
        Get motifs for a PDB from the best available source.
        
        Args:
            pdb_id: PDB structure ID
            source_override: Force use of specific source
            force_refresh: Force re-fetch from API (ignore cache)
            
        Returns:
            Tuple of (motifs dict, source used)
        """
        pdb_id = pdb_id.strip().upper()
        config = get_config()
        
        # If specific source requested, use only that
        if source_override and source_override in self.providers:
            provider = self.providers[source_override]
            if force_refresh and hasattr(provider, 'cache_manager'):
                self.cache_manager.invalidate_cache(pdb_id, source_override)
            
            motifs = provider.get_motifs_for_pdb(pdb_id)
            self._last_source_used = source_override
            print(f"[SOURCE_SELECTOR] Used source_override: {source_override}")
            return motifs, source_override
        
        # Get ordered list of sources to try
        sources_to_try = config.get_source_list()
        print(f"[SOURCE_SELECTOR] No source_override, trying sources in order: {sources_to_try}")
        print(f"[SOURCE_SELECTOR] config.specific_source = {config.specific_source}")
        
        # If using ALL mode, combine results
        if config.source_mode == SourceMode.ALL:
            return self._get_from_all_sources(pdb_id, sources_to_try)
        
        # Try sources in order until one succeeds
        for source_id in sources_to_try:
            if source_id not in self.providers:
                continue
            
            provider = self.providers[source_id]
            
            try:
                # Handle force refresh
                if force_refresh and source_id.endswith('_api'):
                    self.cache_manager.invalidate_cache(pdb_id, source_id)
                
                motifs = provider.get_motifs_for_pdb(pdb_id)
                
                if motifs:
                    self._last_source_used = source_id
                    return motifs, source_id
                    
            except Exception as e:
                print(f"Error getting motifs from {source_id}: {e}")
                continue
        
        # No source had data
        self._last_source_used = None
        return {}, ""
    
    def _get_from_all_sources(
        self, pdb_id: str, source_ids: List[str]
    ) -> Tuple[Dict[str, List[MotifInstance]], str]:
        """
        Get and combine motifs from all available sources.
        
        Args:
            pdb_id: PDB ID
            source_ids: List of source IDs to try
            
        Returns:
            Tuple of (combined motifs dict, comma-separated sources used)
        """
        combined: Dict[str, List[MotifInstance]] = {}
        sources_used: List[str] = []
        
        for source_id in source_ids:
            if source_id not in self.providers:
                continue
            
            try:
                provider = self.providers[source_id]
                motifs = provider.get_motifs_for_pdb(pdb_id)
                
                if motifs:
                    sources_used.append(source_id)
                    
                    for motif_type, instances in motifs.items():
                        # Add source prefix to avoid ID collisions
                        prefixed_type = f"{source_id}:{motif_type}"
                        combined[prefixed_type] = instances
                        
            except Exception as e:
                print(f"Error getting motifs from {source_id}: {e}")
                continue
        
        self._last_source_used = ",".join(sources_used) if sources_used else None
        return combined, self._last_source_used or ""
    
    def get_available_sources(self) -> List[str]:
        """Get list of available source IDs."""
        return list(self.providers.keys())
    
    def get_last_source_used(self) -> Optional[str]:
        """Get the source used for the last query."""
        return self._last_source_used
    
    def check_pdb_availability(self, pdb_id: str) -> Dict[str, bool]:
        """
        Check which sources have data for a PDB.
        
        Args:
            pdb_id: PDB ID to check
            
        Returns:
            Dict mapping source IDs to availability
        """
        pdb_id = pdb_id.upper()
        availability = {}
        
        for source_id, provider in self.providers.items():
            try:
                # For local providers, check if PDB is in their list
                if hasattr(provider, 'get_available_pdb_ids'):
                    availability[source_id] = pdb_id in [
                        p.upper() for p in provider.get_available_pdb_ids()
                    ]
                else:
                    # For API providers, we'd need to actually query
                    # For now, mark as "unknown" by not including
                    pass
            except:
                availability[source_id] = False
        
        return availability
    
    def get_source_info(self) -> Dict[str, Dict]:
        """
        Get information about all registered sources.
        
        Returns:
            Dict mapping source IDs to info dicts
        """
        info = {}
        
        for source_id, provider in self.providers.items():
            info[source_id] = {
                'name': provider.info.name if hasattr(provider, 'info') else source_id,
                'type': 'api' if source_id.endswith('_api') else 'local',
                'motif_types': len(provider.get_available_motif_types()),
            }
            
            # For local providers, include PDB count
            if not source_id.endswith('_api'):
                info[source_id]['pdb_count'] = len(provider.get_available_pdb_ids())
        
        return info
    
    def refresh_from_api(self, pdb_id: str) -> Tuple[Dict[str, List[MotifInstance]], str]:
        """
        Force refresh motifs from API sources.
        
        Args:
            pdb_id: PDB ID to refresh
            
        Returns:
            Tuple of (motifs dict, source used)
        """
        # Invalidate cache for API sources
        pdb_id = pdb_id.upper()
        self.cache_manager.invalidate_cache(pdb_id, "bgsu_api")
        self.cache_manager.invalidate_cache(pdb_id, "rfam_api")
        
        # Try API sources
        for source_id in ["bgsu_api", "rfam_api"]:
            if source_id in self.providers:
                try:
                    motifs = self.providers[source_id].get_motifs_for_pdb(pdb_id)
                    if motifs:
                        self._last_source_used = source_id
                        return motifs, source_id
                except Exception as e:
                    print(f"Error refreshing from {source_id}: {e}")
                    continue
        
        return {}, ""
    
    def get_motifs_for_pdb_and_tool(
        self, pdb_id: str, tool_name: str
    ) -> Dict[str, List[MotifInstance]]:
        """
        Get motifs for a PDB from a specific tool (user annotations).
        
        Args:
            pdb_id: PDB structure ID
            tool_name: Tool name (e.g., 'fr3d', 'rnamotifscan', 'rnamotifscanx')
            
        Returns:
            Dict of motif type -> list of MotifInstance objects
        """
        pdb_id = pdb_id.strip().upper()
        tool_name = tool_name.strip().lower()
        
        # Get user annotations provider
        if 'user' not in self.providers:
            return {}
        
        provider = self.providers['user']
        
        # Set active tool filter before loading
        if hasattr(provider, 'set_active_tool'):
            provider.set_active_tool(tool_name)
        
        # Get motifs for this PDB (will be filtered by active_tool)
        all_motifs = provider.get_motifs_for_pdb(pdb_id)
        if hasattr(provider, 'set_active_tool'):
            provider.set_active_tool(None)
        
        return all_motifs




# Global source selector instance
_source_selector: Optional[SourceSelector] = None


def get_source_selector() -> Optional[SourceSelector]:
    """Get the global source selector (if initialized)."""
    return _source_selector


def initialize_source_selector(
    providers: Dict[str, BaseProvider],
    cache_manager: Optional[CacheManager] = None,
) -> SourceSelector:
    """Initialize the global source selector."""
    global _source_selector
    _source_selector = SourceSelector(providers, cache_manager)
    return _source_selector
