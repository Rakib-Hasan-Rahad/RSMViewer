"""
RSMViewer - Database Registry
Central registry for managing multiple motif database providers.

The registry allows:
- Registering multiple database providers
- Switching between active databases
- Combining results from multiple databases
- Automatic discovery and initialization of providers

Author: Structural Biology Lab
Version: 2.0.0
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Type

from .base_provider import BaseProvider, DatabaseInfo, MotifInstance


class DatabaseRegistry:
    """
    Central registry for managing motif database providers.
    
    Supports registering multiple providers and switching between them.
    Can also combine results from multiple databases for comprehensive analysis.
    """
    
    def __init__(self):
        """Initialize the registry."""
        self._providers: Dict[str, BaseProvider] = {}
        self._active_provider_id: Optional[str] = None
        self._initialized = False
    
    def register_provider(self, provider: BaseProvider, 
                         provider_id: Optional[str] = None) -> bool:
        """
        Register a database provider.
        
        Args:
            provider: Provider instance to register
            provider_id: Optional custom ID (defaults to provider.info.id)
            
        Returns:
            True if registration successful
        """
        try:
            # Initialize provider if needed
            if not provider.is_initialized:
                if not provider.initialize():
                    print(f"Warning: Failed to initialize provider {provider.info.id}")
                    return False
            
            pid = provider_id or provider.info.id
            self._providers[pid] = provider
            
            # Set as active if first provider
            if self._active_provider_id is None:
                self._active_provider_id = pid
            
            return True
            
        except Exception as e:
            print(f"Error registering provider: {e}")
            return False
    
    def unregister_provider(self, provider_id: str) -> bool:
        """
        Unregister a database provider.
        
        Args:
            provider_id: ID of provider to unregister
            
        Returns:
            True if successful
        """
        if provider_id not in self._providers:
            return False
        
        del self._providers[provider_id]
        
        # Clear active if it was the removed provider
        if self._active_provider_id == provider_id:
            self._active_provider_id = next(iter(self._providers.keys()), None)
        
        return True
    
    def set_active_provider(self, provider_id: str) -> bool:
        """
        Set the active database provider.
        
        Args:
            provider_id: ID of provider to activate
            
        Returns:
            True if provider exists and was activated
        """
        if provider_id not in self._providers:
            print(f"Provider '{provider_id}' not found. Available: {list(self._providers.keys())}")
            return False
        
        self._active_provider_id = provider_id
        return True
    
    def get_active_provider(self) -> Optional[BaseProvider]:
        """Get the currently active provider."""
        if self._active_provider_id:
            return self._providers.get(self._active_provider_id)
        return None
    
    def get_provider(self, provider_id: str) -> Optional[BaseProvider]:
        """Get a specific provider by ID."""
        return self._providers.get(provider_id)
    
    def get_all_providers(self) -> Dict[str, BaseProvider]:
        """Get all registered providers."""
        return dict(self._providers)
    
    def get_provider_ids(self) -> List[str]:
        """Get list of all registered provider IDs."""
        return list(self._providers.keys())
    
    def get_database_infos(self) -> List[DatabaseInfo]:
        """Get info for all registered databases."""
        return [p.info for p in self._providers.values()]
    
    # Convenience methods that delegate to active provider
    
    def get_available_motif_types(self) -> List[str]:
        """Get motif types from active provider."""
        provider = self.get_active_provider()
        if provider:
            return provider.get_available_motif_types()
        return []
    
    def get_motifs_for_pdb(self, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """Get motifs from active provider."""
        provider = self.get_active_provider()
        if provider:
            return provider.get_motifs_for_pdb(pdb_id)
        return {}
    
    def get_available_pdb_ids(self) -> List[str]:
        """Get PDB IDs from active provider."""
        provider = self.get_active_provider()
        if provider:
            return provider.get_available_pdb_ids()
        return []
    
    def has_pdb(self, pdb_id: str) -> bool:
        """Check if active provider has motifs for a PDB."""
        provider = self.get_active_provider()
        if provider:
            return provider.has_pdb(pdb_id)
        return False
    
    # Multi-provider methods
    
    def get_all_motifs_for_pdb(self, pdb_id: str) -> Dict[str, Dict[str, List[MotifInstance]]]:
        """
        Get motifs from ALL providers for a PDB structure.
        
        Returns:
            Dict mapping provider_id -> {motif_type -> [instances]}
        """
        results: Dict[str, Dict[str, List[MotifInstance]]] = {}
        
        for pid, provider in self._providers.items():
            motifs = provider.get_motifs_for_pdb(pdb_id)
            if motifs:
                results[pid] = motifs
        
        return results
    
    def get_combined_motifs_for_pdb(self, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """
        Get combined motifs from all providers, with prefixed type IDs.
        
        Motif types are prefixed with provider ID to avoid collisions.
        e.g., 'atlas:HL', 'rfam:GNRA'
        
        Returns:
            Dict mapping prefixed_motif_type -> [instances]
        """
        combined: Dict[str, List[MotifInstance]] = {}
        
        for pid, provider in self._providers.items():
            motifs = provider.get_motifs_for_pdb(pdb_id)
            for motif_type, instances in motifs.items():
                prefixed_type = f"{pid}:{motif_type}"
                combined[prefixed_type] = instances
        
        return combined
    
    def search_all_providers(self, pdb_id: str) -> Dict[str, bool]:
        """
        Check which providers have motifs for a PDB.
        
        Returns:
            Dict mapping provider_id -> has_motifs
        """
        return {
            pid: provider.has_pdb(pdb_id)
            for pid, provider in self._providers.items()
        }
    
    def get_summary(self) -> str:
        """Get summary of all registered databases."""
        lines = ["Registered Databases:"]
        
        for pid, provider in self._providers.items():
            info = provider.info
            active = " (ACTIVE)" if pid == self._active_provider_id else ""
            lines.append(f"  [{pid}]{active}")
            lines.append(f"    Name: {info.name}")
            lines.append(f"    Motif types: {len(info.motif_types)}")
            lines.append(f"    PDB structures: {info.pdb_count}")
        
        return "\n".join(lines)
    
    def print_summary(self) -> None:
        """Print summary to console."""
        print(self.get_summary())


# Global registry instance
_registry_instance: Optional[DatabaseRegistry] = None


def get_registry() -> DatabaseRegistry:
    """Get or create global registry instance."""
    global _registry_instance
    if _registry_instance is None:
        _registry_instance = DatabaseRegistry()
    return _registry_instance


def initialize_registry(motif_database_path: str, enable_api: bool = True) -> DatabaseRegistry:
    """
    Initialize registry with default providers (local + API).
    
    Args:
        motif_database_path: Path to motif_database directory
        enable_api: Whether to enable API providers for online data
        
    Returns:
        Initialized registry
    """
    from .atlas_provider import RNA3DAtlasProvider
    from .rfam_provider import RfamProvider
    from .cache_manager import get_cache_manager
    
    registry = get_registry()
    db_path = Path(motif_database_path)
    cache_manager = get_cache_manager()
    
    # ========================================
    # LOCAL PROVIDERS (bundled data, fast)
    # ========================================
    
    # Try to register RNA 3D Atlas provider (local)
    atlas_path = db_path / 'RNA 3D motif atlas'
    if atlas_path.exists():
        atlas_provider = RNA3DAtlasProvider(str(atlas_path))
        if registry.register_provider(atlas_provider, 'atlas'):
            print(f"Registered RNA 3D Atlas database (local)")
    
    # Try to register Rfam provider (local)
    rfam_path = db_path / 'Rfam motif database'
    if rfam_path.exists():
        rfam_provider = RfamProvider(str(rfam_path))
        if registry.register_provider(rfam_provider, 'rfam'):
            print(f"Registered Rfam motif database (local)")
    
    # ========================================
    # API PROVIDERS (online, comprehensive)
    # ========================================
    
    if enable_api:
        try:
            from .bgsu_api_provider import BGSUAPIProvider
            from .rfam_api_provider import RfamAPIProvider
            
            # Register BGSU RNA 3D Hub API provider
            bgsu_api = BGSUAPIProvider(cache_manager=cache_manager)
            if registry.register_provider(bgsu_api, 'bgsu_api'):
                print(f"Registered BGSU RNA 3D Hub API (~3000+ PDBs)")
            
            # Register Rfam API provider
            rfam_api = RfamAPIProvider(cache_manager=cache_manager)
            if registry.register_provider(rfam_api, 'rfam_api'):
                print(f"Registered Rfam API (named motifs)")
                
        except Exception as e:
            print(f"Note: API providers not available ({e})")
    
    # Initialize source selector with all providers
    try:
        from .source_selector import initialize_source_selector
        initialize_source_selector(registry.get_all_providers(), cache_manager)
    except Exception as e:
        print(f"Note: Source selector not initialized ({e})")
    
    return registry
