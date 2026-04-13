"""
RSMViewer - Configuration Module
Global settings for the plugin including data source preferences.

Author: Structural Biology Lab
Version: 1.0.0
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Dict


# Source ID mapping - centralized registry for all sources
SOURCE_ID_MAP: Dict[int, Dict] = {
    1: {
        'name': 'RNA 3D Atlas',
        'type': 'local',
        'category': 'LOCAL SOURCES',
        'subtype': 'atlas',
        'description': 'RNA 3D Motif Atlas (bundled, offline)',
        'coverage': '759 PDB structures',
        'mode': 'local',
        'command': 'rmv_db 1',
    },
    2: {
        'name': 'Rfam',
        'type': 'local',
        'category': 'LOCAL SOURCES',
        'subtype': 'rfam',
        'description': 'Rfam Database (bundled, offline)',
        'coverage': '173 PDB structures',
        'mode': 'local',
        'command': 'rmv_db 2',
    },
    3: {
        'name': 'BGSU RNA 3D Hub',
        'type': 'web',
        'category': 'ONLINE SOURCES',
        'subtype': 'bgsu',
        'description': 'BGSU RNA 3D Hub API (online, requires internet)',
        'coverage': '~3000+ PDB structures',
        'mode': 'web',
        'command': 'rmv_db 3',
    },
    4: {
        'name': 'Rfam API',
        'type': 'web',
        'category': 'ONLINE SOURCES',
        'subtype': 'rfam_api',
        'description': 'Rfam API (online, requires internet)',
        'coverage': 'All 34 Rfam named motifs',
        'mode': 'web',
        'command': 'rmv_db 4',
    },
    5: {
        'name': 'FR3D Annotations',
        'type': 'user',
        'category': 'USER ANNOTATIONS',
        'tool': 'fr3d',
        'description': 'FR3D analysis output (custom user files)',
        'coverage': 'Custom uploads',
        'mode': 'user',
        'command': 'rmv_db 5',
    },
    6: {
        'name': 'RNAMotifScan (RMS)',
        'type': 'user',
        'category': 'USER ANNOTATIONS',
        'tool': 'rms',
        'description': 'RNAMotifScan analysis output (custom user files)',
        'coverage': 'Custom uploads',
        'mode': 'user',
        'supports_filtering': True,
        'command': 'rmv_db 6',
        'command_with_filtering': 'rmv_db 6 [on|off]',
        'command_with_custom_pvalues': 'rmv_db 6 MOTIF_NAME p_value ...',
    },
    7: {
        'name': 'RNAMotifScanX (RMSX)',
        'type': 'user',
        'category': 'USER ANNOTATIONS',
        'tool': 'rmsx',
        'description': 'RNAMotifScanX analysis output (custom user files)',
        'coverage': 'Custom uploads',
        'mode': 'user',
        'supports_filtering': True,
        'command': 'rmv_db 7',
        'command_with_filtering': 'rmv_db 7 [on|off]',
        'command_with_custom_pvalues': 'rmv_db 7 MOTIF_NAME p_value ...',
    },
    8: {
        'name': 'NoBIAS',
        'type': 'user',
        'category': 'USER ANNOTATIONS',
        'tool': 'nobias',
        'description': 'NoBIAS analysis output (custom user files)',
        'coverage': 'Custom uploads',
        'mode': 'user',
        'supports_filtering': True,
        'command': 'rmv_db 8',
        'command_with_filtering': 'rmv_db 8 [on|off]',
        'command_with_custom_pvalues': 'rmv_db 8 MOTIF_NAME p_value ...',
    },
}


class SourceMode(Enum):
    """Mode for selecting data sources."""
    AUTO = "auto"       # Smart selection: local -> BGSU API -> Rfam API
    LOCAL = "local"     # Only use bundled local database
    BGSU = "bgsu"       # Prefer BGSU RNA 3D Hub API
    RFAM = "rfam"       # Prefer Rfam API
    ALL = "all"         # Combine results from all sources


class FreshnessPolicy(Enum):
    """Policy for handling potentially outdated data."""
    USE_CACHE = "use_cache"         # Always use cache if available
    CHECK_EXPIRY = "check_expiry"   # Use cache only if not expired
    FORCE_REFRESH = "force_refresh" # Always fetch fresh data


@dataclass
class CachePolicy:
    """Cache policy settings."""
    cache_days: int = 30
    policy: FreshnessPolicy = FreshnessPolicy.CHECK_EXPIRY


@dataclass
class PluginConfig:
    """
    Global configuration for the RSMViewer plugin.
    
    Attributes:
        source_mode: How to select data sources
        specific_source: Track which SPECIFIC source was explicitly selected (e.g., 'atlas', 'bgsu_api')
        source_priority: Order to try sources in AUTO mode
        freshness_policy: Cache and freshness settings
        enable_api_fallback: Whether to try API if local fails
        request_timeout: Timeout for API requests (seconds)
        verbose: Enable verbose logging
    """
    
    # Source selection
    source_mode: SourceMode = SourceMode.BGSU
    specific_source: Optional[str] = None  # Track specific source if user explicitly selected one
    source_priority: List[str] = field(default_factory=lambda: [
        "bgsu_api",   # BGSU RNA 3D Hub API (default - ~3000+ PDBs)
        "atlas",      # Local RNA 3D Atlas (bundled)
        "rfam",       # Local Rfam (bundled)
        "rfam_api",   # Rfam API
    ])
    
    # Caching settings
    freshness_policy: CachePolicy = field(default_factory=CachePolicy)
    
    # API settings
    enable_api_fallback: bool = True
    request_timeout: int = 30
    
    # Display settings
    verbose: bool = False
    
    def get_source_list(self) -> List[str]:
        """
        Get ordered list of sources to try based on current mode.
        
        Returns:
            List of source IDs in priority order
        """
        # BUG FIX: If user explicitly selected a specific source, use ONLY that source
        if self.specific_source:
            print(f"[CONFIG] get_source_list() returning ONLY specific_source: [{self.specific_source}]")
            return [self.specific_source]
        
        if self.source_mode == SourceMode.LOCAL:
            print(f"[CONFIG] get_source_list() returning LOCAL sources: ['atlas', 'rfam']")
            return ["atlas", "rfam"]
        elif self.source_mode == SourceMode.BGSU:
            print(f"[CONFIG] get_source_list() returning BGSU priority: ['bgsu_api', 'atlas']")
            return ["bgsu_api", "atlas"]
        elif self.source_mode == SourceMode.RFAM:
            print(f"[CONFIG] get_source_list() returning RFAM priority: ['rfam_api', 'rfam']")
            return ["rfam_api", "rfam"]
        elif self.source_mode == SourceMode.ALL:
            print(f"[CONFIG] get_source_list() returning ALL sources: {list(self.source_priority)}")
            return list(self.source_priority)
        else:  # AUTO
            print(f"[CONFIG] get_source_list() returning AUTO priority: {list(self.source_priority)}")
            return list(self.source_priority)
    
    def set_source_mode(self, mode: str) -> bool:
        """
        Set the source mode.
        
        Args:
            mode: Mode string ('auto', 'local', 'bgsu', 'rfam', 'all')
            
        Returns:
            True if mode was valid and set
        """
        try:
            self.source_mode = SourceMode(mode.lower())
            return True
        except ValueError:
            return False
    
    def to_dict(self) -> dict:
        """Convert config to dictionary."""
        return {
            'source_mode': self.source_mode.value,
            'source_priority': self.source_priority,
            'cache_days': self.freshness_policy.cache_days,
            'freshness_policy': self.freshness_policy.policy.value,
            'enable_api_fallback': self.enable_api_fallback,
            'request_timeout': self.request_timeout,
            'verbose': self.verbose,
        }


# Global config instance
_config: Optional[PluginConfig] = None


def get_config() -> PluginConfig:
    """Get the global plugin configuration."""
    global _config
    if _config is None:
        _config = PluginConfig()
    return _config


def set_config(config: PluginConfig) -> None:
    """Set the global plugin configuration."""
    global _config
    _config = config


def reset_config() -> PluginConfig:
    """Reset configuration to defaults."""
    global _config
    _config = PluginConfig()
    return _config
