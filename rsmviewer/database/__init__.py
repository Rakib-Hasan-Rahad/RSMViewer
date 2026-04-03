"""
RSMViewer - Database Package
Provides scalable database provider abstraction for multiple motif databases.

This package implements a plugin architecture for motif databases:
- BaseProvider: Abstract base class for all database providers
- DatabaseRegistry: Central registry for managing database providers
- Format converters for different data formats
- API providers for online data sources
- Cache management for API responses

Supported databases:
- RNA 3D Motif Atlas (JSON format, local)
- Rfam Motif Database (Stockholm format, local)
- BGSU RNA 3D Hub API (online, ~3000+ PDBs)
- Rfam API (online, named motifs)

Author: CBB LAB KU @Rakib Hasan Rahad
Version: 2.1.0
"""

from .base_provider import (
    BaseProvider,
    MotifInstance,
    MotifType,
    ResidueSpec,
    DatabaseInfo,
)
from .registry import DatabaseRegistry, get_registry, initialize_registry
from .atlas_provider import RNA3DAtlasProvider
from .rfam_provider import RfamProvider

# API providers and utilities
from .cache_manager import CacheManager, get_cache_manager, initialize_cache_manager
from .config import PluginConfig, SourceMode, FreshnessPolicy, CachePolicy, get_config, set_config
from .source_selector import SourceSelector, get_source_selector, initialize_source_selector

# API providers (optional, may fail if network unavailable)
try:
    from .bgsu_api_provider import BGSUAPIProvider
    from .rfam_api_provider import RfamAPIProvider
    _API_AVAILABLE = True
except ImportError:
    _API_AVAILABLE = False

__all__ = [
    # Base classes
    'BaseProvider',
    'MotifInstance',
    'MotifType',
    'ResidueSpec',
    'DatabaseInfo',
    # Registry
    'DatabaseRegistry',
    'get_registry',
    'initialize_registry',
    # Local Providers
    'RNA3DAtlasProvider',
    'RfamProvider',
    # API Providers
    'BGSUAPIProvider',
    'RfamAPIProvider',
    # Cache
    'CacheManager',
    'get_cache_manager',
    'initialize_cache_manager',
    # Config
    'PluginConfig',
    'SourceMode',
    'FreshnessPolicy',
    'CachePolicy',
    'get_config',
    'set_config',
    # Source Selector
    'SourceSelector',
    'get_source_selector',
    'initialize_source_selector',
]
