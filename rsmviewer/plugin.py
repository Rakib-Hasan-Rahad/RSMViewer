"""
RSMViewer - Main Plugin Module
PyMOL plugin for visualizing RNA structural motifs with multi-database support.

Version: 1.0.0
Author: CBB LAB @Rakib Hasan Rahad
License: MIT

This module is the entry point when the plugin is loaded into PyMOL.
Supports 8 data sources:
- RNA 3D Motif Atlas [1], Rfam [2] (bundled, offline)
- BGSU RNA 3D Hub [3], Rfam API [4] (online)
- FR3D [5], RNAMotifScan [6], RNAMotifScanX [7], NoBIAS [8] (user annotations)

Usage in PyMOL:
    rmv_fetch 1S72                # Load PDB structure
    rmv_sources                   # Check available data sources
    rmv_db 3                      # Select data source (BGSU)
    rmv_load_motif                # Fetch motif data
    rmv_summary                   # Show motif types & counts
    rmv_summary HL                # Show hairpin loop instances
    rmv_show HL                   # Render all hairpin loops
    rmv_show HL 1                 # Zoom to specific instance
    rmv_help                      # Full command reference
"""

from pymol import cmd
from pathlib import Path
from .gui import initialize_gui
from .utils import initialize_logger
from .database import initialize_registry, get_registry
from datetime import datetime


def __init_plugin__(app):
    """
    Initialize plugin in PyMOL with multi-database support.
    
    This function is called by PyMOL when the plugin is first loaded.
    
    Initialization steps:
    1. Setup logging
    2. Initialize database registry with all available providers
    3. Register GUI commands
    4. Print welcome message with usage instructions
    
    Args:
        app: PyMOL application instance
    """
    # Initialize logger
    plugin_dir = Path(__file__).parent
    logger = initialize_logger(use_pymol_console=True)
    
    # Print welcome banner first
    last_updated = "02 April 2026"
    print("\n" + "="*80)
    print("\u250c" + " "*78 + "\u2510")
    print("\u2502" + " "*26 + "\U0001f9ec RSMViewer \U0001f9ec" + " "*38 + "\u2502")
    print("\u2502" + " "*78 + "\u2502")
    print("\u2502" + " Version 1.0.0" + " "*63 + "\u2502")
    print("\u2502" + " Last Updated: " + last_updated + " "*42 + "\u2502")
    print("\u2502" + " "*78 + "\u2502")
    print("\u2502" + " RSMViewer \u2014 RNA Structural Motif Viewer for PyMOL" + " "*26 + "\u2502")
    print("\u2502" + " Fast loading: Load PDB first, select source, fetch motifs" + " "*19 + "\u2502")
    print("\u2502" + " "*78 + "\u2502")
    print("\u2514" + " "*78 + "\u2518")
    print("="*80)
    print("\n\U0001f4ca AVAILABLE DATA SOURCES:")
    print("   \u2022 Local:       RNA 3D Atlas [1], Rfam [2] (offline)")
    print("   \u2022 Online:      BGSU RNA 3D Hub [3], Rfam API [4]")
    print("   \u2022 Annotations: FR3D [5], RNAMotifScan [6], RNAMotifScanX [7], NoBIAS [8]")
    print("\n\u26a1 QUICK START:")
    print("   rmv_fetch 1S72            # Step 1: Load PDB structure")
    print("   rmv_sources               # Step 2: Check available data sources")
    print("   rmv_db 3                  # Step 3: Select BGSU API source")
    print("   rmv_load_motif            # Step 4: Fetch motif data")
    print("   rmv_summary               # Step 5: Show motif types & counts")
    print("   rmv_summary HL            # Step 6: Show hairpin loop instances")
    print("   rmv_show HL               # Step 7: Render all hairpin loops")
    print("   rmv_show HL 1             # Step 8: Zoom to specific instance")
    print("\n\U0001f4da COMMANDS & HELP:")
    print("   rmv_help                  # All available commands")
    print("   rmv_sources               # List all data sources")
    print("   rmv_source info           # Show currently selected source")
    print("\n" + "="*80 + "\n")
    
    # Layer 1: Lock chain ID convention to auth_asym_id
    # This ensures PyMOL uses auth_asym_id as the 'chain' property when loading
    # CIF files, which matches the convention used by BGSU, Atlas, and other
    # motif annotation sources. The label_asym_id is stored in 'segi' as fallback.
    try:
        cmd.set("cif_use_auth", 1)
        logger.debug("Chain ID convention locked: cif_use_auth=1 (auth_asym_id)")
    except Exception as e:
        logger.warning(f"Could not set cif_use_auth: {e}")
    
    # Initialize database registry with all available providers
    try:
        database_dir = plugin_dir / 'motif_database'
        logger.debug(f"Initializing database registry from {database_dir}")
        
        # Initialize registry - this registers all available providers
        registry = initialize_registry(str(database_dir))
        
        # Get summary of registered databases
        providers = registry.get_all_providers()
        total_pdbs = sum(len(p.get_available_pdb_ids()) for p in providers.values())
        total_motif_types = sum(len(p.get_available_motif_types()) for p in providers.values())
        
        logger.success(f"Loaded {len(providers)} database(s)")
        for pid, provider in providers.items():
            logger.debug(f"  {pid}: {provider.info.name} - "
                        f"{len(provider.get_available_motif_types())} motif types, "
                        f"{len(provider.get_available_pdb_ids())} PDB structures")
        
    except Exception as e:
        logger.error(f"Error initializing database registry: {e}")
        import traceback
        traceback.print_exc()
    
    # Initialize GUI and register commands
    initialize_gui()


# Module metadata
__all__ = ['__init_plugin__']
