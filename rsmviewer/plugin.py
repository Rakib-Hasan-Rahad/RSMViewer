"""PyMOL entry point for the RSMViewer four-source motif workflow."""

from pymol import cmd
from pathlib import Path
from .gui import initialize_gui
from .utils import initialize_logger
from .database import initialize_registry, get_registry
from datetime import datetime
import sys

# ── Terminal encoding guard ──────────────────────────────────────────────────
# On Linux/macOS with LANG=C or Windows with a legacy code-page, sys.stdout may
# be ASCII-only.  Reconfigure it to UTF-8 with 'replace' so that Unicode output
# from rmv_help / rmv_db never raises UnicodeEncodeError in terminal mode.
try:
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    if hasattr(sys.stderr, 'reconfigure'):
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
# ─────────────────────────────────────────────────────────────────────────────


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
    last_updated = "13 September 2026"

    print("\n" + "=" * 80)
    print("RSMViewer")
    print("RNA Structural Motif Visualization and Comparative Analysis for PyMOL")
    print(f"Version 2.0.0 | Updated: {last_updated} | Compatible with PyMOL 2.x+")
    print("=" * 80)

    print("\nIntegrate motif annotations, Visualize motif instances, Compare structural variation")

    print("\nSUPPORTED ANNOTATION SOURCES:")
    print("   RNA3DMotifAtlas, Rfam, FR3D, RNAMotifScanX")

    print("\nQUICK START:")
    print("   rmv_fetch 1S72")
    print("   rmv_db RNA3DMotifAtlas")
    

    print("\nCOMMANDS & HELP:")
    print("   rmv_help                  # View all available commands")
    print("   rmv_db                    # List / select annotation sources")

    print("\n" + "=" * 80 + "\n")
    # Layer 1: Lock chain ID convention to auth_asym_id
    # This ensures PyMOL uses auth_asym_id as the 'chain' property when loading
    # CIF files, which matches the convention used by BGSU, Atlas, and other
    # motif annotation sources. The label_asym_id is stored in 'segi' as fallback.
    try:
        cmd.set("cif_use_auth", 1)
        logger.debug("Chain ID convention locked: cif_use_auth=1 (auth_asym_id)")
    except Exception as e:
        logger.warning(f"Could not set cif_use_auth: {e}")

    # Query groups ('as GROUP' in rmv_select/rmv_view/rmv_combine_groups)
    # are session-scoped only - clear any left over from a previous PyMOL
    # session so alias names never collide across restarts.
    try:
        from .database.motif_hierarchy_cache import get_hierarchy_cache
        get_hierarchy_cache().reset_all_aliases()
    except Exception as e:
        logger.debug(f"Alias store reset skipped: {e}")
    
    # Initialize database registry with all available providers
    try:
        database_dir = plugin_dir / 'motif_database'
        logger.debug(f"Initializing database registry from {database_dir}")
        
        # Initialize registry - this registers all available providers
        registry = initialize_registry(str(database_dir))
        
        # Build compact summary of registered databases
        providers = registry.get_all_providers()
        from .database.source_registry import get_source_registry
        source_registry = get_source_registry()
        source_names = [
            source_registry.get_public_name_for_provider(provider_id)
            or provider.info.name
            for provider_id, provider in providers.items()
        ]
        logger.success(f"Loaded {len(providers)} sources: {', '.join(source_names)}")
        
    except Exception as e:
        logger.error(f"Error initializing database registry: {e}")
        import traceback
        traceback.print_exc()
    
    # Initialize GUI and register commands. Never let a late registration
    # failure abort the whole plugin: the core commands register first, so the
    # pipeline stays usable even if an optional module fails to import.
    try:
        initialize_gui()
    except Exception as e:
        logger.error(f"GUI initialization incomplete (core commands may still work): {e}")
        import traceback
        traceback.print_exc()


# Module metadata
__all__ = ['__init_plugin__']
