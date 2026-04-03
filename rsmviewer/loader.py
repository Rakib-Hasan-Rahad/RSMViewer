"""
RSMViewer - Loader Module
Handles loading structures and motif annotations using scalable database providers.

This module provides:
- StructureLoader: Loads PDB structures into PyMOL
- UnifiedMotifLoader: Loads motifs from any registered database provider
- VisualizationManager: Coordinates the complete visualization workflow

The loader now uses the database registry to support multiple databases
(RNA 3D Atlas, Rfam, etc.) with a unified interface.

Author: CBB LAB @Rakib Hasan Rahad
Version: 1.0.0
"""

import os
from pathlib import Path
from typing import Dict, List, Optional

from .utils import (
    PDBParser,
    MotifSelector,
    get_logger,
)
from .utils.selectors import sanitize_pymol_name
from . import colors
from .database import (
    get_registry,
    MotifInstance,
)


class StructureLoader:
    """Handles loading RNA structures into PyMOL."""
    
    def __init__(self, cmd):
        """
        Initialize loader.
        
        Args:
            cmd: PyMOL cmd module
        """
        self.cmd = cmd
        self.logger = get_logger()
        self.current_structure = None
        self.current_pdb_id = None
    
    def load_structure(self, pdb_id_or_path):
        """
        Load an RNA structure into PyMOL.
        
        Args:
            pdb_id_or_path (str): Either a PDB ID or local file path
        
        Returns:
            str: Name of loaded structure in PyMOL, or None if failed
        """
        try:
            # Determine if it's a PDB ID or file path
            if os.path.isfile(pdb_id_or_path):
                # Local file
                self.logger.info(f"Loading structure from file: {pdb_id_or_path}")
                structure_name = Path(pdb_id_or_path).stem
                self.cmd.load(pdb_id_or_path, structure_name)
                pdb_id = PDBParser.extract_pdb_id(pdb_id_or_path)
            else:
                # Assume PDB ID
                pdb_id = pdb_id_or_path.strip().upper()
                if not PDBParser.is_valid_pdb_id(pdb_id):
                    self.logger.error(f"Invalid PDB ID format: {pdb_id}")
                    return None
                
                self.logger.info(f"Downloading structure from RCSB: {pdb_id}")
                structure_name = pdb_id
                # Delete any existing object with the same name to avoid "loading mmCIF into existing object" error
                try:
                    self.cmd.delete(structure_name)
                except:
                    pass  # Object doesn't exist, that's fine
                self.cmd.fetch(pdb_id, structure_name)
            
            self.current_structure = structure_name
            self.current_pdb_id = pdb_id
            
            self.logger.success(f"Loaded structure: {structure_name} (PDB: {pdb_id})")
            return structure_name
            
        except Exception as e:
            self.logger.error(f"Failed to load structure: {e}")
            return None
    
    def get_current_structure(self):
        """Get name of currently loaded structure."""
        return self.current_structure
    
    def get_current_pdb_id(self):
        """Get PDB ID of currently loaded structure."""
        return self.current_pdb_id


class UnifiedMotifLoader:
    """
    Unified motif loader that works with any database provider.
    
    Uses the source selector to automatically find the best data source:
    1. Local bundled databases (fast, offline)
    2. BGSU RNA 3D Hub API (comprehensive, ~3000+ PDBs)
    3. Rfam API (named motifs)
    """
    
    def __init__(self, cmd, database_dir: str):
        """
        Initialize motif loader.
        
        Args:
            cmd: PyMOL cmd module
            database_dir (str): Path to motif database directory
        """
        self.cmd = cmd
        self.database_dir = database_dir
        self.logger = get_logger()
        self.selector = MotifSelector(cmd)
        self.loaded_motifs: Dict[str, Dict] = {}  # Track loaded motif objects
        self._last_source_used: Optional[str] = None
        
        # Get registry
        self._registry = get_registry()
    
    def load_motifs(self, structure_name: str, pdb_id: str,
                   provider_id: Optional[str] = None,
                   force_refresh: bool = False) -> Dict:
        """
        Load all available motifs for a structure.
        
        Uses smart source selection: tries local first, then APIs if needed.
        
        Args:
            structure_name (str): Name of structure in PyMOL
            pdb_id (str): PDB ID to look up motifs for
            provider_id (str): Optional specific provider to use (auto-select if None)
            force_refresh (bool): Force re-fetch from API (ignore cache)
        
        Returns:
            dict: Dictionary of loaded motifs: {motif_type: {details}}
        """
        try:
            self.loaded_motifs = {}
            pdb_id = pdb_id.upper()
            
            # Try to use source selector for smart source selection
            from .database import get_source_selector
            source_selector = get_source_selector()
            
            if source_selector and not provider_id:
                # Use smart source selection
                available_motifs, source_used = source_selector.get_motifs_for_pdb(
                    pdb_id, 
                    source_override=provider_id,
                    force_refresh=force_refresh
                )
                self._last_source_used = source_used
                source_name = source_used or "unknown"
            else:
                # Fall back to registry-based provider selection
                if provider_id:
                    provider = self._registry.get_provider(provider_id)
                else:
                    provider = self._registry.get_active_provider()
                
                if not provider:
                    self.logger.error("No database provider available")
                    return {}
                
                available_motifs = provider.get_motifs_for_pdb(pdb_id)
                source_name = provider.info.name if hasattr(provider, 'info') else provider_id
                self._last_source_used = provider_id
            
            if not available_motifs:
                self.logger.warning(f"No motifs found for PDB {pdb_id}")
                self.logger.info("Tip: This PDB may not have RNA motif annotations in any database")
                return {}
            
            total_count = sum(len(instances) for instances in available_motifs.values())
            self.logger.info(f"Found {total_count} motifs in {pdb_id} (source: {source_name})")
            
            # Process each motif type
            for motif_type, instances in available_motifs.items():
                try:
                    # Handle prefixed motif types from combined sources (e.g., "bgsu_api:HL")
                    display_type = motif_type.split(':')[-1] if ':' in motif_type else motif_type
                    self._load_motif_type(structure_name, pdb_id, display_type, instances)
                except Exception as e:
                    self.logger.error(f"Error loading {motif_type} motifs: {e}")
                    continue
            
            return self.loaded_motifs
            
        except Exception as e:
            self.logger.error(f"Failed to load motifs: {e}")
            import traceback
            traceback.print_exc()
            return {}
    
    def get_last_source_used(self) -> Optional[str]:
        """Get the data source used for the last load operation."""
        return self._last_source_used
    
    def _load_motif_type(self, structure_name: str, pdb_id: str,
                        motif_type: str, instances: List,
                        use_direct_coloring: bool = True) -> None:
        """
        Load a specific motif type and visualize in PyMOL.
        
        Creates PyMOL objects (visible in right panel) AND colors residues
        directly on the structure (to avoid z-fighting stripes).
        
        Args:
            structure_name: PyMOL structure name
            pdb_id: PDB ID
            motif_type: Type of motif (HL, IL, GNRA, etc.)
            instances: List of MotifInstance objects
            use_direct_coloring: If True, also color residues directly on structure
        """
        if not instances:
            return
        
        # Build motif_list in format MotifSelector expects
        motif_list: List[Dict] = []
        motif_details = []
        
        for instance in instances:
            if not instance.residues:
                continue
            
            # Convert to legacy format for selector
            legacy_entries = instance.to_legacy_format()
            motif_list.extend(legacy_entries)
            
            # Store with deepcopy of metadata to ensure independence
            from copy import deepcopy
            
            metadata_to_store = deepcopy(instance.metadata) if instance.metadata else {}
            
            motif_details.append({
                'motif_id': instance.motif_id,
                'instance_id': instance.instance_id,
                'residues': [r.to_tuple() for r in instance.residues],
                'annotation': instance.annotation,
                'metadata': metadata_to_store,
            })
        
        if not motif_list:
            self.logger.debug(f"No residues found for {motif_type} motifs in {pdb_id}")
            return
        
        # SORT: Order instances by minimum residue number for consistent numbering across all commands
        def get_min_residue(detail):
            residues = detail.get('residues', [])
            if not residues:
                return float('inf')
            min_resi = float('inf')
            for res in residues:
                if isinstance(res, tuple) and len(res) >= 2:
                    resi = res[1]
                    if isinstance(resi, int):
                        min_resi = min(min_resi, resi)
            return min_resi if min_resi != float('inf') else float('inf')
        
        motif_details.sort(key=get_min_residue)
        
        motif_type_upper = motif_type.upper()
        color_rgb = colors.get_color(motif_type_upper)
        
        # Build combined selection for all residues of this motif type
        from .utils.parser import SelectionParser
        all_selections = []
        for motif in motif_list:
            chain = motif.get('chain')
            residues = motif.get('residues')
            sel = SelectionParser.create_selection_string(chain, residues)
            if sel:
                all_selections.append(f"({sel})")
        
        if all_selections:
            combined_sel = " or ".join(all_selections)
            main_motif_sel = f"(model {structure_name}) and ({combined_sel})"
            # Hide the cartoon on main structure for these residues
            self.cmd.hide('cartoon', main_motif_sel)
        
        # Create PyMOL object (visible in right panel)
        obj_name = self.selector.create_motif_class_object(
            structure_name,
            motif_type_upper,
            motif_list,
            pdb_id=pdb_id,
        )
        
        if obj_name:
            # Show cartoon on the motif object
            self.cmd.show('cartoon', obj_name)
            
            # Apply consistent representation settings for uniform appearance
            self.cmd.set('cartoon_nucleic_acid_mode', 4, obj_name, quiet=1)  # Simple tube mode
            self.cmd.set('cartoon_tube_radius', 0.4, obj_name, quiet=1)
            
            # Color the object with the motif color
            colors.set_motif_color_in_pymol(self.cmd, obj_name, motif_type_upper)
            
            # Object stays ENABLED (visible in panel and on screen)
            # No z-fighting because main structure doesn't render these residues
            
            self.loaded_motifs[motif_type_upper] = {
                'object_name': obj_name,
                'structure_name': structure_name,
                'pdb_id': pdb_id,
                'count': len(instances),
                'visible': True,
                'motifs': motif_list,
                'motif_details': motif_details,
                'color_rgb': color_rgb,
                'main_selection': main_motif_sel if all_selections else None,
            }
            
            self.logger.success(f"Loaded {len(instances)} {motif_type_upper} motifs")
    
    def toggle_motif_type(self, motif_type: str, visible: bool) -> bool:
        """
        Toggle visibility of a motif type.
        
        Shows/hides the motif object.
        
        Args:
            motif_type (str): Motif type (e.g., 'HL', 'IL', 'GNRA')
            visible (bool): True to show, False to hide
        
        Returns:
            bool: True if successful
        """
        motif_type = str(motif_type).upper().strip()
        motif_type = motif_type.replace('-', '_').replace(' ', '_')
        
        if motif_type not in self.loaded_motifs:
            self.logger.warning(f"Motif type {motif_type} not loaded")
            return False
        
        info = self.loaded_motifs[motif_type]
        obj_name = info['object_name']
        
        if visible:
            # Show the motif object
            self.cmd.enable(obj_name)
            self.cmd.show('cartoon', obj_name)
        else:
            # Hide the motif object
            self.cmd.disable(obj_name)
        
        self.loaded_motifs[motif_type]['visible'] = visible
        return True
    
    def get_loaded_motifs(self) -> Dict:
        """Get dictionary of loaded motifs."""
        return self.loaded_motifs
    
    def clear_motifs(self) -> None:
        """Clear all loaded motif objects from PyMOL."""
        try:
            for motif_type, info in self.loaded_motifs.items():
                obj_name = info['object_name']
                
                # Delete the object
                try:
                    self.cmd.delete(obj_name)
                except:
                    pass
            
            self.loaded_motifs = {}
            self.logger.info("Cleared all motif objects")
        except Exception as e:
            self.logger.error(f"Failed to clear motifs: {e}")
    
    def reload_motifs(self, structure_name: str, pdb_id: str,
                     provider_id: Optional[str] = None) -> Dict:
        """
        Reload motifs (clear and reload).
        
        Args:
            structure_name (str): Name of structure in PyMOL
            pdb_id (str): PDB ID
            provider_id (str): Optional provider ID
        
        Returns:
            dict: Loaded motifs
        """
        self.clear_motifs()
        return self.load_motifs(structure_name, pdb_id, provider_id)
    
    def get_available_motif_types(self, pdb_id: str,
                                  provider_id: Optional[str] = None) -> List[str]:
        """Get list of motif types available for a PDB."""
        pdb_id = pdb_id.upper()
        
        if provider_id:
            provider = self._registry.get_provider(provider_id)
        else:
            provider = self._registry.get_active_provider()
        
        if not provider:
            return []
        
        motifs = provider.get_motifs_for_pdb(pdb_id)
        return sorted(list(motifs.keys()))
    
    def get_registry(self):
        """Get the database registry."""
        return self._registry


class VisualizationManager:
    """High-level manager for the entire visualization workflow."""
    
    def __init__(self, cmd, database_dir: str):
        """
        Initialize visualization manager.
        
        Args:
            cmd: PyMOL cmd module
            database_dir (str): Path to motif database directory
        """
        self.cmd = cmd
        self.database_dir = database_dir
        self.structure_loader = StructureLoader(cmd)
        self.motif_loader = UnifiedMotifLoader(cmd, database_dir)
        self.logger = get_logger()
        self._current_provider_id = None
    
    def setup_clean_visualization(self, structure_name: str,
                                 background_color: Optional[str] = None) -> None:
        """
        Set up clean RNA visualization with uniform color.
        
        Workflow:
        1. Hide everything
        2. Select all polymer.nucleic (all RNA chains)
        3. Show cartoon representation
        4. Set cartoon nucleic acid mode
        5. Color uniformly with background_color
        
        Args:
            structure_name (str): Name of structure in PyMOL
            background_color (str): Color for the RNA (default: 'gray80')
        """
        try:
            if background_color is None:
                background_color = colors.NON_MOTIF_COLOR or 'gray80'
            
            # Hide everything first
            self.cmd.hide('everything', 'all')
            self.logger.debug("Hidden all objects")
            
            # Select ALL polymer.nucleic
            rna_selection = f"{structure_name}_rna"
            self.cmd.select(rna_selection, f"model {structure_name} and polymer.nucleic")
            
            # Show cartoon representation
            self.cmd.show('cartoon', rna_selection)
            
            # Set consistent cartoon nucleic acid settings for uniform appearance
            self.cmd.set('cartoon_nucleic_acid_mode', 4)  # Simple tube mode
            self.cmd.set('cartoon_tube_radius', 0.4)
            
            # Color uniformly
            self.cmd.color(background_color, rna_selection)
            self.logger.info(f"Visualization: All RNA shown as {background_color} cartoon")
            
            # Clean up temporary selection
            self.cmd.delete(rna_selection)
            
        except Exception as e:
            self.logger.error(f"Failed to set up visualization: {e}")
    
    def load_and_visualize(self, pdb_id_or_path: str,
                          background_color: Optional[str] = None,
                          provider_id: Optional[str] = None) -> Dict:
        """
        Complete workflow: load structure and visualize all motifs.
        
        Args:
            pdb_id_or_path (str): PDB ID or file path
            background_color (str): Color for RNA backbone (default: 'gray80')
            provider_id (str): Optional database provider ID
        
        Returns:
            dict: Loaded motifs, or empty dict if failed
        """
        # Load structure
        structure_name = self.structure_loader.load_structure(pdb_id_or_path)
        if not structure_name:
            return {}
        
        pdb_id = self.structure_loader.get_current_pdb_id()
        
        # Set up clean visualization
        self.setup_clean_visualization(structure_name, background_color)
        
        # Load motifs from specified or active provider
        motifs = self.motif_loader.load_motifs(structure_name, pdb_id, provider_id)
        
        if provider_id:
            self._current_provider_id = provider_id
        
        # Print detailed summary table to PyMOL console
        if motifs:
            self._print_motif_summary_table(pdb_id, motifs, provider_id)
            
            # NOW visualize all motifs (highlight and create objects)
            self.show_all_motifs()
        
        return motifs
    
    def switch_database(self, provider_id: str) -> bool:
        """
        Switch to a different database provider.
        
        Args:
            provider_id: ID of the provider to switch to
            
        Returns:
            True if successful
        """
        registry = self.motif_loader.get_registry()
        if registry.set_active_provider(provider_id):
            self._current_provider_id = provider_id
            self.logger.info(f"Switched to database: {provider_id}")
            return True
        return False
    
    def reload_with_database(self, provider_id: str,
                            background_color: Optional[str] = None) -> Dict:
        """
        Reload current structure with a different database.
        
        Args:
            provider_id: Database provider to use
            background_color: Optional background color
            
        Returns:
            Loaded motifs
        """
        structure_name = self.structure_loader.get_current_structure()
        pdb_id = self.structure_loader.get_current_pdb_id()
        
        if not structure_name or not pdb_id:
            self.logger.error("No structure loaded")
            return {}
        
        # Clear existing motifs
        self.motif_loader.clear_motifs()
        
        # Set up visualization again
        self.setup_clean_visualization(structure_name, background_color)
        
        # Switch database and load
        self.switch_database(provider_id)
        return self.motif_loader.load_motifs(structure_name, pdb_id, provider_id)
    
    def get_structure_info(self) -> Dict:
        """Get current structure and motif info."""
        registry = self.motif_loader.get_registry()
        active_provider = registry.get_active_provider()
        
        return {
            'structure': self.structure_loader.get_current_structure(),
            'pdb_id': self.structure_loader.get_current_pdb_id(),
            'motifs': self.motif_loader.get_loaded_motifs(),
            'database': active_provider.info.name if active_provider else None,
            'database_id': self._current_provider_id,
        }
    
    def get_available_databases(self) -> List[Dict]:
        """Get list of available database providers."""
        registry = self.motif_loader.get_registry()
        return [
            {
                'id': pid,
                'name': provider.info.name,
                'description': provider.info.description,
                'motif_types': len(provider.get_available_motif_types()),
                'pdb_count': len(provider.get_available_pdb_ids()),
                'active': pid == self._current_provider_id,
            }
            for pid, provider in registry.get_all_providers().items()
        ]
    
    def get_available_motif_summary(self, pdb_id: str) -> str:
        """Get summary of available motifs for a PDB."""
        motif_types = self.motif_loader.get_available_motif_types(pdb_id)
        if not motif_types:
            return f"No motifs found for {pdb_id}"
        return f"Available motifs: {', '.join(motif_types)}"
    
    def _print_motif_summary_table(self, pdb_id: str, motifs: Dict,
                                   provider_id: Optional[str] = None) -> None:
        """
        Print a simplified summary table of loaded motifs to PyMOL console.
        Shows only motif types and their frequencies.
        
        Args:
            pdb_id: PDB ID
            motifs: Dictionary of loaded motifs
            provider_id: Database provider used
        """
        # Map source IDs to friendly names
        source_names = {
            'atlas': 'RNA 3D Motif Atlas (Local)',
            'rfam': 'Rfam (Local)',
            'bgsu_api': 'BGSU RNA 3D Hub (API)',
            'rfam_api': 'Rfam (API)',
            'fr3d': 'FR3D User Annotations',
            'rnamotifscan': 'RNAMotifScan User Annotations',
        }
        
        db_name = "Unknown"
        
        # Use provider_id directly if passed by caller (gui.py already resolves it)
        if provider_id and provider_id != "Unknown":
            db_name = provider_id
        else:
            # Try to get from GUI state first
            try:
                from .gui import gui
                if gui and hasattr(gui, 'current_source_mode'):
                    mode = gui.current_source_mode
                    
                    if mode == 'user' and hasattr(gui, 'current_user_tool'):
                        tool = gui.current_user_tool
                        db_name = source_names.get(tool, f"{tool} User Annotations")
                    
                    elif mode == 'local':
                        if hasattr(gui, 'current_local_source') and gui.current_local_source:
                            db_name = source_names.get(gui.current_local_source, 
                                                      f"{gui.current_local_source} (Local)")
                        else:
                            db_name = 'Local Databases (Atlas + Rfam)'
                    
                    elif mode == 'web' or mode == 'bgsu':
                        if hasattr(gui, 'current_web_source') and gui.current_web_source:
                            db_name = source_names.get(gui.current_web_source, 
                                                      f"{gui.current_web_source} (API)")
                        else:
                            db_name = 'Online APIs'
                    
                    elif mode == 'auto':
                        db_name = 'Auto-selected (Local First → API)'
                    
                    elif mode == 'all':
                        db_name = 'Combined (All Sources)'
                    
                    elif mode == 'combine':
                        db_name = 'Combined (Multiple Sources)'
            except:
                # Fallback if GUI not available
                pass
        
        # If still unknown, try last_source
        if db_name == "Unknown":
            last_source = self.motif_loader.get_last_source_used()
            if last_source:
                if ',' in last_source:
                    sources = [source_names.get(s.strip(), s.strip()) for s in last_source.split(',')]
                    db_name = ' + '.join(sources)
                else:
                    db_name = source_names.get(last_source, last_source)
        
        # Build the table
        print("\n" + "=" * 50)
        print(f"  MOTIF SUMMARY - {pdb_id}")
        print("=" * 50)
        print(f"  Database: {db_name}")
        print("-" * 50)
        
        # Header
        print(f"  {'MOTIF TYPE':<20} {'INSTANCES':>12}")
        print("-" * 50)
        
        total_motifs = 0
        
        # Sort by count (descending) - highest to lowest
        for motif_type, info in sorted(motifs.items(), key=lambda x: x[1].get('count', 0), reverse=True):
            count = info.get('count', 0)
            total_motifs += count
            print(f"  {motif_type:<20} {count:>12}")
        
        print("-" * 50)
        print(f"  {'TOTAL':<20} {total_motifs:>12}")
        print("=" * 50)
        print("\n  Next steps:")
        if total_motifs > 0:
            # Find the motif type with highest count to suggest
            first_motif = None
            sorted_by_count = sorted(motifs.items(), key=lambda x: x[1].get('count', 0), reverse=True)
            if sorted_by_count:
                first_motif = sorted_by_count[0][0]
            if first_motif:
                print(f"    rmv_summary {first_motif:<18}  Show {first_motif} instances")
                print(f"    rmv_show {first_motif:<20}  Highlight & render {first_motif}")
                print(f"    rmv_super {first_motif:<19}  Superimpose {first_motif} instances")
            print(f"    rmv_show ALL             Show all motif types with objects")
        print("=" * 50 + "\n")
    
    def show_motif_type(self, motif_type: str,
                        filter_pdb: str = '', filter_suffix: str = '',
                        padding: int = 0) -> bool:
        """
        Show only a specific motif type highlighted, with full structure visible in gray.
        
        Workflow:
        1. Create PyMOL object if needed (for visibility in object panel)
        2. Hide ALL separate motif objects (avoid overlap/stripe artifacts)
        3. Show full PDB structure uniformly
        4. Color the entire structure gray80
        5. Color ONLY the selected motif residues in their color (within the structure)
        
        Args:
            motif_type: Motif type to show (e.g., 'GNRA', 'HL')
            filter_pdb: When set, only show instances matching this PDB ID
            filter_suffix: When set, only show instances matching this source suffix
            
        Returns:
            True if successful
        """
        motif_type = motif_type.upper().strip()
        loaded_motifs = self.motif_loader.get_loaded_motifs()
        
        if motif_type not in loaded_motifs:
            self.logger.error(f"Motif type '{motif_type}' not loaded")
            self.logger.info(f"Available: {', '.join(loaded_motifs.keys())}")
            return False
        
        # Get structure name and motif details
        info = loaded_motifs[motif_type]
        structure_name = info.get('structure_name')
        motif_details = info.get('motif_details', [])

        # Apply source filter when requested (show only current source data)
        if filter_pdb:
            motif_details = [
                d for d in motif_details
                if d.get('_pdb_id', info.get('pdb_id', '')) == filter_pdb
                and d.get('_source_suffix', info.get('source_suffix', '')) == filter_suffix
            ]
            if not motif_details:
                self.logger.error(f"No {motif_type} instances for current source")
                return False

        main_selection = info.get('main_selection')
        motif_list = info.get('motifs', [])  # For PyMOL object creation
        # Use filter_suffix for naming when provided — correct for combine mode
        source_suffix = filter_suffix if filter_suffix else info.get('source_suffix', '')
        
        if not structure_name:
            self.logger.error("No structure name found")
            return False
        
        # Step 0: Create PyMOL object if it doesn't exist (needed for object panel visibility)
        obj_name = info.get('object_name')
        pdb_id = info.get('pdb_id', self.structure_loader.get_current_pdb_id() or '')
        if not obj_name and motif_list:
            # Create the PyMOL object (like HL_ALL_1S72_S3)
            obj_name = self.motif_loader.selector.create_motif_class_object(
                structure_name,
                motif_type,
                motif_list,
                source_suffix=source_suffix,
                pdb_id=pdb_id,
            )
            if obj_name:
                # Color the object with the motif color
                colors.set_motif_color_in_pymol(self.cmd, obj_name, motif_type)
                # Update the info dict
                info['object_name'] = obj_name
                loaded_motifs[motif_type] = info
                self.logger.debug(f"Created PyMOL object: {obj_name}")

        # When padding is active, build a _P object from padded selections
        if padding > 0 and structure_name and motif_details:
            pad_obj = (obj_name + "_P") if obj_name else sanitize_pymol_name(
                f"{motif_type}_ALL_{pdb_id}{source_suffix}_P")
            from .utils.parser import SelectionParser as _SP
            pad_sels = []
            for _det in motif_details:
                _res = _det.get('residues', [])
                if not _res:
                    continue
                _cr = {}
                for _r in _res:
                    if isinstance(_r, tuple) and len(_r) >= 3:
                        _cr.setdefault(_r[2], []).append(_r[1])
                for _ch, _rl in _cr.items():
                    expanded = set()
                    for _v in _rl:
                        for _off in range(-padding, padding + 1):
                            expanded.add(_v + _off)
                    _s = _SP.create_selection_string(_ch, sorted(expanded))
                    if _s:
                        pad_sels.append(f"({_s})")
            if pad_sels:
                try:
                    existing = self.cmd.get_object_list()
                    if pad_obj in existing:
                        self.cmd.delete(pad_obj)
                    pad_full = f"(model {structure_name}) and (" + " or ".join(pad_sels) + ")"
                    self.cmd.create(pad_obj, pad_full)
                    self.cmd.set('cartoon_nucleic_acid_mode', 4, pad_obj, quiet=1)
                    self.cmd.set('cartoon_tube_radius', 0.4, pad_obj, quiet=1)
                    colors.set_motif_color_in_pymol(self.cmd, pad_obj, motif_type)
                except Exception as e:
                    self.logger.debug(f"Could not create padded object: {e}")
        
        # Hide ALL separate motif objects (prevents overlap/stripes)
        for mt, mt_info in loaded_motifs.items():
            obj_name_check = mt_info.get('object_name')
            if obj_name_check:
                self.cmd.disable(obj_name_check)
        
        # Hide any previously created instance objects
        for obj in self.cmd.get_object_list():
            for mt in loaded_motifs.keys():
                mt_sanitized = sanitize_pymol_name(mt)
                suffix = loaded_motifs[mt].get('source_suffix', '')
                mt_pdb_id = loaded_motifs[mt].get('pdb_id', '')
                pdb_tag = f"_{mt_pdb_id}" if mt_pdb_id else ""
                prefix = sanitize_pymol_name(f"{mt_sanitized}{pdb_tag}{suffix}")
                if obj.startswith(f"{prefix}_") and obj[len(prefix)+1:].isdigit():
                    self.cmd.disable(obj)
        
        # Show the full structure with uniform representation
        self.cmd.enable(structure_name)
        self.cmd.show('cartoon', f"model {structure_name} and polymer.nucleic")
        self.cmd.set('cartoon_nucleic_acid_mode', 4, f"model {structure_name}")
        self.cmd.set('cartoon_tube_radius', 0.4, f"model {structure_name}")
        
        # Step 3: Color the ENTIRE structure gray80 first
        self.cmd.color('gray80', f"model {structure_name} and polymer.nucleic")
        
        # Color the selected motif residues in their color
        if motif_details:
            # Show the motif residues (in case they were hidden)
            if main_selection:
                self.cmd.show('cartoon', main_selection)
            
            # Color each instance individually to avoid PyMOL selection string length limits
            # (Large "or" selections with 100+ instances can exceed PyMOL's parsing limits)
            from .utils.parser import SelectionParser
            for detail in motif_details:
                residues = detail.get('residues', [])
                if not residues:
                    continue
                
                # Build selection for this individual instance
                chain_residues = {}
                for res in residues:
                    if isinstance(res, tuple) and len(res) >= 3:
                        nucleotide, resi, chain = res[0], res[1], res[2]
                        if chain not in chain_residues:
                            chain_residues[chain] = []
                        chain_residues[chain].append(resi)
                
                # Create selection for this instance and color it
                selections = []
                for chain, resi_list in chain_residues.items():
                    if padding > 0:
                        expanded = set()
                        for r in resi_list:
                            for offset in range(-padding, padding + 1):
                                expanded.add(r + offset)
                        resi_list = sorted(expanded)
                    sel = SelectionParser.create_selection_string(chain, sorted(resi_list))
                    if sel:
                        selections.append(f"({sel})")
                
                if selections:
                    combined_sel = " or ".join(selections)
                    instance_sel = f"(model {structure_name}) and ({combined_sel})"
                    # Color this instance with the motif color
                    colors.set_motif_color_in_pymol(self.cmd, instance_sel, motif_type)
        
        # Print instance table
        self._print_motif_instance_table(motif_type, motif_details)
        
        pad_note = f" (padding=±{padding})" if padding > 0 else ""
        self.logger.success(f"Showing {len(motif_details)} {motif_type} instances{pad_note}")
        
        # Print follow-up suggestions
        print("  Next steps:")
        print(f"    rmv_show {motif_type} <NO>         Zoom to specific instance (1-{len(motif_details)})")
        print(f"    rmv_super {motif_type}             Superimpose all {motif_type} instances")
        print(f"    rmv_show <OTHER_MOTIF>       Show different motif type")
        print(f"    rmv_show ALL                 Show all motif types")
        print()
        return True
    
    def _create_single_instance_object(self, motif_type: str, instance_no: int,
                                         detail: Dict, structure_name: str, source_suffix: str = '') -> bool:
        """
        Create a single PyMOL object for one motif instance (on-demand).
        
        Args:
            motif_type: Motif type (e.g., 'HL', 'GNRA')
            instance_no: 1-indexed instance number
            detail: Single motif instance detail dict
            structure_name: Name of the structure in PyMOL
            source_suffix: Source suffix for object naming (e.g., '_S3')
            
        Returns:
            True if successful
        """
        from .utils.parser import SelectionParser
        
        residues = detail.get('residues', [])
        
        if not residues:
            return False
        
        # Build selection for this instance
        chain_residues = {}
        for res in residues:
            if isinstance(res, tuple) and len(res) >= 3:
                nucleotide, resi, chain = res[0], res[1], res[2]
                if chain not in chain_residues:
                    chain_residues[chain] = []
                chain_residues[chain].append(resi)
        
        # Create selection
        selections = []
        for chain, resi_list in chain_residues.items():
            sel = SelectionParser.create_selection_string(chain, sorted(resi_list))
            if sel:
                selections.append(f"({sel})")
        
        if not selections:
            return False
        
        combined_sel = " or ".join(selections)
        instance_sel = f"(model {structure_name}) and ({combined_sel})"
        
        # Create object name: MOTIF_NO_PDB_S3 (e.g., GNRA_1_1S72_S3, GNRA_2_4V9F_S7)
        pdb_id = self.structure_loader.get_current_pdb_id() or ''
        pdb_tag = f"_{pdb_id}" if pdb_id else ""
        obj_name = sanitize_pymol_name(f"{motif_type}_{instance_no}{pdb_tag}{source_suffix}")
        
        # Create the object
        try:
            self.cmd.create(obj_name, instance_sel)
            self.cmd.show('cartoon', obj_name)
            self.cmd.set('cartoon_nucleic_acid_mode', 4, obj_name)
            self.cmd.set('cartoon_tube_radius', 0.4, obj_name)
            colors.set_motif_color_in_pymol(self.cmd, obj_name, motif_type)
            return True
        except Exception as e:
            self.logger.debug(f"Could not create object {obj_name}: {e}")
            return False

    def _get_ranges_from_chainbreak(self, residues, chainbreak, metadata):
        """
        Get residue ranges from official chainbreak metadata.
        Uses alignment positions to determine strand boundaries.
        
        Args:
            residues: List of (nuc, resi, chain) tuples in ALIGNMENT ORDER
            chainbreak: List of strings like ["10"] or ["4", "8"], or empty []
            metadata: Dict that may contain 'alignment' key
            
        Returns:
            String like "0:965-968, 0:1000-1003" for IL with 2 strands
        """
        # Handle empty chainbreak as single strand (e.g., HL motifs)
        if isinstance(chainbreak, list) and len(chainbreak) == 0:
            # Single strand - compute min/max for all residues
            chains_dict = {}
            for res in residues:
                if isinstance(res, tuple) and len(res) >= 3:
                    nuc, resi, chain = res[0], res[1], res[2]
                    if chain not in chains_dict:
                        chains_dict[chain] = []
                    chains_dict[chain].append(resi)
            
            range_parts = []
            for chain in sorted(chains_dict.keys()):
                if chains_dict[chain]:
                    min_resi = min(chains_dict[chain])
                    max_resi = max(chains_dict[chain])
                    range_parts.append(f"{chain}:{min_resi}-{max_resi}")
            
            return ', '.join(range_parts) if range_parts else None
        
        if not chainbreak:
            return None
        
        try:
            # Parse chainbreak positions (1-indexed alignment positions)
            break_positions = sorted([int(x) for x in chainbreak])
        except (ValueError, TypeError):
            return None
        
        # Create strand groups based on alignment positions
        # Chainbreak position N means: positions 1..N in strand 1, N+1..M in strand 2, etc.
        groups = []
        start_pos = 1
        
        for break_pos in break_positions:
            # Collect residues with alignment positions in range [start_pos, break_pos]
            # Since residues list is indexed 0-based but alignment is 1-indexed
            end_pos = break_pos
            group = []
            for idx in range(len(residues)):
                align_pos = idx + 1  # 1-indexed alignment position
                if start_pos <= align_pos <= end_pos:
                    group.append(residues[idx])
            if group:
                groups.append(group)
            start_pos = break_pos + 1
        
        # Add remaining residues as final strand
        if start_pos <= len(residues):
            group = residues[start_pos - 1:]
            if group:
                groups.append(group)
        
        # For each strand group, find min/max residue numbers per chain
        range_parts = []
        for group in groups:
            if not group:
                continue
            
            # Group by chain in this strand
            chains_in_strand = {}
            for res in group:
                if isinstance(res, tuple) and len(res) >= 3:
                    nuc, resi, chain = res[0], res[1], res[2]
                    if chain not in chains_in_strand:
                        chains_in_strand[chain] = []
                    chains_in_strand[chain].append(resi)
            
            # For each chain, calculate min-max (IGNORE gaps within strand)
            for chain in sorted(chains_in_strand.keys()):
                resi_values = chains_in_strand[chain]
                if resi_values:
                    min_resi = min(resi_values)
                    max_resi = max(resi_values)
                    range_parts.append(f"{chain}:{min_resi}-{max_resi}")
        
        return ', '.join(range_parts) if range_parts else None
    
    def _get_nucleotides_for_strands(self, residues, metadata, pos_to_nuc):
        """
        Build nucleotide string for all strands using chainbreak metadata.
        
        Args:
            residues: List of (nuc, resi, chain) tuples in alignment order
            metadata: Dict that may contain 'chainbreak'
            pos_to_nuc: Dict mapping (chain, resi) to nucleotide
            
        Returns:
            String like "AUNAA...CUCGG..." (up to 25 chars) with per-strand nucleotides
        """
        all_nucs = []
        
        # If chainbreak exists, use it to group strands
        if metadata and 'chainbreak' in metadata:
            chainbreak = metadata['chainbreak']
            try:
                break_positions = sorted([int(x) for x in chainbreak])
            except (ValueError, TypeError):
                break_positions = []
            
            if break_positions:
                # Split into groups
                groups = []
                start_idx = 0
                for break_pos in break_positions:
                    end_idx = break_pos
                    if start_idx < len(residues):
                        groups.append(residues[start_idx:end_idx])
                    start_idx = end_idx
                if start_idx < len(residues):
                    groups.append(residues[start_idx:])
                
                # Extract nucleotides per strand
                for group in groups:
                    for res in group:
                        if isinstance(res, tuple) and len(res) >= 3:
                            nuc, resi, chain = res[0], res[1], res[2]
                            all_nucs.append(nuc if nuc and nuc != '-' else 'N')
        
        # Fallback: use all residues in order
        if not all_nucs:
            for res in residues:
                if isinstance(res, tuple) and len(res) >= 1:
                    nuc = res[0]
                    all_nucs.append(nuc if nuc and nuc != '-' else 'N')
        
        # Format output
        nucs_str = ''.join(all_nucs[:25])
        if len(all_nucs) > 25:
            nucs_str += '...'
        
        return nucs_str if nucs_str else '-'
    
    def _identify_strands(self, residues_with_pos):
        """
        Fallback: Identify separate strands by detecting gaps in residue numbers.
        Used when chainbreak metadata is not available.
        
        Args:
            residues_with_pos: List of (nucleotide, resi, chain) tuples
            
        Returns:
            Dict mapping chain to list of strand ranges [(start, end), ...]
        """
        strands_by_chain = {}
        
        for res in residues_with_pos:
            if isinstance(res, tuple) and len(res) >= 3:
                nucleotide, resi, chain = res[0], res[1], res[2]
                
                if chain not in strands_by_chain:
                    strands_by_chain[chain] = []
                
                strands_by_chain[chain].append(resi)
        
        # For each chain, identify contiguous ranges (strands)
        strand_ranges = {}
        for chain in strands_by_chain:
            positions = sorted(set(strands_by_chain[chain]))
            
            if not positions:
                continue
            
            ranges = []
            current_start = positions[0]
            current_end = positions[0]
            
            for pos in positions[1:]:
                if pos == current_end + 1:
                    # Contiguous, extend current range
                    current_end = pos
                else:
                    # Gap detected, save current range and start new one
                    ranges.append((current_start, current_end))
                    current_start = pos
                    current_end = pos
            
            # Add final range
            ranges.append((current_start, current_end))
            strand_ranges[chain] = ranges
        
        return strand_ranges

    def _format_source_label(self, metadata: Dict) -> str:
        """Format a display label combining _source_label and _also_found_in."""
        src = metadata.get('_source_label', '')
        also = metadata.get('_also_found_in', [])
        if also:
            return ' + '.join([src] + also)
        return src

    def _print_motif_instance_table(self, motif_type: str, motif_details: List[Dict]) -> None:
        """
        Print detailed instance table for a motif type.
        Uses official chainbreak metadata to identify separate strands.
        
        Args:
            motif_type: Motif type
            motif_details: List of motif instance details
        """
        # Data is already sorted at storage time (_load_motif_type)
        # so motif_details is in ascending order by minimum residue number
        
        print("\n" + "=" * 100)
        print(f"  {motif_type} MOTIF INSTANCES")
        print("=" * 100)
        print(f"  Total Instances: {len(motif_details)}")
        print("-" * 100)
        
        # Check if any instance has source metadata (combine mode)
        has_source = any(
            detail.get('metadata', {}).get('_source_label')
            for detail in motif_details
        )
        
        if has_source:
            # NUCLEOTIDES column hidden - nucleotide info available via rmv_show <TYPE> <NO>
            print(f"  {'NO.':<6} {'CHAIN':<10} {'RESIDUE RANGES':<50} {'SOURCE':<30}")
        else:
            print(f"  {'NO.':<6} {'CHAIN':<10} {'RESIDUE RANGES':<50}")
        print("-" * 100)
        
        for idx, detail in enumerate(motif_details, 1):
            residues = detail.get('residues', [])
            metadata = detail.get('metadata', {})
            
            if not residues:
                if has_source:
                    src = self._format_source_label(metadata)
                    print(f"  {idx:<6} {'-':<10} {'-':<50} {src:<30}")
                else:
                    print(f"  {idx:<6} {'-':<10} {'-':<50}")
                continue
            
            # Build position map: (chain, pos) -> nucleotide
            pos_to_nuc = {}
            for res in residues:
                if isinstance(res, tuple) and len(res) >= 3:
                    nucleotide, resi, chain = res[0], res[1], res[2]
                    pos_to_nuc[(chain, resi)] = nucleotide
            
            # Get unique chains
            chains = sorted(set(res[2] for res in residues if isinstance(res, tuple) and len(res) >= 3))
            chains_str = ', '.join(chains)
            
            # Residue ranges - prefer specific region metadata
            residue_range = None
            
            # First try: Check for RNAMotifScan 'regions' metadata
            if metadata and 'regions' in metadata and metadata['regions']:
                regions = metadata['regions']
                if isinstance(regions, list) and len(regions) > 0:
                    range_parts = []
                    for r in regions:
                        if isinstance(r, tuple) and len(r) >= 3:
                            chain, start, end = r[0], r[1], r[2]
                            range_parts.append(f"{chain}:{start}-{end}")
                    if range_parts:
                        residue_range = ', '.join(range_parts)
            
            # Second try: Check for RMSX 'aligned_regions' metadata
            elif metadata and 'aligned_regions' in metadata and metadata['aligned_regions']:
                aligned_regions = metadata['aligned_regions']
                if isinstance(aligned_regions, list) and len(aligned_regions) > 0:
                    range_parts = []
                    for r in aligned_regions:
                        if isinstance(r, tuple) and len(r) >= 2:
                            range_parts.append(f"{r[0]}-{r[1]}")
                    if range_parts:
                        residue_range = ', '.join(range_parts)
            
            # Third try: Use official chainbreak metadata (BGSU format)
            elif metadata and 'chainbreak' in metadata:
                chainbreak = metadata['chainbreak']
                residue_range = self._get_ranges_from_chainbreak(residues, chainbreak, metadata)
            
            # Fallback: Identify strands by detecting gaps in residue numbers
            if not residue_range:
                strand_ranges = self._identify_strands(residues)
                range_parts = []
                
                for chain in sorted(strand_ranges.keys()):
                    for start, end in strand_ranges[chain]:
                        range_parts.append(f"{chain}:{start}-{end}")
                
                residue_range = ', '.join(range_parts) if range_parts else '-'
            
            # Nucleotide string computation kept but not displayed in table
            # nucs_str = self._get_nucleotides_for_strands(residues, metadata, pos_to_nuc)
            
            if has_source:
                src = self._format_source_label(metadata)
                print(f"  {idx:<6} {chains_str:<10} {residue_range:<50} {src:<30}")
            else:
                print(f"  {idx:<6} {chains_str:<10} {residue_range:<50}")
        
        print("-" * 100)
        print("\n  To view a specific instance:")
        print(f"    rmv_show {motif_type} <NO>")
        print(f"    Example: rmv_show {motif_type} 1")
        print("=" * 100 + "\n")
    
    def show_motif_instance(self, motif_type: str, instance_no: int,
                            filter_pdb: str = '', filter_suffix: str = '',
                            padding: int = 0) -> bool:
        """
        Show only a specific instance of a motif type highlighted, with full structure in gray.
        Creates a separate PyMOL object for the instance.
        
        Workflow:
        1. Hide all separate motif objects
        2. Show full structure in gray80
        3. Create a separate object for the instance and highlight it
        
        Args:
            motif_type: Motif type (e.g., 'GNRA')
            instance_no: Instance number (1-indexed)
            filter_pdb: When set, only show instances matching this PDB ID
            filter_suffix: When set, only show instances matching this source suffix
            
        Returns:
            True if successful
        """
        motif_type = motif_type.upper().strip()
        loaded_motifs = self.motif_loader.get_loaded_motifs()
        
        if motif_type not in loaded_motifs:
            self.logger.error(f"Motif type '{motif_type}' not loaded")
            return False
        
        info = loaded_motifs[motif_type]
        motif_details = info.get('motif_details', [])
        structure_name = info.get('structure_name')
        # Use filter_suffix for naming when provided — correct for combine mode
        # where info-level source_suffix may be stale from a previous single-source load
        source_suffix = filter_suffix if filter_suffix else info.get('source_suffix', '')
        
        # Apply source filter when requested (show only current source data)
        if filter_pdb:
            motif_details = [
                d for d in motif_details
                if d.get('_pdb_id', info.get('pdb_id', '')) == filter_pdb
                and d.get('_source_suffix', info.get('source_suffix', '')) == filter_suffix
            ]
        
        if instance_no < 1 or instance_no > len(motif_details):
            self.logger.error(f"Instance {instance_no} not found. Valid range: 1-{len(motif_details)}")
            return False
        
        if not structure_name:
            self.logger.error("No structure name found")
            return False
        
        # Hide ALL separate motif objects
        for mt, mt_info in loaded_motifs.items():
            obj_name = mt_info.get('object_name')
            if obj_name:
                self.cmd.disable(obj_name)
        
        # Hide any previously created instance objects
        for obj in self.cmd.get_object_list():
            for mt in loaded_motifs.keys():
                mt_sanitized = sanitize_pymol_name(mt)
                sfx = loaded_motifs[mt].get('source_suffix', '')
                mt_pdb_id = loaded_motifs[mt].get('pdb_id', '')
                pdb_tag = f"_{mt_pdb_id}" if mt_pdb_id else ""
                prefix = sanitize_pymol_name(f"{mt_sanitized}{pdb_tag}{sfx}")
                if obj.startswith(f"{prefix}_") and obj[len(prefix)+1:].isdigit():
                    self.cmd.disable(obj)
        
        # Show the full structure with uniform representation in gray80
        self.cmd.enable(structure_name)
        self.cmd.show('cartoon', f"model {structure_name} and polymer.nucleic")
        self.cmd.set('cartoon_nucleic_acid_mode', 4, f"model {structure_name}")
        self.cmd.set('cartoon_tube_radius', 0.4, f"model {structure_name}")
        self.cmd.color('gray80', f"model {structure_name} and polymer.nucleic")
        
        # Color the instance residues WITHIN the main structure
        detail = motif_details[instance_no - 1]
        residues = detail.get('residues', [])
        inst_pdb_id = info.get('pdb_id', self.structure_loader.get_current_pdb_id() or '')
        inst_pdb_tag = f"_{inst_pdb_id}" if inst_pdb_id else ""
        # Compute both object name variants (with and without _P)
        base_obj = sanitize_pymol_name(f"{motif_type}_{instance_no}{inst_pdb_tag}{source_suffix}")
        pad_obj = base_obj + "_P"
        instance_obj = pad_obj if padding > 0 else base_obj
        
        if residues:
            from .utils.parser import SelectionParser
            
            # Build chain-residue mapping
            chain_residues = {}
            for res in residues:
                if isinstance(res, tuple) and len(res) >= 3:
                    nucleotide, resi, chain = res[0], res[1], res[2]
                    if chain not in chain_residues:
                        chain_residues[chain] = []
                    chain_residues[chain].append(resi)
            
            # Create selection for this instance
            selections = []
            for chain, resi_list in chain_residues.items():
                if padding > 0:
                    expanded = set()
                    for r in resi_list:
                        for offset in range(-padding, padding + 1):
                            expanded.add(r + offset)
                    resi_list = sorted(expanded)
                sel = SelectionParser.create_selection_string(chain, sorted(resi_list))
                if sel:
                    selections.append(f"({sel})")
            
            if selections:
                combined_sel = " or ".join(selections)
                instance_sel = f"(model {structure_name}) and ({combined_sel})"
                
                # Color the instance residues WITHIN the main structure (no overlap)
                colors.set_motif_color_in_pymol(self.cmd, instance_sel, motif_type)
                
                # Clean up conflicting objects to avoid visual overlap
                existing_objects = self.cmd.get_object_list()
                # Delete the OTHER variant (base vs padded) so both aren't visible
                other_obj = pad_obj if padding == 0 else base_obj
                if other_obj in existing_objects:
                    self.cmd.delete(other_obj)
                # Always delete and recreate target so padding value changes
                # are reflected (e.g., padding=10 → padding=2)
                if instance_obj in existing_objects:
                    self.cmd.delete(instance_obj)
                try:
                    self.cmd.create(instance_obj, instance_sel)
                    self.cmd.set('cartoon_nucleic_acid_mode', 4, instance_obj, quiet=1)
                    self.cmd.set('cartoon_tube_radius', 0.4, instance_obj, quiet=1)
                    colors.set_motif_color_in_pymol(self.cmd, instance_obj, motif_type)
                    self.cmd.enable(instance_obj)
                except Exception as e:
                    self.logger.debug(f"Could not create instance object: {e}")
                
                # Zoom to the instance (using the selection, not the object)
                self.cmd.zoom(instance_sel, 5)
        
        # Print instance details
        self._print_single_instance_info(motif_type, instance_no, detail,
                                         source_suffix=source_suffix,
                                         pdb_id=inst_pdb_id)
        
        # Print follow-up suggestions
        print("  Next steps:")
        if instance_no > 1:
            print(f"    rmv_show {motif_type} {instance_no-1}             View previous instance")
        if instance_no < len(motif_details):
            print(f"    rmv_show {motif_type} {instance_no+1}             View next instance")
        print(f"    rmv_show {motif_type}               Show all {motif_type} instances")
        print(f"    rmv_super {motif_type}              Superimpose all {motif_type} instances")
        print(f"    rmv_save {motif_type} {instance_no}              Save image of this instance")
        print(f"    rmv_show ALL                 Show all motif types")
        print()
        
        return True
    
    def _print_single_instance_info(self, motif_type: str, instance_no: int, 
                                     detail: Dict, source_suffix: str = '',
                                     pdb_id: str = '',
                                     show_object: bool = True) -> None:
        """Print information about a single motif instance."""
        residues = detail.get('residues', [])
        instance_id = detail.get('instance_id', f'{motif_type}_{instance_no}')
        annotation = detail.get('annotation', '')
        metadata = detail.get('metadata', {})
        
        print("\n" + "=" * 50)
        print(f"  {motif_type} INSTANCE #{instance_no}")
        print("=" * 50)
        print(f"  Instance ID: {instance_id}")
        if annotation:
            print(f"  Annotation: {annotation}")
        
        # Display numeric metadata if available (RMSX specific)
        if metadata:
            alignment_score = metadata.get('alignment_score')
            p_value = metadata.get('p_value')
            if alignment_score is not None or p_value is not None:
                score_str = f"{alignment_score:.1f}" if alignment_score is not None else "N/A"
                pval_str = f"{p_value:.6f}" if p_value is not None else "N/A"
                print(f"  Alignment Score: {score_str} | P-value: {pval_str}")
        
        print(f"  Residues: {len(residues)}")
        print("-" * 50)
        
        # List all residues (nucleotide column hidden)
        print(f"  {'CHAIN':<8} {'RESI':<8}")
        print("-" * 50)
        
        for res in residues:
            if isinstance(res, tuple) and len(res) >= 3:
                nucleotide, resi, chain = res[0], res[1], res[2]
                print(f"  {chain:<8} {resi:<8}")
        
        print("=" * 50)
        if show_object:
            pdb_tag = f"_{pdb_id}" if pdb_id else ""
            obj_display = sanitize_pymol_name(f"{motif_type}_{instance_no}{pdb_tag}{source_suffix}")
            print(f"  Object: {obj_display}")
            print("=" * 50)
        print()
    
    def show_all_motifs(self, filter_pdb: str = '', filter_suffix: str = '') -> None:
        """
        Show ALL loaded motif types — creates PyMOL objects for each type
        (like calling rmv_show <TYPE> for every loaded motif type at once).
        
        Workflow:
        1. Create PyMOL objects for every loaded motif type (if not already created)
        2. Hide individual instance objects (avoid clutter)
        3. Show full PDB structure uniformly in gray80
        4. Enable all motif-type objects and color their residues

        Args:
            filter_pdb: When set, only show instances matching this PDB ID
            filter_suffix: When set, only show instances matching this source suffix
        """
        loaded_motifs = self.motif_loader.get_loaded_motifs()
        
        if not loaded_motifs:
            self.logger.error("No motifs loaded. Use 'rmv_load_motif' first.")
            return
        
        # Get structure name
        structure_name = None
        for info in loaded_motifs.values():
            if info.get('structure_name'):
                structure_name = info.get('structure_name')
                break
        
        if not structure_name:
            self.logger.error("No structure loaded")
            return
        
        # Step 1: Create PyMOL objects for each motif type (like show_motif_type does)
        for motif_type, info in loaded_motifs.items():
            obj_name = info.get('object_name')
            motif_list = info.get('motifs', [])
            source_suffix = info.get('source_suffix', '')
            pdb_id = info.get('pdb_id', self.structure_loader.get_current_pdb_id() or '')
            
            if not obj_name and motif_list:
                obj_name = self.motif_loader.selector.create_motif_class_object(
                    structure_name,
                    motif_type,
                    motif_list,
                    source_suffix=source_suffix,
                    pdb_id=pdb_id,
                )
                if obj_name:
                    colors.set_motif_color_in_pymol(self.cmd, obj_name, motif_type)
                    info['object_name'] = obj_name
                    loaded_motifs[motif_type] = info
                    self.logger.debug(f"Created PyMOL object: {obj_name}")
        
        # Step 2: Hide individual instance objects (avoid clutter)
        for obj in self.cmd.get_object_list():
            for mt, mt_info in loaded_motifs.items():
                mt_sanitized = sanitize_pymol_name(mt)
                suffix = mt_info.get('source_suffix', '')
                mt_pdb_id = mt_info.get('pdb_id', '')
                pdb_tag = f"_{mt_pdb_id}" if mt_pdb_id else ""
                prefix = sanitize_pymol_name(f"{mt_sanitized}{pdb_tag}{suffix}")
                if obj.startswith(f"{prefix}_") and obj[len(prefix)+1:].isdigit():
                    self.cmd.disable(obj)
        
        # Step 3: Show the full structure with uniform representation
        self.cmd.enable(structure_name)
        self.cmd.show('cartoon', f"model {structure_name} and polymer.nucleic")
        self.cmd.set('cartoon_nucleic_acid_mode', 4, f"model {structure_name}")
        self.cmd.set('cartoon_tube_radius', 0.4, f"model {structure_name}")
        
        # Color the ENTIRE structure gray80 first
        self.cmd.color('gray80', f"model {structure_name} and polymer.nucleic")
        
        # Step 4: Enable all motif-type objects and color residues on the structure
        total_instances = 0
        from .utils.parser import SelectionParser
        for motif_type, info in loaded_motifs.items():
            # Enable the motif-type PyMOL object
            obj_name = info.get('object_name')
            if obj_name:
                self.cmd.enable(obj_name)
            
            motif_details = info.get('motif_details', [])
            # Apply source filter when requested
            if filter_pdb:
                motif_details = [
                    d for d in motif_details
                    if d.get('_pdb_id', info.get('pdb_id', '')) == filter_pdb
                    and d.get('_source_suffix', info.get('source_suffix', '')) == filter_suffix
                ]
            if not motif_details:
                continue
            total_instances += len(motif_details)
            
            # Color each instance's residues on the main structure
            for detail in motif_details:
                residues = detail.get('residues', [])
                if not residues:
                    continue
                
                chain_residues = {}
                for res in residues:
                    if isinstance(res, tuple) and len(res) >= 3:
                        nucleotide, resi, chain = res[0], res[1], res[2]
                        if chain not in chain_residues:
                            chain_residues[chain] = []
                        chain_residues[chain].append(resi)
                
                selections = []
                for chain, resi_list in chain_residues.items():
                    sel = SelectionParser.create_selection_string(chain, sorted(resi_list))
                    if sel:
                        selections.append(f"({sel})")
                
                if selections:
                    combined_sel = " or ".join(selections)
                    instance_sel = f"(model {structure_name}) and ({combined_sel})"
                    try:
                        self.cmd.show('cartoon', instance_sel)
                        colors.set_motif_color_in_pymol(self.cmd, instance_sel, motif_type)
                    except Exception as e:
                        self.logger.debug(f"Could not color instance {motif_type}: {e}")
        
        self.logger.success(f"Showing all {len(loaded_motifs)} motif types ({total_instances} total instances)")
        
        # Print summary of created objects
        print(f"\n  PyMOL objects created for {len(loaded_motifs)} motif types:")
        for motif_type, info in sorted(loaded_motifs.items()):
            obj_name = info.get('object_name', 'N/A')
            count = len(info.get('motif_details', []))
            print(f"    {obj_name:<25} {motif_type} ({count} instances)")
        
        # Print follow-up suggestions
        print("\n  Next steps:")
        if loaded_motifs:
            first_motif = next(iter(sorted(loaded_motifs.keys())), None)
            if first_motif:
                print(f"    rmv_show {first_motif:<20}  Focus on specific motif type")
        print(f"    rmv_summary              View motif summary table")
        print(f"    rmv_save ALL             Save all motif images")
        print(f"    rmv_save HL              Save specific motif type images")
        print(f"    rmv_fetch <PDB_ID>       Load a different structure")
        print()
    
    def save_all_motif_images(self, representation: str = 'cartoon') -> bool:
        """
        Save images of all loaded motif instances.
        
        Creates folder structure: plugin_dir/motif_images/pdb_id/MOTIF_TYPE/instance_*_info.png
        
        Args:
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        
        Returns:
            True if successful, False otherwise
        """
        try:
            from .image_saver import MotifImageSaver
            
            loaded_motifs = self.motif_loader.get_loaded_motifs()
            
            if not loaded_motifs:
                self.logger.warning("No motifs loaded to save")
                return False
            
            # Extract structure name and PDB ID from loaded motifs
            structure_name = None
            pdb_id = None
            
            for motif_type, info in loaded_motifs.items():
                if info.get('structure_name'):
                    structure_name = info.get('structure_name')
                    break
            
            # If not found in motifs, try structure_loader
            if not structure_name:
                structure_name = self.structure_loader.get_current_structure()
            
            # Get PDB ID
            if not pdb_id:
                pdb_id = self.structure_loader.get_current_pdb_id()
            
            if not structure_name or not pdb_id:
                self.logger.error("No structure loaded or motif data incomplete")
                return False
            
            saver = MotifImageSaver(self.cmd)
            stats = saver.save_all_motifs(loaded_motifs, structure_name, pdb_id,
                                         representation=representation)
            
            return True
        
        except Exception as e:
            self.logger.error(f"Failed to save motif images: {e}")
            return False
    
    def save_motif_type_images(self, motif_type: str, representation: str = 'cartoon') -> bool:
        """
        Save images of all instances of a specific motif type.
        
        Creates folder structure: plugin_dir/motif_images/pdb_id/MOTIF_TYPE/instance_*_info.png
        
        Args:
            motif_type: Motif type to save (e.g., 'HL', 'IL')
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        
        Returns:
            True if successful, False otherwise
        """
        try:
            from .image_saver import MotifImageSaver
            
            loaded_motifs = self.motif_loader.get_loaded_motifs()
            
            if not loaded_motifs:
                self.logger.warning("No motifs loaded to save")
                return False
            
            # Extract structure name and PDB ID from loaded motifs
            structure_name = None
            pdb_id = None
            
            for motif_type_check, info in loaded_motifs.items():
                if info.get('structure_name'):
                    structure_name = info.get('structure_name')
                    break
            
            # If not found in motifs, try structure_loader
            if not structure_name:
                structure_name = self.structure_loader.get_current_structure()
            
            # Get PDB ID
            if not pdb_id:
                pdb_id = self.structure_loader.get_current_pdb_id()
            
            if not structure_name or not pdb_id:
                self.logger.error("No structure loaded or motif data incomplete")
                return False
            
            saver = MotifImageSaver(self.cmd)
            stats = saver.save_motif_type_images(loaded_motifs, motif_type, 
                                               structure_name, pdb_id,
                                               representation=representation)
            
            if stats['saved'] == 0 and stats['failed'] == 0:
                return False
            
            return True
        
        except Exception as e:
            self.logger.error(f"Failed to save {motif_type} images: {e}")
            return False    
    def save_motif_instance_by_id(self, motif_type: str, instance_id: int, 
                                  representation: str = 'cartoon') -> bool:
        """
        Save image of a specific motif instance by ID.
        
        Creates file: plugin_dir/motif_images/pdb_id/MOTIF_TYPE/instance_ID_info.png
        
        Args:
            motif_type: Motif type (e.g., 'HL', 'IL')
            instance_id: Instance number (1-based, as shown in rmv_summary)
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        
        Returns:
            True if successful, False otherwise
        """
        try:
            from .image_saver import MotifImageSaver
            
            loaded_motifs = self.motif_loader.get_loaded_motifs()
            
            if not loaded_motifs:
                self.logger.warning("No motifs loaded to save")
                return False
            
            if motif_type not in loaded_motifs:
                self.logger.error(f"Motif type '{motif_type}' not found")
                return False
            
            # Get instances for this motif type
            motif_data = loaded_motifs[motif_type]
            instances = motif_data.get('motif_details', [])
            
            if not instances:
                self.logger.error(f"No instances found for motif type '{motif_type}'")
                return False
            
            # Validate instance ID
            if instance_id < 1 or instance_id > len(instances):
                self.logger.error(f"Instance ID {instance_id} out of range (1-{len(instances)})")
                return False
            
            # Extract structure name and PDB ID
            structure_name = motif_data.get('structure_name')
            pdb_id = None
            
            # If not found in motif data, try structure_loader
            if not structure_name:
                structure_name = self.structure_loader.get_current_structure()
            
            # Get PDB ID
            if not pdb_id:
                pdb_id = self.structure_loader.get_current_pdb_id()
            
            if not structure_name or not pdb_id:
                self.logger.error("No structure loaded or motif data incomplete")
                return False
            
            # Save the specific instance
            saver = MotifImageSaver(self.cmd)
            success = saver.save_motif_instance(loaded_motifs, motif_type, instance_id,
                                              structure_name, pdb_id,
                                              representation=representation)
            
            return success
        
        except Exception as e:
            self.logger.error(f"Failed to save {motif_type} instance #{instance_id}: {e}")
            return False
    
    def save_current_view(self, filename: str) -> bool:
        """
        Save the current PyMOL view to PNG (like PyMOL's png command).
        
        Saves the current view exactly as displayed in the PyMOL window,
        without any modifications or automatic zooming.
        
        Args:
            filename: Output filename (e.g., 'my_structure.png')
        
        Returns:
            True if successful, False otherwise
        """
        try:
            from .image_saver import MotifImageSaver
            
            saver = MotifImageSaver(self.cmd)
            success = saver.save_current_view(filename)
            
            return success
        
        except Exception as e:
            self.logger.error(f"Failed to save current view: {e}")
            return False

    # ------------------------------------------------------------------
    # Structure export (mmCIF) methods
    # ------------------------------------------------------------------

    def export_all_motif_structures(self) -> bool:
        """Export all loaded motif instances as mmCIF files.

        Extracts ORIGINAL coordinates from the on-disk CIF file.

        Returns:
            True if any instances were exported successfully.
        """
        try:
            from .structure_exporter import MotifStructureExporter

            loaded_motifs = self.motif_loader.get_loaded_motifs()
            if not loaded_motifs:
                self.logger.warning("No motifs loaded to export")
                return False

            pdb_id = self.structure_loader.get_current_pdb_id()
            if not pdb_id:
                self.logger.error("No structure loaded")
                return False

            exporter = MotifStructureExporter(self.cmd)
            stats = exporter.export_all_motifs(loaded_motifs, pdb_id)
            return stats.get('total_saved', 0) > 0

        except Exception as e:
            self.logger.error(f"Failed to export motif structures: {e}")
            return False

    def export_motif_type_structures(self, motif_type: str) -> bool:
        """Export all instances of a specific motif type as mmCIF files.

        Args:
            motif_type: Motif type to export (e.g., 'HL', 'IL').

        Returns:
            True if any instances were exported successfully.
        """
        try:
            from .structure_exporter import MotifStructureExporter

            loaded_motifs = self.motif_loader.get_loaded_motifs()
            if not loaded_motifs:
                self.logger.warning("No motifs loaded to export")
                return False

            pdb_id = self.structure_loader.get_current_pdb_id()
            if not pdb_id:
                self.logger.error("No structure loaded")
                return False

            exporter = MotifStructureExporter(self.cmd)
            stats = exporter.export_motif_type(loaded_motifs, motif_type, pdb_id)
            return stats.get('saved', 0) > 0

        except Exception as e:
            self.logger.error(f"Failed to export {motif_type} structures: {e}")
            return False

    def export_motif_instance_structure(self, motif_type: str,
                                        instance_id: int) -> bool:
        """Export a specific motif instance as mmCIF file.

        Args:
            motif_type: Motif type (e.g., 'HL', 'IL').
            instance_id: Instance number (1-based).

        Returns:
            True if exported successfully.
        """
        try:
            from .structure_exporter import MotifStructureExporter

            loaded_motifs = self.motif_loader.get_loaded_motifs()
            if not loaded_motifs:
                self.logger.warning("No motifs loaded to export")
                return False

            pdb_id = self.structure_loader.get_current_pdb_id()
            if not pdb_id:
                self.logger.error("No structure loaded")
                return False

            exporter = MotifStructureExporter(self.cmd)
            return exporter.export_motif_instance(
                loaded_motifs, motif_type, instance_id, pdb_id)

        except Exception as e:
            self.logger.error(f"Failed to export {motif_type} #{instance_id}: {e}")
            return False

    # ------------------------------------------------------------------
    # rmv_view helpers — zoom / select on the BASE structure (no objects)
    # ------------------------------------------------------------------

    def view_motif_type(self, motif_type: str,
                        filter_pdb: str = '', filter_suffix: str = '') -> bool:
        """Zoom to all instances of a motif type on the base structure.

        No new PyMOL objects are created.  The structure is colored gray80
        and only the matching residues are highlighted with the motif color.
        """
        motif_type = motif_type.upper().strip()
        loaded_motifs = self.motif_loader.get_loaded_motifs()

        if motif_type not in loaded_motifs:
            self.logger.error(f"Motif type '{motif_type}' not loaded")
            self.logger.info(f"Available: {', '.join(loaded_motifs.keys())}")
            return False

        info = loaded_motifs[motif_type]
        structure_name = info.get('structure_name')
        motif_details = info.get('motif_details', [])

        if filter_pdb:
            motif_details = [
                d for d in motif_details
                if d.get('_pdb_id', info.get('pdb_id', '')) == filter_pdb
                and d.get('_source_suffix', info.get('source_suffix', '')) == filter_suffix
            ]
            if not motif_details:
                self.logger.error(f"No {motif_type} instances for current source")
                return False

        if not structure_name:
            self.logger.error("No structure name found")
            return False

        from .utils.parser import SelectionParser

        # Show structure in gray
        self.cmd.enable(structure_name)
        self.cmd.show('cartoon', f"model {structure_name} and polymer.nucleic")
        self.cmd.set('cartoon_nucleic_acid_mode', 4, f"model {structure_name}")
        self.cmd.set('cartoon_tube_radius', 0.4, f"model {structure_name}")
        self.cmd.color('gray80', f"model {structure_name} and polymer.nucleic")

        # Collect all residue selections for zoom target
        all_selections = []
        for detail in motif_details:
            residues = detail.get('residues', [])
            if not residues:
                continue
            chain_residues = {}
            for res in residues:
                if isinstance(res, tuple) and len(res) >= 3:
                    chain_residues.setdefault(res[2], []).append(res[1])
            for chain, resi_list in chain_residues.items():
                sel = SelectionParser.create_selection_string(chain, sorted(resi_list))
                if sel:
                    all_selections.append(f"({sel})")
                    # Color on structure
                    instance_sel = f"(model {structure_name}) and ({sel})"
                    colors.set_motif_color_in_pymol(self.cmd, instance_sel, motif_type)

        if all_selections:
            combined = " or ".join(all_selections)
            zoom_sel = f"(model {structure_name}) and ({combined})"
            self.cmd.zoom(zoom_sel, 8)

        self._print_motif_instance_table(motif_type, motif_details)
        self.logger.success(f"Viewing {len(motif_details)} {motif_type} instances")

        print("  Next steps:")
        print(f"    rmv_view {motif_type} <NO>         Zoom to specific instance")
        print(f"    rmv_show {motif_type}              Create objects & render")
        print()
        return True

    def view_motif_instance(self, motif_type: str, instance_no: int,
                            filter_pdb: str = '', filter_suffix: str = '') -> bool:
        """Zoom to a specific motif instance and create a named SELECTION (not object).

        The base structure stays as the only object; a PyMOL *selection*
        (sele_<TYPE>_<NO>) is created so the user can inspect it.
        """
        motif_type = motif_type.upper().strip()
        loaded_motifs = self.motif_loader.get_loaded_motifs()

        if motif_type not in loaded_motifs:
            self.logger.error(f"Motif type '{motif_type}' not loaded")
            return False

        info = loaded_motifs[motif_type]
        motif_details = info.get('motif_details', [])
        structure_name = info.get('structure_name')

        if filter_pdb:
            motif_details = [
                d for d in motif_details
                if d.get('_pdb_id', info.get('pdb_id', '')) == filter_pdb
                and d.get('_source_suffix', info.get('source_suffix', '')) == filter_suffix
            ]

        if instance_no < 1 or instance_no > len(motif_details):
            self.logger.error(f"Instance {instance_no} not found. Valid: 1-{len(motif_details)}")
            return False

        if not structure_name:
            self.logger.error("No structure name found")
            return False

        from .utils.parser import SelectionParser

        # Show structure in gray
        self.cmd.enable(structure_name)
        self.cmd.show('cartoon', f"model {structure_name} and polymer.nucleic")
        self.cmd.set('cartoon_nucleic_acid_mode', 4, f"model {structure_name}")
        self.cmd.set('cartoon_tube_radius', 0.4, f"model {structure_name}")
        self.cmd.color('gray80', f"model {structure_name} and polymer.nucleic")

        detail = motif_details[instance_no - 1]
        residues = detail.get('residues', [])

        if residues:
            chain_residues = {}
            for res in residues:
                if isinstance(res, tuple) and len(res) >= 3:
                    chain_residues.setdefault(res[2], []).append(res[1])

            selections = []
            for chain, resi_list in chain_residues.items():
                sel = SelectionParser.create_selection_string(chain, sorted(resi_list))
                if sel:
                    selections.append(f"({sel})")

            if selections:
                combined_sel = " or ".join(selections)
                instance_sel = f"(model {structure_name}) and ({combined_sel})"
                colors.set_motif_color_in_pymol(self.cmd, instance_sel, motif_type)

                # Create a named SELECTION (not object)
                sel_name = sanitize_pymol_name(f"sele_{motif_type}_{instance_no}")
                self.cmd.select(sel_name, instance_sel)

                self.cmd.zoom(instance_sel, 5)

        # Print instance details
        pdb_id = info.get('pdb_id', self.structure_loader.get_current_pdb_id() or '')
        source_suffix = info.get('source_suffix', '')
        self._print_single_instance_info(motif_type, instance_no, detail,
                                         source_suffix=source_suffix,
                                         pdb_id=pdb_id, show_object=False)

        sel_display = sanitize_pymol_name(f"sele_{motif_type}_{instance_no}")
        print(f"  Selection: {sel_display}")
        print()
        print("  Next steps:")
        if instance_no > 1:
            print(f"    rmv_view {motif_type} {instance_no-1}             Previous instance")
        if instance_no < len(motif_details):
            print(f"    rmv_view {motif_type} {instance_no+1}             Next instance")
        print(f"    rmv_show {motif_type} {instance_no}              Create object & render")
        print()
        return True