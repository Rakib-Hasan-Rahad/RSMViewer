"""
RSMViewer - Selection Module
Handles creation and management of PyMOL selections for motifs.
"""

import re
from .parser import SelectionParser, validate_motif_data
from .logger import get_logger


def sanitize_pymol_name(name):
    """Sanitize a string for use as a PyMOL object/selection name.
    
    Replaces spaces, parentheses, and other invalid characters with underscores.
    Strips leading/trailing underscores and collapses multiple underscores.
    
    Args:
        name (str): Raw name string
    
    Returns:
        str: Sanitized name safe for PyMOL
    """
    # Replace invalid characters with underscore
    clean = re.sub(r'[^A-Za-z0-9_]', '_', name)
    # Collapse multiple underscores
    clean = re.sub(r'_+', '_', clean)
    # Strip leading/trailing underscores
    clean = clean.strip('_')
    return clean


class MotifSelector:
    """Manages PyMOL selections for RNA motifs."""
    
    def __init__(self, cmd):
        """
        Initialize selector.
        
        Args:
            cmd: PyMOL cmd module
        """
        self.cmd = cmd
        self.logger = get_logger()
    
    def create_motif_object(self, structure_name, motif_type, motif_id, chain, residues, source_suffix=''):
        """
        Create a PyMOL object for a single motif.
        
        Includes Layer 3 protection: after creating the selection, verifies
        that atoms were actually found. If not, tries the 'segi' field as
        a fallback (which stores label_asym_id when cif_use_auth=1).
        
        Args:
            structure_name (str): Name of the loaded structure
            motif_type (str): Type of motif (e.g., 'KTURN')
            motif_id (str): Unique ID for this motif instance
            chain (str): Chain identifier
            residues (list): List of residue numbers
            source_suffix (str): Source suffix for object naming (e.g., '_S3')
        
        Returns:
            str: Name of created PyMOL object
        """
        if not validate_motif_data({'chain': chain, 'residues': residues, 'motif_id': motif_id}):
            return None
        
        try:
            # Create object name (sanitized for PyMOL)
            obj_name = sanitize_pymol_name(f"{motif_type}_{motif_id}_{chain}{source_suffix}")
            
            # Create selection with Layer 2 validation
            selection = SelectionParser.create_selection_string(
                chain, residues, structure_name=structure_name
            )
            if not selection:
                return None
            
            # Create object in PyMOL
            full_selection = f"(model {structure_name}) and {selection}"
            self.cmd.create(obj_name, full_selection)
            
            # Layer 3: Verify atoms were found, try segi fallback if not
            atom_count = self.cmd.count_atoms(obj_name)
            if atom_count == 0:
                self.logger.debug(
                    f"No atoms found for {motif_type} using chain '{chain}', "
                    f"trying segi fallback..."
                )
                # Try segi fallback (label_asym_id)
                segi_selection = SelectionParser.create_selection_string(
                    chain, residues, structure_name=structure_name, use_segi=True
                )
                if segi_selection:
                    segi_full = f"(model {structure_name}) and {segi_selection}"
                    segi_count = self.cmd.count_atoms(segi_full)
                    if segi_count > 0:
                        self.cmd.delete(obj_name)
                        self.cmd.create(obj_name, segi_full)
                        self.logger.info(
                            f"Recovered {segi_count} atoms for {motif_type} "
                            f"using segi fallback (chain '{chain}')"
                        )
                    else:
                        self.logger.warning(
                            f"Chain '{chain}' not found in '{structure_name}' "
                            f"via chain or segi. Available chains: "
                            f"{self.cmd.get_chains(structure_name)}"
                        )
            
            return obj_name
        except Exception as e:
            self.logger.error(f"Failed to create motif object: {e}")
            return None
    
    def create_motif_class_object(self, structure_name, motif_type, motif_list, source_suffix='', pdb_id=''):
        """
        Create a combined PyMOL object for all motifs of a type.
        
        Includes Layer 3 protection: after creating the combined object,
        verifies atoms were found. If not, retries with segi fallback.
        
        Args:
            structure_name (str): Name of the loaded structure
            motif_type (str): Type of motif (e.g., 'KTURN')
            motif_list (list): List of motif dictionaries with keys: chain, residues, motif_id
            source_suffix (str): Source suffix for object naming (e.g., '_S3')
            pdb_id (str): PDB ID to include in object name (e.g., '1S72')
        
        Returns:
            str: Name of created PyMOL object (e.g., 'KTURN_ALL_1S72_S3')
        """
        if not motif_list:
            return None
        
        try:
            pdb_tag = f"_{pdb_id.upper()}" if pdb_id else ""
            obj_name = sanitize_pymol_name(f"{motif_type}_ALL{pdb_tag}{source_suffix}")
            
            # Collect all selection strings with Layer 2 validation
            selections = []
            for motif in motif_list:
                if not validate_motif_data(motif):
                    self.logger.warning(f"Skipping invalid motif: {motif}")
                    continue
                
                chain = motif.get('chain')
                residues = motif.get('residues')
                
                selection = SelectionParser.create_selection_string(
                    chain, residues, structure_name=structure_name
                )
                if selection:
                    selections.append((selection, chain, residues))
            
            if not selections:
                self.logger.warning(f"No valid selections found for {motif_type}")
                return None
            
            # Combine all selections with OR
            combined_selection = " or ".join([f"({s})" for s, _, _ in selections])
            full_selection = f"(model {structure_name}) and ({combined_selection})"
            
            # Create combined object
            self.cmd.create(obj_name, full_selection)
            
            # Layer 3: Verify atoms were found, try segi fallback if not
            atom_count = self.cmd.count_atoms(obj_name)
            if atom_count == 0:
                self.logger.debug(
                    f"No atoms found for {motif_type}_ALL, trying segi fallback..."
                )
                segi_selections = []
                for _, chain, residues in selections:
                    segi_sel = SelectionParser.create_selection_string(
                        chain, residues, structure_name=structure_name, use_segi=True
                    )
                    if segi_sel:
                        segi_selections.append(segi_sel)
                
                if segi_selections:
                    segi_combined = " or ".join([f"({s})" for s in segi_selections])
                    segi_full = f"(model {structure_name}) and ({segi_combined})"
                    segi_count = self.cmd.count_atoms(segi_full)
                    if segi_count > 0:
                        self.cmd.delete(obj_name)
                        self.cmd.create(obj_name, segi_full)
                        self.logger.info(
                            f"Recovered {segi_count} atoms for {motif_type}_ALL "
                            f"using segi fallback"
                        )
                    else:
                        self.logger.warning(
                            f"No atoms found for {motif_type}_ALL via chain or segi. "
                            f"Available chains: {self.cmd.get_chains(structure_name)}"
                        )
            else:
                self.logger.info(f"Created motif object: {obj_name}")
            
            return obj_name
        except Exception as e:
            self.logger.error(f"Failed to create motif class object {motif_type}: {e}")
            return None
    
    def color_motif_residues(self, structure_name, motif_type, motif_list, color_rgb, source_suffix=''):
        """
        Color residues directly on the structure without creating overlapping objects.
        This avoids z-fighting/striping artifacts.
        
        Includes Layer 3 protection: verifies atoms were selected and tries
        segi fallback if primary chain selection fails.
        
        Args:
            structure_name (str): Name of the loaded structure
            motif_type (str): Type of motif (e.g., 'KTURN')
            motif_list (list): List of motif dictionaries with keys: chain, residues
            color_rgb (tuple): RGB color tuple (0-1 range)
            source_suffix (str): Source suffix for selection naming (e.g., '_S3')
        
        Returns:
            str: Selection name for the colored residues
        """
        if not motif_list:
            return None
        
        try:
            selection_name = sanitize_pymol_name(f"{motif_type}_sel{source_suffix}")
            
            # Collect all selection strings with Layer 2 validation
            selections = []
            for motif in motif_list:
                if not validate_motif_data(motif):
                    continue
                
                chain = motif.get('chain')
                residues = motif.get('residues')
                
                selection = SelectionParser.create_selection_string(
                    chain, residues, structure_name=structure_name
                )
                if selection:
                    selections.append((selection, chain, residues))
            
            if not selections:
                return None
            
            # Combine all selections
            combined_selection = " or ".join([f"({s})" for s, _, _ in selections])
            full_selection = f"(model {structure_name}) and ({combined_selection})"
            
            # Create a named selection (not a new object)
            self.cmd.select(selection_name, full_selection)
            
            # Layer 3: Verify atoms were selected, try segi fallback if not
            atom_count = self.cmd.count_atoms(selection_name)
            if atom_count == 0:
                self.logger.debug(
                    f"No atoms selected for {motif_type} coloring, trying segi fallback..."
                )
                segi_selections = []
                for _, chain, residues in selections:
                    segi_sel = SelectionParser.create_selection_string(
                        chain, residues, structure_name=structure_name, use_segi=True
                    )
                    if segi_sel:
                        segi_selections.append(segi_sel)
                
                if segi_selections:
                    segi_combined = " or ".join([f"({s})" for s in segi_selections])
                    segi_full = f"(model {structure_name}) and ({segi_combined})"
                    segi_count = self.cmd.count_atoms(segi_full)
                    if segi_count > 0:
                        self.cmd.select(selection_name, segi_full)
                        self.logger.info(
                            f"Recovered {segi_count} atoms for {motif_type} coloring "
                            f"using segi fallback"
                        )
            
            # Define and apply color
            color_name = f"motif_{motif_type}"
            self.cmd.set_color(color_name, color_rgb)
            self.cmd.color(color_name, selection_name)
            
            # Hide the selection indicator (the pink squares)
            self.cmd.deselect()
            
            self.logger.info(f"Colored {motif_type} residues directly on structure")
            
            return selection_name
        except Exception as e:
            self.logger.error(f"Failed to color motif residues {motif_type}: {e}")
            return None
    
    def toggle_object_visibility(self, obj_name, visible):
        """
        Toggle visibility of a PyMOL object.
        
        Args:
            obj_name (str): Name of object
            visible (bool): True to show, False to hide
        """
        try:
            if visible:
                self.cmd.show('sticks', obj_name)
                self.cmd.show('cartoon', obj_name)
            else:
                self.cmd.hide('everything', obj_name)
            self.logger.debug(f"Set visibility of {obj_name} to {visible}")
        except Exception as e:
            self.logger.error(f"Failed to toggle visibility of {obj_name}: {e}")
    
    def delete_object(self, obj_name):
        """
        Delete a PyMOL object.
        
        Args:
            obj_name (str): Name of object to delete
        """
        try:
            self.cmd.delete(obj_name)
            self.logger.debug(f"Deleted object: {obj_name}")
        except Exception as e:
            self.logger.error(f"Failed to delete object {obj_name}: {e}")
    
    def highlight_object(self, obj_name):
        """
        Highlight a motif object (increase stick radius, show surfaces, etc.).
        
        Args:
            obj_name (str): Name of object to highlight
        """
        try:
            self.cmd.show('sticks', obj_name)
            self.cmd.show('cartoon', obj_name)
            self.cmd.set('stick_radius', 0.3, obj_name)
            self.cmd.set('cartoon_fancy_helices', 1)
        except Exception as e:
            self.logger.error(f"Failed to highlight object {obj_name}: {e}")
    
    def get_all_motif_objects(self):
        """
        Get all motif-related objects in the scene.
        
        Returns:
            list: Names of all objects
        """
        try:
            all_objects = self.cmd.get_names('objects')
            motif_objects = [obj for obj in all_objects if '_ALL' in obj or '_K' in obj or '_A' in obj]
            return motif_objects
        except Exception as e:
            self.logger.error(f"Failed to get motif objects: {e}")
            return []
