"""
RSMViewer - Image Saver Module
Handles saving motif instance images to organized folder structures.

Features:
- Save all motif images organized by type and instance
- Save specific motif type images
- Create organized folder hierarchy: PDB_ID/MOTIF_TYPE/instance_*_info.png
- Automatic image naming with instance ID, chain, and residue info

Author: CBB LAB @Rakib Hasan Rahad
Version: 1.0.0
"""

import os
from pathlib import Path
from typing import Dict, List, Optional

from .utils import get_logger


class MotifImageSaver:
    """Handles saving motif instance images to organized folders."""
    
    def __init__(self, cmd):
        """
        Initialize the image saver.
        
        Args:
            cmd: PyMOL cmd module
        """
        self.cmd = cmd
        self.logger = get_logger()
        self.default_image_width = 800
        self.default_image_height = 600
        self.default_image_format = 'png'
    
    def create_folder_hierarchy(self, pdb_id: str, output_base_dir: Optional[str] = None) -> Path:
        """
        Create folder hierarchy for saving images.
        
        Structure: output_dir/PDB_ID/
        
        Args:
            pdb_id: PDB ID (used as folder name)
            output_base_dir: Base directory for output (defaults to plugin_dir/motif_images)
        
        Returns:
            Path to the PDB folder
        """
        if output_base_dir is None:
            # Use plugin directory by default
            plugin_dir = Path(__file__).parent.parent  # rsmviewer is a package
            output_base_dir = plugin_dir / "motif_images"
        else:
            output_base_dir = Path(output_base_dir)
        
        pdb_folder = output_base_dir / pdb_id.lower()
        pdb_folder.mkdir(parents=True, exist_ok=True)
        
        self.logger.info(f"Created PDB folder: {pdb_folder}")
        return pdb_folder
    
    def create_motif_type_folder(self, pdb_folder: Path, motif_type: str) -> Path:
        """
        Create subfolder for a specific motif type.
        
        Args:
            pdb_folder: Path to PDB folder
            motif_type: Motif type (e.g., 'HL', 'IL', 'GNRA')
        
        Returns:
            Path to the motif type folder
        """
        motif_folder = pdb_folder / motif_type.upper()
        motif_folder.mkdir(parents=True, exist_ok=True)
        
        return motif_folder
    
    def generate_instance_filename(self, instance_no: int, motif_details: Dict, 
                                  motif_type: str = '') -> str:
        """
        Generate filename for a motif instance based on ID and residue info.
        
        Format: <type>-<no>-<chain>-<residues>.<ext>
        Example: HL-1-0-55_64.png (HL motif, instance 1, chain 0, residues 55-64)
        
        Args:
            instance_no: Instance number (1-based, as shown in rmv_summary)
            motif_details: Dictionary with instance info (residues, annotation, etc.)
            motif_type: Motif type (e.g., 'HL', 'IL') - optional for backward compatibility
        
        Returns:
            Filename string
        """
        residues = motif_details.get('residues', [])
        chain_residues = {}
        
        # Extract chains and residues
        for res in residues:
            if isinstance(res, tuple) and len(res) >= 3:
                nuc, resi, chain = res[0], res[1], res[2]
                if chain not in chain_residues:
                    chain_residues[chain] = []
                chain_residues[chain].append(str(resi))
        
        # Build filename: TYPE-NO-CHAIN-RESIDUES
        # Handle multiple chains by joining with underscores
        chain_parts = []
        if chain_residues:
            for chain in sorted(chain_residues.keys()):
                resis = chain_residues[chain]
                
                # Condense consecutive residues
                condensed = self._condense_residues(resis)
                chain_parts.append(f"{chain}_{condensed}")
        
        chain_residue_str = "_".join(chain_parts) if chain_parts else "unknown"
        
        # Build the full filename: TYPE-NO-CHAIN-RESIDUE
        if motif_type:
            filename = f"{motif_type}-{instance_no}-{chain_residue_str}.{self.default_image_format}"
        else:
            # Fallback for backward compatibility (old format)
            filename = f"instance_{instance_no}_{chain_residue_str}.{self.default_image_format}"
        
        return filename
    
    def _condense_residues(self, residue_list: List[str]) -> str:
        """
        Condense list of residue numbers into ranges.
        
        Example: ['1', '2', '3', '5', '6'] -> '1-3_5-6'
        
        Args:
            residue_list: List of residue numbers as strings
        
        Returns:
            Condensed residue string
        """
        try:
            residues = sorted([int(r) for r in residue_list])
        except (ValueError, TypeError):
            return "_".join(residue_list[:3])  # Fallback for non-numeric residues
        
        if not residues:
            return ""
        
        ranges = []
        start = residues[0]
        end = residues[0]
        
        for resi in residues[1:]:
            if resi == end + 1:
                end = resi
            else:
                if start == end:
                    ranges.append(str(start))
                else:
                    ranges.append(f"{start}-{end}")
                start = resi
                end = resi
        
        if start == end:
            ranges.append(str(start))
        else:
            ranges.append(f"{start}-{end}")
        
        return "_".join(ranges)
    
    def _apply_representation(self, object_name: str, representation: str = 'cartoon') -> None:
        """
        Apply display representation to a PyMOL object.
        
        Args:
            object_name: Name of the PyMOL object
            representation: Type of representation to apply
                Available: 'cartoon', 'sticks', 'spheres', 'ribbon', 'lines', 'cartoon+sticks'
                Default: 'cartoon'
        """
        try:
            rep = representation.lower().strip()
            
            # Handle combined representations
            if '+' in rep:
                # e.g., 'cartoon+sticks'
                parts = rep.split('+')
                for part in parts:
                    self.cmd.show(part.strip(), object_name)
            else:
                # Single representation
                # Validate representation type
                valid_reps = ['cartoon', 'sticks', 'spheres', 'ribbon', 'lines', 'licorice', 'surface']
                if rep not in valid_reps:
                    self.logger.warning(f"Unknown representation '{rep}', using 'cartoon'")
                    rep = 'cartoon'
                
                self.cmd.show(rep, object_name)
        except Exception as e:
            self.logger.error(f"Failed to apply representation '{representation}': {e}")
            # Fallback to cartoon
            try:
                self.cmd.show('cartoon', object_name)
            except:
                pass
    
    def save_instance_image(self, motif_folder: Path, instance_no: int, 
                           motif_type: str, motif_details: Dict,
                           structure_name: str, representation: str = 'cartoon') -> bool:
        """
        Save image of a single motif instance with motif type coloring.
        
        Args:
            motif_folder: Path to motif type folder
            instance_no: Instance number (1-based)
            motif_type: Motif type (e.g., 'HL')
            motif_details: Dictionary with instance residues and info
            structure_name: Name of the PDB structure in PyMOL
            representation: Display representation ('cartoon', 'sticks', 'spheres', 'ribbon', 'lines', 'cartoon+sticks')
                          Default: 'cartoon'
        
        Returns:
            True if successful, False otherwise
        """
        try:
            filename = self.generate_instance_filename(instance_no, motif_details, motif_type)
            filepath = motif_folder / filename
            
            residues = motif_details.get('residues', [])
            if not residues:
                self.logger.warning(f"No residues for {motif_type} instance {instance_no}, skipping")
                return False
            
            # Build selection for this instance
            from .utils.parser import SelectionParser
            chain_residues = {}
            
            for res in residues:
                if isinstance(res, tuple) and len(res) >= 3:
                    nuc, resi, chain = res[0], res[1], res[2]
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
                self.logger.warning(f"Could not create selection for {motif_type} instance {instance_no}")
                return False
            
            combined_sel = " or ".join(selections)
            instance_sel = f"(model {structure_name}) and ({combined_sel})"
            
            # Temporary object name
            temp_motif_obj = f"_tmp_motif_{instance_no}"
            
            try:
                # HIDE original structure temporarily
                self.cmd.disable(structure_name)
                
                # Create temporary object with ONLY the motif residues
                temp_instance_sel = f"(model {structure_name}) and ({combined_sel})"
                self.cmd.create(temp_motif_obj, temp_instance_sel)
                self.cmd.enable(temp_motif_obj)
                
                # Set the color using the proper colors module function
                from . import colors
                colors.set_motif_color_in_pymol(self.cmd, temp_motif_obj, motif_type)
                
                # Apply the representation (default: cartoon)
                self._apply_representation(temp_motif_obj, representation)
                
                # Zoom to motif instance
                self.cmd.zoom(temp_motif_obj, 5)
                
                # Refresh and capture
                self.cmd.refresh()
                self.cmd.png(str(filepath), 
                            width=self.default_image_width, 
                            height=self.default_image_height)
                
            finally:
                # Delete temp object
                try:
                    self.cmd.delete(temp_motif_obj)
                except:
                    pass
                
                # RESTORE original structure visibility
                self.cmd.enable(structure_name)
                
                # Deselect
                try:
                    self.cmd.deselect()
                except:
                    pass
            
            self.logger.success(f"Saved: {filename}")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to save image for {motif_type} instance {instance_no}: {e}")
            return False
    
    def save_all_motifs(self, loaded_motifs: Dict, structure_name: str, 
                       pdb_id: str, output_base_dir: Optional[str] = None,
                       representation: str = 'cartoon') -> Dict:
        """
        Save images of all motif instances.
        
        Creates folder structure: plugin_dir/motif_images/pdb_id/MOTIF_TYPE/instance_*_info.png
        
        Args:
            loaded_motifs: Dictionary of loaded motifs from UnifiedMotifLoader
            structure_name: Name of structure in PyMOL
            pdb_id: PDB ID for folder naming
            output_base_dir: Base output directory (defaults to plugin_dir/motif_images)
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        
        Returns:
            Dictionary with save statistics
        """
        pdb_folder = self.create_folder_hierarchy(pdb_id, output_base_dir)
        
        stats = {
            'pdb_id': pdb_id,
            'output_dir': str(pdb_folder),
            'motif_types': {},
            'total_saved': 0,
            'total_failed': 0,
        }
        
        if not loaded_motifs:
            self.logger.warning("No motifs loaded to save")
            return stats
        
        self.logger.info(f"\nSaving all motif images for {pdb_id}")
        print("=" * 60)
        print(f"  SAVING MOTIF IMAGES - {pdb_id}")
        print("=" * 60)
        
        for motif_type, info in sorted(loaded_motifs.items()):
            motif_folder = self.create_motif_type_folder(pdb_folder, motif_type)
            
            motif_details = info.get('motif_details', [])
            count_saved = 0
            count_failed = 0
            
            print(f"\n  {motif_type}: Saving {len(motif_details)} instances...")
            
            for instance_no, detail in enumerate(motif_details, 1):
                success = self.save_instance_image(motif_folder, instance_no, 
                                                  motif_type, detail, structure_name,
                                                  representation=representation)
                if success:
                    count_saved += 1
                    stats['total_saved'] += 1
                else:
                    count_failed += 1
                    stats['total_failed'] += 1
            
            stats['motif_types'][motif_type] = {
                'saved': count_saved,
                'failed': count_failed,
                'total': len(motif_details),
                'folder': str(motif_folder),
            }
            
            print(f"    ✅ Saved {count_saved}/{len(motif_details)} instances")
            if count_failed > 0:
                print(f"    ⚠️  Failed: {count_failed}")
        
        print("\n" + "=" * 60)
        print(f"  SUMMARY")
        print("=" * 60)
        print(f"  Total saved: {stats['total_saved']}")
        print(f"  Total failed: {stats['total_failed']}")
        print(f"  Output folder: {pdb_folder}")
        print("=" * 60 + "\n")
        
        self.logger.success(f"Saved {stats['total_saved']} motif images to {pdb_folder}")
        
        return stats
    
    def save_motif_type_images(self, loaded_motifs: Dict, motif_type: str,
                              structure_name: str, pdb_id: str,
                              output_base_dir: Optional[str] = None,
                              representation: str = 'cartoon') -> Dict:
        """
        Save images for a specific motif type.
        
        Creates folder structure: plugin_dir/motif_images/pdb_id/MOTIF_TYPE/instance_*_info.png
        
        Args:
            loaded_motifs: Dictionary of loaded motifs from UnifiedMotifLoader
            motif_type: Motif type to save (e.g., 'HL', 'IL')
            structure_name: Name of structure in PyMOL
            pdb_id: PDB ID for folder naming
            output_base_dir: Base output directory (defaults to plugin_dir/motif_images)
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        
        Returns:
            Dictionary with save statistics
        """
        pdb_folder = self.create_folder_hierarchy(pdb_id, output_base_dir)
        
        stats = {
            'pdb_id': pdb_id,
            'motif_type': motif_type,
            'output_dir': str(pdb_folder),
            'saved': 0,
            'failed': 0,
        }
        
        motif_type_upper = motif_type.upper()
        
        if motif_type_upper not in loaded_motifs:
            self.logger.error(f"Motif type '{motif_type}' not found in loaded motifs")
            self.logger.info(f"Available: {', '.join(loaded_motifs.keys())}")
            return stats
        
        info = loaded_motifs[motif_type_upper]
        motif_details = info.get('motif_details', [])
        
        if not motif_details:
            self.logger.warning(f"No instances found for {motif_type}")
            return stats
        
        motif_folder = self.create_motif_type_folder(pdb_folder, motif_type_upper)
        
        self.logger.info(f"\nSaving {motif_type} images for {pdb_id}")
        print("=" * 60)
        print(f"  SAVING MOTIF IMAGES - {pdb_id}")
        print("=" * 60)
        print(f"  Motif Type: {motif_type_upper}")
        print(f"  Instances: {len(motif_details)}")
        print("-" * 60)
        
        for instance_no, detail in enumerate(motif_details, 1):
            success = self.save_instance_image(motif_folder, instance_no,
                                             motif_type_upper, detail, structure_name,
                                             representation=representation)
            if success:
                stats['saved'] += 1
            else:
                stats['failed'] += 1
        
        print("\n" + "=" * 60)
        print(f"  SUMMARY")
        print("=" * 60)
        print(f"  Total saved: {stats['saved']}/{len(motif_details)}")
        if stats['failed'] > 0:
            print(f"  Failed: {stats['failed']}")
        print(f"  Output folder: {motif_folder}")
        print("=" * 60 + "\n")
        
        self.logger.success(f"Saved {stats['saved']} {motif_type} images to {motif_folder}")
        
        return stats    
    def save_motif_instance(self, loaded_motifs: Dict, motif_type: str, 
                           instance_id: int, structure_name: str, pdb_id: str,
                           output_base_dir: Optional[str] = None,
                           representation: str = 'cartoon') -> bool:
        """
        Save image for a specific motif instance.
        
        Creates file: plugin_dir/motif_images/pdb_id/MOTIF_TYPE/instance_ID_info.png
        
        Args:
            loaded_motifs: Dictionary of loaded motifs from UnifiedMotifLoader
            motif_type: Motif type (e.g., 'HL', 'IL')
            instance_id: Instance number (1-based, as shown in rmv_summary)
            structure_name: Name of structure in PyMOL
            pdb_id: PDB ID for folder naming
            output_base_dir: Base output directory (defaults to plugin_dir/motif_images)
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        
        Returns:
            True if successful, False otherwise
        """
        try:
            pdb_folder = self.create_folder_hierarchy(pdb_id, output_base_dir)
            motif_type_upper = motif_type.upper()
            
            if motif_type_upper not in loaded_motifs:
                self.logger.error(f"Motif type '{motif_type}' not found in loaded motifs")
                return False
            
            info = loaded_motifs[motif_type_upper]
            motif_details = info.get('motif_details', [])
            
            if not motif_details:
                self.logger.error(f"No instances found for {motif_type}")
                return False
            
            # Validate instance ID
            if instance_id < 1 or instance_id > len(motif_details):
                self.logger.error(f"Instance ID {instance_id} out of range (1-{len(motif_details)})")
                return False
            
            motif_folder = self.create_motif_type_folder(pdb_folder, motif_type_upper)
            
            # Get the specific instance (instance_id is 1-based)
            detail = motif_details[instance_id - 1]
            
            self.logger.info(f"\nSaving {motif_type} instance #{instance_id}")
            print("=" * 60)
            print(f"  SAVING MOTIF INSTANCE - {pdb_id}")
            print("=" * 60)
            print(f"  Motif Type: {motif_type_upper}")
            print(f"  Instance: {instance_id}/{len(motif_details)}")
            print("-" * 60)
            
            success = self.save_instance_image(motif_folder, instance_id,
                                             motif_type_upper, detail, structure_name,
                                             representation=representation)
            
            if success:
                print(f"  ✅ Saved successfully")
                print(f"  Output folder: {motif_folder}")
                print("=" * 60 + "\n")
                self.logger.success(f"Saved {motif_type} instance #{instance_id}")
                return True
            else:
                print(f"  ❌ Failed to save instance")
                print("=" * 60 + "\n")
                self.logger.error(f"Failed to save {motif_type} instance #{instance_id}")
                return False
        
        except Exception as e:
            self.logger.error(f"Error saving motif instance: {e}")
            return False
    
    def save_current_view(self, filename: str, width: int = 2400, 
                          height: int = 1800, dpi: int = 300, ray: int = 0) -> bool:
        """
        Save the current PyMOL view to a high-resolution PNG file.
        
        Captures the exact current view (same rotation, angle, zoom) at
        high resolution without any modifications or automatic zooming.
        ray=0 preserves the exact OpenGL appearance on screen.
        
        Args:
            filename: Output filename (e.g., 'my_structure.png')
            width: Image width in pixels (default: 2400)
            height: Image height in pixels (default: 1800)
            dpi: Dots per inch for print quality (default: 300)
            ray: 0 = exact screen appearance, 1 = ray-traced (publication quality)
        
        Returns:
            True if successful, False otherwise
        """
        try:
            filepath = Path(filename)
            
            # Ensure output directory exists
            filepath.parent.mkdir(parents=True, exist_ok=True)
            
            # Save current view as high-resolution PNG
            # ray=0 preserves the exact on-screen appearance (no ray-tracing changes)
            self.cmd.png(str(filepath), 
                        width=width, 
                        height=height,
                        dpi=dpi,
                        ray=ray)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to save current view: {e}")
            return False