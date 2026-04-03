"""
User Annotation Provider

Loads motif data from user-uploaded annotation files (FR3D, RNAMotifScan, RNAMotifScanX, etc.).
Uses converters to transform external formats to standard MotifInstance format.
"""

import os
from copy import deepcopy
from pathlib import Path
from typing import Dict, List, Optional
from ..base_provider import BaseProvider, MotifInstance, DatabaseInfo, DatabaseSourceType, ResidueSpec
from .converters import FR3DConverter, RNAMotifScanConverter, RNAMotifScanXConverter, NoBIASConverter, MotifInstanceSimple


class UserAnnotationProvider(BaseProvider):
    """
    Provides motifs from user-uploaded annotation files.
    
    Supports:
    - FR3D output CSV files
    - RNAMotifScan output files (RMS)
    - RNAMotifScanX output files (RMSX)
    """
    
    def __init__(self, user_annotations_dir: str):
        """
        Initialize user annotation provider.
        
        Args:
            user_annotations_dir: Path to directory containing user annotation subdirectories
                                 (fr3d/, RNAMotifScan/, RNAMotifScanX/, etc.)
        """
        self.user_annotations_dir = Path(user_annotations_dir)
        self.user_annotations_dir.mkdir(parents=True, exist_ok=True)
        
        self._info = DatabaseInfo(
            id='user_annotations',
            name='User Annotations',
            description='User-uploaded motif annotation files from external tools',
            version='1.0.0',
            source_type=DatabaseSourceType.LOCAL_DIRECTORY,
        )
        
        # Supported tool formats
        self.supported_tools = ['fr3d', 'RNAMotifScan', 'RNAMotifScanX', 'NoBIAS']
        self.active_tool: Optional[str] = None  # Filter to specific tool if set
        
        # Filtering control for RMS, RMSX, and NoBIAS (added for on/off toggle)
        self.apply_rms_filtering: bool = True    # Default: filters ON for RMS
        self.apply_rmsx_filtering: bool = True   # Default: filters ON for RMSX
        self.apply_nobias_filtering: bool = True  # Default: filters ON for NoBIAS
        
        # Custom P-values for RMS, RMSX, and NoBIAS per motif type
        self.rms_custom_pvalues: Dict[str, float] = {}      # {motif_name: p_value}
        self.rmsx_custom_pvalues: Dict[str, float] = {}     # {motif_name: p_value}
        self.nobias_custom_pvalues: Dict[str, float] = {}   # {motif_name: p_value}
        
        # Override tool directory: when set, maps tool_name -> Path
        # Used for custom data paths (e.g., rmv_db 7 /path/to/rmsx/data)
        self.override_tool_dirs: Dict[str, Path] = {}
        
        # Track loaded data
        self._loaded_motifs_cache: Dict[str, Dict] = {}
        self._available_pdbs: List[str] = []
        self._motif_types: Dict[str, List[MotifInstance]] = {}
        self._initialized = False
    
    def set_active_tool(self, tool_name: Optional[str]):
        """Set active tool filter (fr3d, rms, rmsx, or None for all)."""
        if tool_name is not None:
            tool_name_lower = tool_name.lower()
            # Map display names to internal names (with aliases)
            tool_map = {
                'fr3d': 'fr3d',
                'rnamotifscan': 'RNAMotifScan',
                'rms': 'RNAMotifScan',  # Alias for RNAMotifScan
                'rnamotifscanx': 'RNAMotifScanX',
                'rmsx': 'RNAMotifScanX',  # Alias for RNAMotifScanX
                'nobias': 'NoBIAS',  # NoBIAS analysis
            }
            self.active_tool = tool_map.get(tool_name_lower)
        else:
            self.active_tool = None
    
    @property
    def info(self) -> DatabaseInfo:
        """Get database metadata."""
        return self._info
    
    def initialize(self) -> bool:
        """
        Initialize the provider by scanning for annotation files.
        
        Returns:
            True if initialization successful
        """
        try:
            # Scan for available PDB files
            pdbs = set()
            
            for tool_name in self.supported_tools:
                tool_dir = self.user_annotations_dir / tool_name
                if not tool_dir.exists():
                    continue
                
                if tool_name.lower() == 'fr3d':
                    # FR3D: Scan files directly in fr3d/ folder
                    for file_path in tool_dir.glob('*'):
                        if file_path.is_file() and file_path.suffix in ['.csv', '.tsv', '.txt']:
                            # Extract PDB ID from filename (usually first 4 chars)
                            filename = file_path.stem.lower()
                            pdb_id = filename.split('_')[0]
                            if len(pdb_id) >= 4:
                                pdbs.add(pdb_id.upper())
                
                elif tool_name == 'RNAMotifScan':
                    # RNAMotifScan: Scan motif-type folders for Res_* files
                    for motif_folder in tool_dir.iterdir():
                        if not motif_folder.is_dir():
                            continue
                        
                        for res_file in motif_folder.glob('Res_*'):
                            # Extract PDB ID from filename: Res_1s72 -> 1s72
                            filename = res_file.stem.lower()
                            if filename.startswith('res_'):
                                pdb_id = filename.replace('res_', '')
                                if len(pdb_id) >= 4:
                                    pdbs.add(pdb_id.upper())
                
                elif tool_name == 'RNAMotifScanX':
                    # RNAMotifScanX: Scan motif-type folders for result_*.log files
                    for motif_folder in tool_dir.iterdir():
                        if not motif_folder.is_dir():
                            continue
                        
                        for result_file in motif_folder.glob('result_*.log'):
                            # Try to extract PDB ID from file content
                            try:
                                with open(result_file, 'r') as f:
                                    first_line = f.readline()
                                    # Look for pattern like "1S72_0:..." in fragment IDs
                                    import re
                                    match = re.search(r'(\d\w{3})_', first_line.upper())
                                    if match:
                                        pdb_id = match.group(1)
                                        pdbs.add(pdb_id.upper())
                            except Exception:
                                continue
                
                elif tool_name == 'NoBIAS':
                    # NoBIAS: Scan flat files named <pdb>_<motif>_nobias.txt
                    for file_path in tool_dir.glob('*_nobias.txt'):
                        if file_path.is_file():
                            # Extract PDB ID from filename: 1s72_k-turn_nobias.txt -> 1s72
                            filename = file_path.stem.lower()
                            pdb_id = filename.split('_')[0]
                            if len(pdb_id) >= 4:
                                pdbs.add(pdb_id.upper())
            
            self._available_pdbs = sorted(list(pdbs))
            self._initialized = True
            return True
        except Exception as e:
            print(f"Warning: Could not initialize UserAnnotationProvider: {e}")
            return False
    
    def set_rms_custom_pvalues(self, pvalues: Dict[str, float]):
        """Set custom P-value thresholds for RMS motif types.
        
        Args:
            pvalues: Dict mapping motif names to P-value thresholds
                    e.g., {'C-LOOP': 0.05, 'KINK-TURN': 0.02}
        """
        self.rms_custom_pvalues = pvalues.copy() if pvalues else {}
    
    def set_rmsx_custom_pvalues(self, pvalues: Dict[str, float]):
        """Set custom P-value thresholds for RMSX motif types.
        
        Args:
            pvalues: Dict mapping motif names to P-value thresholds
        """
        self.rmsx_custom_pvalues = pvalues.copy() if pvalues else {}
    
    def set_nobias_custom_pvalues(self, pvalues: Dict[str, float]):
        """Set custom P-value thresholds for NoBIAS motif types.
        
        Args:
            pvalues: Dict mapping motif names to P-value thresholds
        """
        self.nobias_custom_pvalues = pvalues.copy() if pvalues else {}
    
    def get_available_motif_types(self) -> List[str]:
        """Get all available motif types across all loaded files."""
        all_types = set()
        for motif_list in self._motif_types.values():
            for instance in motif_list:
                all_types.add(instance.motif_id)
        return sorted(list(all_types))
    
    def get_motif_type(self, type_id: str) -> Optional[Dict]:
        """Get all instances of a specific motif type."""
        all_instances = []
        for instances in self._motif_types.values():
            for inst in instances:
                if inst.motif_id == type_id:
                    all_instances.append(inst)
        
        if all_instances:
            return {
                'type_id': type_id,
                'instances': all_instances,
                'count': len(all_instances)
            }
        return None
    
    def get_motifs_for_pdb(self, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """
        Get motifs for a PDB ID from user annotation files.
        
        Searches for files matching PDB_ID pattern in tool subdirectories.
        
        For RNAMotifScan: Looks in motif-type folders (Kturn/, c_loop/, etc.) for Res_<pdb_id> files
        For RNAMotifScanX: Looks in motif-type folders (k-turn_consensus/, etc.) for result_*.log files
        For FR3D: Looks for <pdb_id>*.csv files directly in fr3d/ folder
        
        Args:
            pdb_id: PDB ID to search for
            
        Returns:
            Dict mapping motif types to lists of MotifInstance objects
        """
        pdb_id_lower = pdb_id.lower()
        all_motifs = {}
        
        # If active_tool is set, only load from that tool
        tools_to_load = [self.active_tool] if self.active_tool else self.supported_tools
        
        # Search each tool subdirectory
        for tool_name in tools_to_load:
            # Use override directory if set for this tool, otherwise default
            if tool_name in self.override_tool_dirs:
                tool_dir = self.override_tool_dirs[tool_name]
            else:
                tool_dir = self.user_annotations_dir / tool_name
            if not tool_dir.exists():
                continue
            
            # Handle single-file override: load the file directly
            if tool_dir.is_file():
                try:
                    motif_type = tool_dir.stem  # e.g., "1S72_0_kturn"
                    motifs = self._load_file(tool_dir, tool_name, pdb_id, motif_type)
                    all_motifs.update(motifs)
                except Exception as e:
                    print(f"Warning: Could not load {tool_dir}: {e}")
                continue
            
            if tool_name.lower() == 'fr3d':
                # FR3D: Look for files directly in fr3d/ folder
                # Search case-insensitively for PDB files
                for file_path in tool_dir.iterdir():
                    if file_path.is_file() and file_path.suffix in ['.csv', '.tsv', '.txt']:
                        # Check if filename starts with PDB ID (case-insensitive)
                        if file_path.stem.lower().startswith(pdb_id_lower):
                            try:
                                motifs = self._load_file(file_path, tool_name, pdb_id)
                                all_motifs.update(motifs)
                            except Exception as e:
                                print(f"Warning: Could not load {file_path}: {e}")
                                continue
            
            elif tool_name == 'RNAMotifScan':
                # RNAMotifScan: Look in motif-type subfolders for Res_<pdb_id> files
                found_any = False
                for motif_folder in tool_dir.iterdir():
                    if not motif_folder.is_dir():
                        continue
                    
                    # Look for Res_<pdb_id> file
                    res_file = motif_folder / f"Res_{pdb_id_lower}"
                    if res_file.exists():
                        found_any = True
                        try:
                            motif_type = motif_folder.name  # e.g., "Kturn", "c_loop"
                            motifs = self._load_file(res_file, tool_name, pdb_id, motif_type)
                            all_motifs.update(motifs)
                        except Exception as e:
                            print(f"Warning: Could not load {res_file}: {e}")
                            continue
                
                # Fallback: scan flat files directly in the custom directory
                # Supports users who place RMS output files directly in their
                # custom path (e.g., Res_1s72.txt, 1S72_0_kturn.txt)
                if not found_any:
                    for flat_file in tool_dir.iterdir():
                        if not flat_file.is_file():
                            continue
                        if flat_file.suffix not in ['.txt', '.log', '.tsv', '']:
                            continue
                        if flat_file.stat().st_size == 0:
                            continue
                        try:
                            with open(flat_file, 'r') as f:
                                content = f.read(500)
                                if pdb_id.upper() not in content and pdb_id.lower() not in content:
                                    continue
                            motif_type = flat_file.stem
                            motifs = self._load_file(flat_file, tool_name, pdb_id, motif_type)
                            all_motifs.update(motifs)
                        except Exception as e:
                            print(f"Warning: Could not load {flat_file}: {e}")
                            continue
            
            elif tool_name == 'RNAMotifScanX':
                # RNAMotifScanX: Look in motif-type subfolders for result_*.log files containing pdb_id
                found_any = False
                for motif_folder in tool_dir.iterdir():
                    if not motif_folder.is_dir():
                        continue
                    
                    # Prioritize result_0_100_withbs.log (with binding site info)
                    # Then fall back to others: result_0_100.log, result_0_withbs.log, etc.
                    priority_files = [
                        'result_0_100_withbs.log',
                        'result_0_100.log',
                        'result_0_withbs.log',
                        'result_0.log'
                    ]
                    
                    result_file = None
                    for priority_filename in priority_files:
                        candidate = motif_folder / priority_filename
                        if candidate.exists() and candidate.stat().st_size > 0:
                            # Check if file contains the PDB ID
                            try:
                                with open(candidate, 'r') as f:
                                    content = f.read(500)
                                    if pdb_id.upper() in content or pdb_id.lower() in content:
                                        result_file = candidate
                                        break
                            except:
                                continue
                    
                    # If priority files not found, use first available
                    if not result_file:
                        for result_file_candidate in motif_folder.glob("result_*.log"):
                            if result_file_candidate.stat().st_size > 0:
                                try:
                                    with open(result_file_candidate, 'r') as f:
                                        content = f.read(500)
                                        if pdb_id.upper() in content or pdb_id.lower() in content:
                                            result_file = result_file_candidate
                                            break
                                except:
                                    continue
                    
                    # Load if file was found
                    if result_file:
                        found_any = True
                        try:
                            motif_type = motif_folder.name  # e.g., "k-turn_consensus"
                            motifs = self._load_file(result_file, tool_name, pdb_id, motif_type)
                            all_motifs.update(motifs)
                        except Exception as e:
                            print(f"Warning: Could not load {result_file}: {e}")
                            continue
                
                # Fallback: scan flat .txt/.log files directly in the custom directory
                # This supports users who place RMSX output files directly in their
                # custom path (e.g., Res_1s72.txt, 1S72_0_kturn.txt)
                if not found_any:
                    for flat_file in tool_dir.iterdir():
                        if not flat_file.is_file():
                            continue
                        if flat_file.suffix not in ['.txt', '.log', '.tsv']:
                            continue
                        if flat_file.stat().st_size == 0:
                            continue
                        try:
                            with open(flat_file, 'r') as f:
                                content = f.read(500)
                                if pdb_id.upper() not in content and pdb_id.lower() not in content:
                                    continue
                            # Infer motif type from filename
                            # e.g., "1S72_0_kturn.txt" -> "kturn", "Res_1s72.txt" -> filename stem
                            stem = flat_file.stem  # e.g., "Res_1s72" or "1S72_0_kturn"
                            motif_type = stem  # converter will clean this up
                            motifs = self._load_file(flat_file, tool_name, pdb_id, motif_type)
                            all_motifs.update(motifs)
                        except Exception as e:
                            print(f"Warning: Could not load {flat_file}: {e}")
                            continue
            
            elif tool_name == 'NoBIAS':
                # NoBIAS: Flat files named <pdb>_<motif>_nobias.txt in NoBIAS/ directory
                found_any = False
                for nobias_file in tool_dir.glob(f'*_nobias.txt'):
                    if not nobias_file.is_file() or nobias_file.stat().st_size == 0:
                        continue
                    # Check if file is for this PDB (filename starts with pdb_id)
                    if not nobias_file.stem.lower().startswith(pdb_id_lower):
                        continue
                    found_any = True
                    try:
                        motif_type = nobias_file.stem  # e.g., "1s72_k-turn_nobias"
                        motifs = self._load_file(nobias_file, tool_name, pdb_id, motif_type)
                        all_motifs.update(motifs)
                    except Exception as e:
                        print(f"Warning: Could not load {nobias_file}: {e}")
                        continue
                
                # Fallback: scan all files for PDB ID content (custom path support)
                if not found_any:
                    for flat_file in tool_dir.iterdir():
                        if not flat_file.is_file():
                            continue
                        if flat_file.suffix not in ['.txt', '.log', '.tsv']:
                            continue
                        if flat_file.stat().st_size == 0:
                            continue
                        try:
                            with open(flat_file, 'r') as f:
                                content = f.read(500)
                                if pdb_id.upper() not in content and pdb_id.lower() not in content:
                                    continue
                            motif_type = flat_file.stem
                            motifs = self._load_file(flat_file, tool_name, pdb_id, motif_type)
                            all_motifs.update(motifs)
                        except Exception as e:
                            print(f"Warning: Could not load {flat_file}: {e}")
                            continue
        
        # Convert MotifInstanceSimple to standard MotifInstance
        result = {}
        for motif_type, instances in all_motifs.items():
            result[motif_type] = [self._convert_instance(inst, pdb_id) for inst in instances]
        
        total_motifs = sum(len(motifs) for motifs in result.values())
        
        # Cache for later reference
        self._motif_types[pdb_id] = sum(result.values(), [])
        
        return result
    
    def get_available_pdb_ids(self) -> List[str]:
        """Get list of all PDB IDs with annotation files."""
        if not self._available_pdbs:
            self.initialize()
        return self._available_pdbs
    
    def get_motif_residues(self, pdb_id: str, motif_type: str, 
                          instance_id: str) -> List[ResidueSpec]:
        """
        Get residues for a specific motif instance.
        
        Args:
            pdb_id: PDB structure identifier
            motif_type: Motif type identifier
            instance_id: Instance identifier
            
        Returns:
            List of ResidueSpec objects
        """
        motifs = self.get_motifs_for_pdb(pdb_id)
        
        if motif_type not in motifs:
            return []
        
        for instance in motifs[motif_type]:
            if instance.instance_id == instance_id:
                return instance.residues
        
        return []
    
    def _load_file(self, file_path: Path, tool_name: str, pdb_id: str, motif_type: str = None) -> Dict[str, List[MotifInstanceSimple]]:
        """
        Load motifs from a specific file using appropriate converter.
        
        Args:
            file_path: Path to annotation file
            tool_name: Tool name ('fr3d', 'RNAMotifScan', 'RNAMotifScanX')
            pdb_id: PDB ID
            motif_type: Motif type (optional, used for RMS/RMSX/NoBIAS to pass folder name)
            
        Returns:
            Dict of motifs keyed by type
        """
        if tool_name.lower() == 'fr3d':
            return FR3DConverter.convert_file(str(file_path))
        elif tool_name == 'RNAMotifScan':
            return RNAMotifScanConverter.convert_file(str(file_path), motif_type, 
                                                      apply_filters=self.apply_rms_filtering,
                                                      custom_pvalues=self.rms_custom_pvalues)
        elif tool_name == 'RNAMotifScanX':
            return RNAMotifScanXConverter.convert_file(str(file_path), motif_type,
                                                       apply_filters=self.apply_rmsx_filtering,
                                                       custom_pvalues=self.rmsx_custom_pvalues)
        elif tool_name == 'NoBIAS':
            return NoBIASConverter.convert_file(str(file_path), motif_type,
                                                apply_filters=self.apply_nobias_filtering,
                                                custom_pvalues=self.nobias_custom_pvalues)
        else:
            raise ValueError(f"Unknown tool format: {tool_name}")
    
    def _convert_instance(self, simple_instance: MotifInstanceSimple, pdb_id: str) -> MotifInstance:
        """
        Convert MotifInstanceSimple to standard MotifInstance format.
        
        Args:
            simple_instance: Simple instance from converter
            pdb_id: PDB ID
            
        Returns:
            Standard MotifInstance object
        """
        # Convert residue tuples to ResidueSpec objects
        residues = []
        for nucleotide, res_num, chain in simple_instance.residues:
            residue = ResidueSpec(
                nucleotide=nucleotide or 'N',
                residue_number=res_num,
                chain=chain
            )
            residues.append(residue)
        
        # Create MotifInstance - use deepcopy for metadata to ensure complete independence
        instance = MotifInstance(
            motif_id=simple_instance.motif_id,
            instance_id=simple_instance.instance_id,
            residues=residues,
            pdb_id=pdb_id,
            annotation=simple_instance.annotation,
            metadata=deepcopy(simple_instance.metadata) if simple_instance.metadata else {}
        )
        
        return instance
    
    def is_available(self) -> bool:
        """Check if any user annotation files exist."""
        if not self.user_annotations_dir.exists():
            return False
        
        for tool_name in self.supported_tools:
            tool_dir = self.user_annotations_dir / tool_name
            if tool_dir.exists() and any(tool_dir.glob('*.csv')) or \
               any(tool_dir.glob('*.tsv')) or \
               any(tool_dir.glob('*.txt')):
                return True
        
        return False

