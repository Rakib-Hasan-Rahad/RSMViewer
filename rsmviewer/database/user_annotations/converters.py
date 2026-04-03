"""
User Annotation Converters

Converts external tool formats (FR3D, RNAMotifScan, etc.) to the standard MotifInstance format.

Each converter follows this pattern:
    1. Parse tool-specific format
    2. Extract motif name, residue positions, and metadata
    3. Convert to standard MotifInstance objects
    4. Return dict: {motif_type: [MotifInstance, ...]}
"""

import csv
from typing import Dict, List, Tuple
from pathlib import Path


# P-value thresholds for RMSX by motif family (from RNAMotifScanX paper Table 6)
RMSX_PVALUE_THRESHOLDS = {
    'KINK-TURN': 0.066,
    'K-TURN': 0.066,  # Alternative name
    'C-LOOP': 0.044,
    'SARCIN-RICIN': 0.040,
    'REVERSE KINK-TURN': 0.018,
    'REVERSE-KINK-TURN': 0.018,  # Alternative name
    'E-LOOP': 0.018,
}

# P-value thresholds for NoBIAS (same as RMSX)
NOBIAS_PVALUE_THRESHOLDS = dict(RMSX_PVALUE_THRESHOLDS)

# P-value thresholds for RNAMotifScan by motif family (from RNAMotifScan paper Table 3, page 9)
RMS_PVALUE_THRESHOLDS = {
    'KINK-TURN': 0.07,
    'K-TURN': 0.07,  # Alternative name
    'KTURN': 0.07,   # Folder name variant
    'C-LOOP': 0.04,
    'C_LOOP': 0.04,  # Folder name variant
    'SARCIN-RICIN': 0.02,
    'SARCIN_RICIN': 0.02,  # Folder name variant
    'SARCIN–RICIN': 0.02,  # With en-dash
    'REVERSE KINK-TURN': 0.14,
    'REVERSE-KINK-TURN': 0.14,  # With hyphen
    'REVERSE_KINK_TURN': 0.14,  # With underscore
    'REVERSE KTURN': 0.14,  # Alternative
    'REVERSE_KTURN': 0.14,  # Alternative
    'REVERSE-KTURN': 0.14,  # Alternative
    'E-LOOP': 0.13,
    'E_LOOP': 0.13,  # Folder name variant
}



class MotifInstanceSimple:
    """Lightweight MotifInstance for user annotations (before standardization)."""
    
    def __init__(self, motif_id: str, instance_id: str, residues: List[Tuple], 
                 annotation: str = "", metadata: Dict = None):
        self.motif_id = motif_id
        self.instance_id = instance_id
        self.residues = residues  # List of (nucleotide, residue_number, chain)
        self.annotation = annotation
        self.metadata = metadata or {}  # Store numeric fields: p_value, alignment_score, etc.

    
    def to_legacy_format(self) -> List[Dict]:
        """Convert to legacy format for PyMOL selector."""
        result = []
        current_chain = None
        residue_list = []
        
        for nucleotide, res_num, chain in self.residues:
            if chain != current_chain:
                if residue_list:
                    result.append({
                        'motif_id': self.motif_id,
                        'residues': residue_list,
                        'chain': current_chain,
                    })
                current_chain = chain
                residue_list = []
            residue_list.append(res_num)
        
        if residue_list:
            result.append({
                'motif_id': self.motif_id,
                'residues': residue_list,
                'chain': current_chain,
            })
        
        return result


class FR3DConverter:
    """Convert FR3D output CSV format to standard motif format.
    
    FR3D CSV format (comma-delimited):
    - Motif order: Sequential number
    - Motif type: Type of motif (e.g., "Hairpin", "Internal loop", "Bulge")
    - Resolution: Resolution value (e.g., "NA" or numeric)
    - Positions: Format "PDB_ID|model|chain|start-end"
    - Sequence: The nucleotide sequence
    - cWW/other: Count or descriptor
    - Description: Human-readable description
    
    Example:
    1,Hairpin,NA,"1S72|1|0|13-530","GCCAGCUGGUUGCG...",278,"Hairpin with 10 base pairs"
    """
    
    @staticmethod
    def parse_positions(positions_str: str) -> tuple:
        """
        Parse FR3D positions format: "PDB_ID|model|chain|start-end"
        Example: "1S72|1|0|13-530"
        
        Returns: (pdb_id, model, chain, start, end)
        """
        parts = positions_str.split('|')
        if len(parts) != 4:
            raise ValueError(f"Invalid FR3D positions format: {positions_str}")
        
        pdb_id, model, chain, range_str = parts
        
        # Parse range
        range_parts = range_str.split('-')
        if len(range_parts) != 2:
            raise ValueError(f"Invalid range format: {range_str}")
        
        try:
            start = int(range_parts[0])
            end = int(range_parts[1])
        except ValueError:
            raise ValueError(f"Invalid residue numbers: {range_str}")
        
        return pdb_id, model, chain, start, end
    
    @staticmethod
    def convert_file(csv_path: str) -> Dict[str, List[MotifInstanceSimple]]:
        """
        Convert FR3D CSV file to motif instances.
        
        Args:
            csv_path: Path to FR3D CSV file (comma-delimited)
            
        Returns:
            Dict mapping motif types to lists of MotifInstanceSimple
        """
        motifs_by_type = {}
        
        try:
            with open(csv_path, 'r', encoding='utf-8') as f:
                # FR3D files are comma-delimited
                reader = csv.DictReader(f)
                
                for row in reader:
                    try:
                        # Parse key fields
                        motif_order = row.get('Motif order', '').strip()
                        motif_type = row.get('Motif type', '').strip()
                        positions_str = row.get('Positions', '').strip().strip('"')
                        sequence = row.get('Sequence', '').strip().strip('"')
                        description = row.get('Description', '').strip().strip('"')
                        
                        if not motif_type or not positions_str:
                            continue
                        
                        # Parse positions
                        pdb_id, model, chain, start, end = FR3DConverter.parse_positions(positions_str)
                        
                        # Generate residues list from start to end
                        residues = []
                        for res_num in range(start, end + 1):
                            # Assign nucleotide from sequence if available
                            seq_idx = res_num - start
                            if 0 <= seq_idx < len(sequence):
                                nucleotide = sequence[seq_idx]
                            else:
                                nucleotide = 'N'  # Unknown
                            residues.append((nucleotide, res_num, chain))
                        
                        if not residues:
                            continue
                        
                        # Create instance ID
                        instance_id = f"FR3D_{pdb_id}_{chain}_{start}_{end}"
                        
                        # Build annotation with metadata
                        annotation = f"{description} | Range: {start}-{end}"
                        
                        # Create metadata dict
                        metadata = {
                            'positions': f"{start}-{end}",
                            'pdb_id': pdb_id,
                            'chain': chain,
                            'sequence_length': len(sequence),
                            'residue_count': len(residues),
                        }
                        
                        instance = MotifInstanceSimple(
                            motif_id=motif_type,
                            instance_id=instance_id,
                            residues=residues,
                            annotation=annotation,
                            metadata=metadata
                        )
                        
                        if motif_type not in motifs_by_type:
                            motifs_by_type[motif_type] = []
                        motifs_by_type[motif_type].append(instance)
                        
                    except Exception as e:
                        # Log but continue processing
                        continue
            
            return motifs_by_type
            
        except FileNotFoundError:
            raise FileNotFoundError(f"FR3D CSV file not found: {csv_path}")
        except Exception as e:
            raise Exception(f"Error parsing FR3D CSV file: {e}")


class RNAMotifScanConverter:
    """Convert RNAMotifScan (RMS) output format to standard motif format.
    
    RNAMotifScan output format (tab-separated):
    pdb_id\tlocation: 'chain'start-'chain'end/'chain'start-'chain'end\tScore: X\tP-value: Y\tFPR: Z
    
    Example:
    1s72_09	location: '0'71-'0'83/'0'91-'0'106	Score: 60.800	P-value: 0.00928	FPR: 0.00000
    
    The location field contains two regions separated by '/', each with format:
    'chain'start-'chain'end
    
    Features:
    - Stores p_value and alignment_score as numeric metadata
    - Filters by family-specific P-value thresholds (Table 3, RNAMotifScan paper)
    - Ranks by alignment_score (highest score first)
    """
    
    @staticmethod
    def parse_location(location_str: str) -> List[Tuple[str, int, int]]:
        """
        Parse RNAMotifScan location string.
        
        Format: 'chain'start-'chain'end/'chain'start-'chain'end
        Example: '0'71-'0'83/'0'91-'0'106
        
        Returns:
            List of (chain, start, end) tuples
        """
        regions = []
        
        # Split by '/' to get separate regions
        parts = location_str.split('/')
        
        for part in parts:
            part = part.strip()
            
            # Match pattern: 'chain'start-'chain'end
            # Extract chain and positions
            import re
            matches = re.findall(r"'([^']+)'(\d+)", part)
            
            if len(matches) >= 2:
                chain_start, start_pos = matches[0]
                chain_end, end_pos = matches[1]
                
                # Use the first chain (they should be the same)
                chain = chain_start
                start = int(start_pos)
                end = int(end_pos)
                
                regions.append((chain, start, end))
        
        return regions
    

    
    @staticmethod
    def convert_file(file_path: str, motif_type: str = None, apply_filters: bool = True, custom_pvalues: dict = None) -> Dict[str, List[MotifInstanceSimple]]:
        """
        Convert RNAMotifScan output file to motif instances.
        
        Args:
            file_path: Path to RNAMotifScan output file (e.g., Res_1s72)
            motif_type: Motif type name (inferred from folder name if not provided)
            apply_filters: Whether to apply P-value filtering (default: True)
            custom_pvalues: Optional dict of custom P-value thresholds {motif_name: p_value}
            
        Returns:
            Dict mapping motif types to lists of MotifInstanceSimple (filtered & deduplicated)
        """
        motifs_by_type = {}
        raw_instances = []
        custom_pvalues = custom_pvalues or {}
        
        # Infer motif type from file path if not provided
        if not motif_type:
            # Get parent folder name (e.g., "Kturn", "c_loop")
            motif_type_raw = Path(file_path).parent.name
            # Normalize: "Kturn" -> "KINK-TURN", "c_loop" -> "C-LOOP", etc.
            if motif_type_raw.lower() in ['kturn', 'kink_turn', 'kink-turn']:
                motif_type = 'KINK-TURN'
            elif motif_type_raw.lower() in ['c_loop', 'c-loop']:
                motif_type = 'C-LOOP'
            elif 'sarcin' in motif_type_raw.lower():
                motif_type = 'SARCIN-RICIN'
            elif 'reverse' in motif_type_raw.lower():
                motif_type = 'REVERSE-KINK-TURN'
            elif motif_type_raw.lower() in ['e_loop', 'e-loop']:
                motif_type = 'E-LOOP'
            else:
                motif_type = motif_type_raw.upper().replace('_', '-')
        else:
            motif_type = motif_type.upper().replace('_', '-')
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line_idx, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    try:
                        # Parse tab-separated fields
                        parts = line.split('\t')
                        if len(parts) < 4:
                            continue
                        
                        pdb_id = parts[0].strip()
                        location_field = parts[1].strip()
                        score_str = parts[2].strip() if len(parts) > 2 else "0.0"
                        pvalue_str = parts[3].strip() if len(parts) > 3 else "1.0"
                        
                        # Extract Score value
                        try:
                            score_val = float(score_str.replace("Score:", "").strip())
                        except ValueError:
                            score_val = 0.0
                        
                        # Extract P-value
                        try:
                            pvalue_val = float(pvalue_str.replace("P-value:", "").strip())
                        except ValueError:
                            pvalue_val = 1.0
                        
                        # Extract location (remove "location: " prefix)
                        if location_field.startswith("location:"):
                            location_str = location_field.replace("location:", "").strip()
                        else:
                            location_str = location_field
                        
                        # Parse location to get regions
                        regions = RNAMotifScanConverter.parse_location(location_str)
                        
                        if not regions:
                            continue
                        
                        # Collect all residues from all regions
                        residues = []
                        for chain, start, end in regions:
                            for res_num in range(start, end + 1):
                                residues.append(('N', res_num, chain))
                        
                        if not residues:
                            continue
                        
                        # Create instance ID
                        instance_id = f"RMS_{pdb_id}_{line_idx}"
                        
                        # Build annotation with metrics
                        annotation = f"Score: {score_val}, P-value: {pvalue_val}"
                        
                        # Create metadata dict with numeric fields
                        metadata = {
                            'p_value': pvalue_val,
                            'alignment_score': score_val,
                            'regions': regions,
                            'pdb_id': pdb_id,
                        }
                        
                        instance = MotifInstanceSimple(
                            motif_id=motif_type,
                            instance_id=instance_id,
                            residues=residues,
                            annotation=annotation,
                            metadata=metadata
                        )
                        
                        raw_instances.append(instance)
                        
                    except Exception:
                        continue
            
            if not raw_instances:
                return {}
            
            # Phase 1: Apply P-value filtering (if enabled)
            if apply_filters:
                # Use custom P-value if provided, otherwise use default threshold
                if motif_type in custom_pvalues:
                    threshold = custom_pvalues[motif_type]
                else:
                    threshold = RMS_PVALUE_THRESHOLDS.get(motif_type, 0.10)
                
                filtered_instances = [
                    inst for inst in raw_instances 
                    if inst.metadata.get('p_value', 1.0) <= threshold
                ]
                
                # Sort by alignment score (highest first)
                filtered_instances.sort(
                    key=lambda x: x.metadata.get('alignment_score', 0.0),
                    reverse=True
                )
                
                print(f"[RMS] {motif_type}: {len(raw_instances)} total → {len(filtered_instances)} after P-value filter (threshold={threshold})")
            else:
                # No filtering: use all raw instances, sorted by alignment score
                print(f"[RMS] {motif_type}: {len(raw_instances)} total (filtering disabled - showing raw data)")
                filtered_instances = raw_instances
                filtered_instances.sort(
                    key=lambda x: x.metadata.get('alignment_score', 0.0),
                    reverse=True
                )
            
            if filtered_instances:
                motifs_by_type[motif_type] = filtered_instances
            
            return motifs_by_type
            
        except FileNotFoundError:
            raise FileNotFoundError(f"RNAMotifScan output file not found: {file_path}")
        except Exception as e:
            raise Exception(f"Error parsing RNAMotifScan output file: {e}")


class RNAMotifScanXConverter:
    """Convert RNAMotifScanX (RMSX) output format to standard motif format.
    
    RNAMotifScanX output format (tab-separated with header):
    #fragment_ID	aligned_regions	alignment_score	P-value
    
    Example:
    1S72_0:75-85_89-98_58-60	0:'0'77-4:'0'81,13:'0'93-20:'0'100	144.8	0.00733485
    
    fragment_ID format: PDB_chain:residue_ranges (underscore-separated)
    aligned_regions format: index:'chain'start-index:'chain'end (comma-separated pairs)
    
    Features:
    - Stores p_value and alignment_score as numeric metadata
    - Parses aligned_regions; falls back to fragment_id if empty
    - Filters by family-specific P-value thresholds
    - Ranks by alignment_score (highest score first)
    """
    
    @staticmethod
    def parse_fragment_id(fragment_id: str) -> Tuple[str, str, List[Tuple[int, int]]]:
        """
        Parse RNAMotifScanX fragment ID.
        
        Format: PDB_chain:range1_range2_range3
        Example: 1S72_0:75-85_89-98_58-60
        
        Returns:
            (pdb_id, chain, [(start, end), ...])
        """
        import re
        
        # Split by ':'
        parts = fragment_id.split(':')
        if len(parts) != 2:
            return None, None, []
        
        pdb_chain = parts[0]
        ranges_str = parts[1]
        
        # Extract PDB ID and chain
        # Format: 1S72_0:RANGE_RANGE... or PDB_CHAIN
        pdb_parts = pdb_chain.split('_')
        if len(pdb_parts) >= 2:
            pdb_id = pdb_parts[0]
            chain = pdb_parts[1]
        else:
            pdb_id = pdb_parts[0]
            chain = '0'
        
        # Parse ranges - ranges_str has format: RANGE_RANGE_RANGE
        ranges = []
        for range_str in ranges_str.split('_'):
            match = re.match(r'(\d+)-(\d+)', range_str)
            if match:
                start = int(match.group(1))
                end = int(match.group(2))
                ranges.append((start, end))
        
        return pdb_id, chain, ranges
    
    @staticmethod
    def parse_aligned_regions(aligned_regions: str) -> List[Tuple[int, int]]:
        """
        Parse RMSX aligned_regions format.
        
        Format: motif_idx:'chain'res-motif_idx:'chain'res,comma-separated
        Example: 2:'0'1436-5:'0'1439,6:'0'1687-13:'0'1694,14:'0'1425-19:'0'1430
        
        Returns:
            List of (start_res, end_res) tuples from structure coordinates
        """
        import re
        
        if not aligned_regions or aligned_regions.strip() == '':
            return []
        
        ranges = []
        try:
            # Split by comma to get region pairs
            region_pairs = aligned_regions.split(',')
            
            for pair in region_pairs:
                # Match pattern: NUMBER:'DIGIT'NUMBER-NUMBER:'DIGIT'NUMBER
                # Example: 2:'0'1436-5:'0'1439
                match = re.search(r"\d+:'[^']*'(\d+)-\d+:'[^']*'(\d+)", pair)
                if match:
                    start_res = int(match.group(1))
                    end_res = int(match.group(2))
                    ranges.append((start_res, end_res))
        except Exception:
            return []
        
        return ranges
    
    @staticmethod
    def convert_file(file_path: str, motif_type: str = None, apply_filters: bool = True, custom_pvalues: dict = None) -> Dict[str, List[MotifInstanceSimple]]:
        """
        Convert RNAMotifScanX output file to motif instances with filtering.
        
        Args:
            file_path: Path to RNAMotifScanX output file (e.g., result_0_100.log)
            motif_type: Motif type name (inferred from folder name if not provided)
            apply_filters: Whether to apply P-value filtering (default: True)
            custom_pvalues: Optional dict of custom P-value thresholds {motif_name: p_value}
            
        Returns:
            Dict mapping motif types to lists of filtered MotifInstanceSimple
        """
        motifs_by_type = {}
        raw_instances = []  # Store all instances before filtering
        custom_pvalues = custom_pvalues or {}
        
        # Clean motif type name
        if not motif_type:
            # Get parent folder name (e.g., "k-turn_consensus")
            folder_name = Path(file_path).parent.name
            motif_type = folder_name
        
        # Always clean: remove _consensus suffix and convert to uppercase
        # "k-turn_consensus" → "K-TURN"
        # "c-loop_consensus" → "C-LOOP"
        # "sarcin-ricin_consensus" → "SARCIN-RICIN"
        motif_type_clean = motif_type.replace('_consensus', '').upper()
        
        # For flat files the motif_type may be a filename stem like
        # "1S72_0_kturn" or "Res_1s72".  Try to extract a known motif
        # keyword from the parts.
        _KNOWN_MOTIFS = {
            'KTURN': 'K-TURN', 'K-TURN': 'K-TURN', 'KINK-TURN': 'K-TURN',
            'CLOOP': 'C-LOOP', 'C-LOOP': 'C-LOOP', 'C_LOOP': 'C-LOOP',
            'SARCIN': 'SARCIN-RICIN', 'SARCIN-RICIN': 'SARCIN-RICIN',
            'REVERSE-KTURN': 'REVERSE-K-TURN', 'REVERSE_KTURN': 'REVERSE-K-TURN',
            'ELOOP': 'E-LOOP', 'E-LOOP': 'E-LOOP', 'E_LOOP': 'E-LOOP',
        }
        if motif_type_clean not in _KNOWN_MOTIFS:
            # Scan individual tokens for a known keyword
            for token in motif_type_clean.replace('-', '_').split('_'):
                if token in _KNOWN_MOTIFS:
                    motif_type_clean = _KNOWN_MOTIFS[token]
                    break
        else:
            motif_type_clean = _KNOWN_MOTIFS[motif_type_clean]
        
        motif_type = motif_type_clean
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                # Skip header lines and empty lines
                for line in f:
                    line = line.strip()
                    
                    # Skip header, comments, and empty lines
                    if not line or line.startswith('#') or line.startswith('No base-stacking'):
                        continue
                    
                    try:
                        # Parse tab-separated fields
                        parts = line.split('\t')
                        if len(parts) < 3:
                            continue
                        
                        fragment_id = parts[0].strip()
                        aligned_regions = parts[1].strip()
                        score_str = parts[2].strip()
                        pvalue_str = parts[3].strip() if len(parts) > 3 else "1.0"
                        
                        # Parse numeric values
                        try:
                            alignment_score = float(score_str)
                        except ValueError:
                            alignment_score = 0.0
                        
                        try:
                            p_value = float(pvalue_str)
                        except ValueError:
                            p_value = 1.0  # Default to worst value if unparseable
                        
                        # Parse fragment ID to get PDB ID, chain, and ranges
                        pdb_id, chain, ranges = RNAMotifScanXConverter.parse_fragment_id(fragment_id)
                        
                        if not ranges:
                            continue
                        
                        # Try to parse aligned_regions first; fallback to fragment_id
                        aligned_ranges = RNAMotifScanXConverter.parse_aligned_regions(aligned_regions)
                        ranges_to_use = aligned_ranges if aligned_ranges else ranges
                        
                        # Collect all residues from ranges
                        residues = []
                        for start, end in ranges_to_use:
                            for res_num in range(start, end + 1):
                                residues.append(('N', res_num, chain))
                        
                        # Create instance ID from fragment_id
                        instance_id = f"RMSX_{fragment_id.replace(':', '_').replace('-', '_')}"
                        
                        # Build annotation
                        annotation = f"Score: {alignment_score}, P-value: {p_value}"
                        
                        # Create metadata dict with numeric fields
                        metadata = {
                            'p_value': p_value,
                            'alignment_score': alignment_score,
                            'aligned_regions': aligned_ranges,  # Store parsed tuples, not raw string
                            'fragment_id': fragment_id,
                            'pdb_id': pdb_id,
                            'chain': chain,
                        }
                        
                        instance = MotifInstanceSimple(
                            motif_id=motif_type,
                            instance_id=instance_id,
                            residues=residues,
                            annotation=annotation,
                            metadata=metadata
                        )
                        
                        raw_instances.append(instance)
                        
                    except Exception as e:
                        continue
            
            # Phase 4: Apply P-value filtering (if enabled)
            if apply_filters:
                # Use custom P-value if provided, otherwise use default threshold
                if motif_type in custom_pvalues:
                    threshold = custom_pvalues[motif_type]
                else:
                    threshold = RMSX_PVALUE_THRESHOLDS.get(motif_type, 0.05)
                
                filtered_instances = [
                    inst for inst in raw_instances 
                    if inst.metadata.get('p_value', 1.0) <= threshold
                ]
                
                # Sort by alignment score (highest first)
                filtered_instances.sort(
                    key=lambda x: x.metadata.get('alignment_score', 0.0),
                    reverse=True
                )
                
                print(f"[RMSX] {motif_type}: {len(raw_instances)} total → {len(filtered_instances)} after P-value filter (threshold={threshold})")
            else:
                # No filtering: use all raw instances, sorted by alignment score
                print(f"[RMSX] {motif_type}: {len(raw_instances)} total (filtering disabled - showing raw data)")
                filtered_instances = raw_instances
                filtered_instances.sort(
                    key=lambda x: x.metadata.get('alignment_score', 0.0),
                    reverse=True
                )
            
            # Build result
            motifs_by_type[motif_type] = filtered_instances
            
            # Motifs processed: raw → filtered (by P-value) → sorted by alignment_score
            
            return motifs_by_type
            
        except FileNotFoundError:
            raise FileNotFoundError(f"RNAMotifScanX output file not found: {file_path}")
        except Exception as e:
            raise Exception(f"Error parsing RNAMotifScanX output file: {e}")
    


class NoBIASConverter:
    """Convert NoBIAS output format to standard motif format.
    
    NoBIAS uses the same tab-separated format as RNAMotifScanX:
    #fragment_ID	aligned_regions	alignment_score	P-value
    
    Files are flat in the NoBIAS/ directory, named <pdb>_<motif>_nobias.txt.
    Delegates parsing to RNAMotifScanXConverter, but uses NOBIAS_PVALUE_THRESHOLDS
    and prints [NoBIAS] prefixed messages.
    """
    
    @staticmethod
    def convert_file(file_path: str, motif_type: str = None, apply_filters: bool = True, custom_pvalues: dict = None) -> Dict[str, List[MotifInstanceSimple]]:
        """
        Convert NoBIAS output file to motif instances with filtering.
        
        Args:
            file_path: Path to NoBIAS output file (e.g., 1s72_k-turn_nobias.txt)
            motif_type: Motif type name (inferred from filename if not provided)
            apply_filters: Whether to apply P-value filtering (default: True)
            custom_pvalues: Optional dict of custom P-value thresholds {motif_name: p_value}
            
        Returns:
            Dict mapping motif types to lists of filtered MotifInstanceSimple
        """
        custom_pvalues = custom_pvalues or {}
        
        # Infer motif type from filename if not provided, or clean it if it
        # looks like a full NoBIAS filename stem (e.g., "1s72_k-turn_nobias")
        if not motif_type:
            stem = Path(file_path).stem
        else:
            stem = motif_type
        
        # Clean: remove _nobias suffix and PDB prefix
        cleaned = stem.lower()
        if cleaned.endswith('_nobias'):
            cleaned = cleaned[:-len('_nobias')]
        # Remove leading PDB ID (4 chars + underscore)
        if len(cleaned) > 5 and cleaned[4] == '_':
            cleaned = cleaned[5:]
        motif_type = cleaned
        
        # Delegate to RNAMotifScanXConverter for the actual parsing
        # but intercept the result to replace prefixes and use NoBIAS thresholds
        result = RNAMotifScanXConverter.convert_file(
            file_path, motif_type,
            apply_filters=False,  # We'll filter ourselves with NoBIAS thresholds
            custom_pvalues=custom_pvalues
        )
        
        if not apply_filters:
            # Re-label instances with NoBIAS prefix
            for motif_name, instances in result.items():
                for inst in instances:
                    inst.instance_id = inst.instance_id.replace('RMSX_', 'NoBIAS_', 1)
                print(f"[NoBIAS] {motif_name}: {len(instances)} total (filtering disabled - showing raw data)")
            return result
        
        # Apply NoBIAS-specific P-value filtering
        filtered_result = {}
        for motif_name, instances in result.items():
            if motif_name in custom_pvalues:
                threshold = custom_pvalues[motif_name]
            else:
                threshold = NOBIAS_PVALUE_THRESHOLDS.get(motif_name, 0.05)
            
            filtered = [
                inst for inst in instances
                if inst.metadata.get('p_value', 1.0) <= threshold
            ]
            filtered.sort(key=lambda x: x.metadata.get('alignment_score', 0.0), reverse=True)
            
            # Re-label instances with NoBIAS prefix
            for inst in filtered:
                inst.instance_id = inst.instance_id.replace('RMSX_', 'NoBIAS_', 1)
            
            print(f"[NoBIAS] {motif_name}: {len(instances)} total → {len(filtered)} after P-value filter (threshold={threshold})")
            
            if filtered:
                filtered_result[motif_name] = filtered
        
        return filtered_result
