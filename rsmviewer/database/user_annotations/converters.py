"""
User Annotation Converters

Converts FR3D and RNAMotifScanX formats to the standard MotifInstance format.

Each converter follows this pattern:
    1. Parse tool-specific format
    2. Extract motif name, residue positions, and metadata
    3. Convert to standard MotifInstance objects
    4. Return dict: {motif_type: [MotifInstance, ...]}
"""

import csv
import re
from typing import Dict, List, Tuple
from pathlib import Path


# P-value thresholds for RMSX by motif family (from RNAMotifScanX paper Table 6).
# Keys are the canonical family names produced by RNAMotifScanXConverter.
RMSX_PVALUE_THRESHOLDS = {
    'K-TURN': 0.066,
    'C-LOOP': 0.044,
    'SARCIN-RICIN': 0.040,
    'REVERSE-K-TURN': 0.018,
    'E-LOOP': 0.018,
}
RMSX_DEFAULT_PVALUE = 0.05

# Spelling-insensitive aliases (letters only, upper-case) -> canonical family.
_RMSX_FAMILY_ALIASES = {
    'KTURN': 'K-TURN', 'KINKTURN': 'K-TURN',
    'CLOOP': 'C-LOOP',
    'SARCIN': 'SARCIN-RICIN', 'SARCINRICIN': 'SARCIN-RICIN',
    'REVERSEKTURN': 'REVERSE-K-TURN', 'REVERSEKINKTURN': 'REVERSE-K-TURN',
    'ELOOP': 'E-LOOP',
}


def canonical_rmsx_family(name: str) -> str:
    """Map any spelling of an RMSX family (e.g. 'reverse-kturn_consensus',
    'Reverse Kink-Turn', 'K-TURN') to its canonical name, or the upper-cased
    input when it is not a known family."""
    text = str(name or '').strip().upper()
    if text.endswith('_CONSENSUS'):
        text = text[:-len('_CONSENSUS')]
    letters = re.sub(r'[^A-Z]', '', text)
    return _RMSX_FAMILY_ALIASES.get(letters, text)


def rmsx_pvalue_threshold(motif_type: str, custom_pvalues: dict = None) -> float:
    """P-value cutoff for a family: the config/custom value if given, else the
    paper default, else RMSX_DEFAULT_PVALUE. Family names match regardless of
    spelling, so 'REVERSE-KINK-TURN' in a config applies to REVERSE-K-TURN."""
    family = canonical_rmsx_family(motif_type)
    for name, value in (custom_pvalues or {}).items():
        if canonical_rmsx_family(name) == family:
            return float(value)
    return RMSX_PVALUE_THRESHOLDS.get(family, RMSX_DEFAULT_PVALUE)


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
    """Convert FR3D output formats to standard motif format.

    Supported formats:
    1. BGSU loops CSV (downloaded from rna.bgsu.edu/rna3dhub/loops/download/<PDB>):
       "HL_1S72_001","1S72|1|0|U|55,1S72|1|0|G|56,..."
       Loop-ID prefix: HL = Hairpin Loop, IL = Internal Loop, J3/J4/... = Junction

    2. FR3D motif CSV (range-based, comma-delimited):
       1,Hairpin,NA,"1S72|1|0|13-530","GCCAGCUGGUUGCG...",278,"Hairpin with 10 base pairs"

    3. FR3D pairwise TXT (basepair annotations from NA_pairwise_interactions):
       1S72|1|A|G|71\tcWW\t1S72|1|A|C|83\t0
    """

    # Map BGSU loop-ID prefix to human-readable motif type name
    _LOOP_TYPE_NAMES = {
        'HL': 'Hairpin Loop',
        'IL': 'Internal Loop',
        'J3': '3-Way Junction',
        'J4': '4-Way Junction',
        'J5': '5-Way Junction',
        'J6': '6-Way Junction',
        'J7': '7-Way Junction',
        'J8': '8-Way Junction',
    }

    @staticmethod
    def _detect_bgsu_loops(file_path: str) -> bool:
        """Return True if the file looks like a BGSU loops download CSV."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
            # BGSU loops CSV starts with a quoted loop ID like "HL_" or "IL_" or "J3_"
            return bool(re.match(r'^"?(HL|IL|J\d)_', first_line))
        except Exception:
            return False

    @staticmethod
    def _detect_fr3d_search_csv(file_path: str) -> bool:
        """Return True if the file looks like FR3D search candidate CSV output."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                header = f.readline().strip()
            header_lower = header.lower()
            return 'similarity order' in header_lower and 'position 1' in header_lower
        except Exception:
            return False

    @staticmethod
    def _motif_type_from_search_filename(file_path: str) -> str:
        """Infer motif type from a FR3D query-search CSV filename."""
        stem = Path(file_path).stem
        lowered = stem.lower()
        prefix = '_fr3d_query_'
        if prefix in lowered:
            idx = lowered.find(prefix)
            suffix = stem[idx + len(prefix):]
            if suffix:
                return suffix.replace('_', ' ').strip().upper()
        return stem.replace('_', ' ').strip().upper()

    @staticmethod
    def _convert_fr3d_search_csv(file_path: str) -> Dict[str, List[MotifInstanceSimple]]:
        """Convert FR3D search candidate CSV output to motif instances."""
        motifs_by_type: Dict[str, List[MotifInstanceSimple]] = {}
        motif_type = FR3DConverter._motif_type_from_search_filename(file_path)
        instances: List[MotifInstanceSimple] = []

        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            headers = reader.fieldnames or []
            position_columns = [h for h in headers if h and h.lower().startswith('position ')]

            for row in reader:
                residues = []
                unit_ids = []

                for col in position_columns:
                    unit_id = str(row.get(col, '') or '').strip()
                    if not unit_id:
                        continue
                    unit_ids.append(unit_id)
                    try:
                        _pdb_id, chain, residue_number, nucleotide = FR3DConverter._parse_unit_id(unit_id)
                    except ValueError:
                        continue
                    residues.append((nucleotide, residue_number, chain))

                if not residues:
                    continue

                similarity_order = str(row.get('Similarity order', '') or '').strip()
                discrepancy_text = str(row.get('Discrepancy', '') or '').strip()
                discrepancy_value = None
                if discrepancy_text:
                    try:
                        discrepancy_value = float(discrepancy_text)
                    except ValueError:
                        discrepancy_value = None

                instance_suffix = similarity_order if similarity_order else str(len(instances) + 1)
                instance_id = f"FR3D_SEARCH_{motif_type.replace(' ', '_')}_{instance_suffix}"

                metadata = {
                    'source_format': 'fr3d_search_csv',
                    'similarity_order': similarity_order,
                    'discrepancy': discrepancy_value,
                    'unit_ids': unit_ids,
                }

                annotation = f"FR3D search candidate {instance_suffix}"
                if discrepancy_value is not None:
                    annotation += f" | discrepancy={discrepancy_value:.4f}"

                instances.append(
                    MotifInstanceSimple(
                        motif_id=motif_type,
                        instance_id=instance_id,
                        residues=residues,
                        annotation=annotation,
                        metadata=metadata,
                    )
                )

        if instances:
            motifs_by_type[motif_type] = instances
        return motifs_by_type

    @staticmethod
    def _convert_bgsu_loops_csv(file_path: str) -> Dict[str, List[MotifInstanceSimple]]:
        """Parse the BGSU loops download CSV into structural motif instances.

        Format per line:
            "HL_1S72_001","1S72|1|0|U|55,1S72|1|0|G|56,..."

        Each residue token: pdb_id|model|chain|nucleotide|residue_number
        """
        motifs_by_type: Dict[str, List[MotifInstanceSimple]] = {}

        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                stripped = line.strip()
                if not stripped:
                    continue

                # Strip surrounding quotes and split at first ","
                # Line format: "HL_1S72_001","res1,res2,..."
                parts = re.split(r'","', stripped.strip('"'))
                if len(parts) < 2:
                    continue

                loop_id = parts[0].strip('"')
                residues_str = parts[1].strip('"')

                # Determine motif type from loop_id prefix (HL_, IL_, J3_, …)
                prefix_match = re.match(r'^(HL|IL|J\d+)_', loop_id)
                if not prefix_match:
                    continue
                prefix = prefix_match.group(1)
                motif_type = FR3DConverter._LOOP_TYPE_NAMES.get(prefix, prefix)

                # Parse each residue token
                residues = []
                for token in residues_str.split(','):
                    token = token.strip()
                    if not token:
                        continue
                    token_parts = token.split('|')
                    if len(token_parts) < 5:
                        continue
                    nucleotide = token_parts[3].strip() or 'N'
                    try:
                        res_num = int(re.match(r'^(-?\d+)', token_parts[4].strip()).group(1))
                    except (AttributeError, ValueError):
                        continue
                    chain = token_parts[2].strip()
                    residues.append((nucleotide, res_num, chain))

                if not residues:
                    continue

                instance = MotifInstanceSimple(
                    motif_id=motif_type,
                    instance_id=f"FR3D_{loop_id}",
                    residues=residues,
                    annotation=f"FR3D/BGSU {motif_type} ({loop_id})",
                    metadata={'loop_id': loop_id, 'source_format': 'bgsu_loops_csv'},
                )
                motifs_by_type.setdefault(motif_type, []).append(instance)

        return motifs_by_type

    @staticmethod
    def _parse_unit_id(unit_id: str) -> Tuple[str, str, int, str]:
        """Parse FR3D unit ID into (pdb_id, chain, residue_number, nucleotide)."""
        parts = unit_id.split('|')
        if len(parts) < 5:
            raise ValueError(f"Invalid FR3D unit ID: {unit_id}")

        pdb_id = parts[0].strip().upper()
        chain = parts[2].strip()
        nucleotide = parts[3].strip() or 'N'
        residue_token = parts[4].strip()

        match = re.match(r'^(-?\d+)', residue_token)
        if not match:
            raise ValueError(f"Invalid FR3D residue token: {residue_token}")
        residue_number = int(match.group(1))

        return pdb_id, chain, residue_number, nucleotide

    @staticmethod
    def _convert_pairwise_txt(file_path: str) -> Dict[str, List[MotifInstanceSimple]]:
        """Convert FR3D pairwise interaction TXT output to motif instances."""
        motifs_by_type: Dict[str, List[MotifInstanceSimple]] = {}

        with open(file_path, 'r', encoding='utf-8') as f:
            for idx, line in enumerate(f, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith('#'):
                    continue

                cols = stripped.split('\t')
                if len(cols) < 3:
                    continue

                unit1 = cols[0].strip()
                interaction = cols[1].strip()
                unit2 = cols[2].strip()
                if not unit1 or not interaction or not unit2:
                    continue

                try:
                    pdb1, chain1, res1, nt1 = FR3DConverter._parse_unit_id(unit1)
                    pdb2, chain2, res2, nt2 = FR3DConverter._parse_unit_id(unit2)
                except ValueError:
                    continue

                # Keep a predictable motif key for summaries and filtering.
                motif_type = interaction.upper()
                instance_id = f"FR3D_{pdb1}_{chain1}{res1}_{chain2}{res2}_{idx}"
                residues = [(nt1, res1, chain1), (nt2, res2, chain2)]
                annotation = f"FR3D pairwise interaction {interaction}"

                metadata = {
                    'source_format': 'fr3d_pairwise_txt',
                    'unit_1': unit1,
                    'unit_2': unit2,
                    'interaction': interaction,
                    'pdb_id_1': pdb1,
                    'pdb_id_2': pdb2,
                }

                instance = MotifInstanceSimple(
                    motif_id=motif_type,
                    instance_id=instance_id,
                    residues=residues,
                    annotation=annotation,
                    metadata=metadata,
                )

                motifs_by_type.setdefault(motif_type, []).append(instance)

        return motifs_by_type
    
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

        file_ext = Path(csv_path).suffix.lower()
        if file_ext == '.txt':
            return FR3DConverter._convert_pairwise_txt(csv_path)

        # Detect BGSU loops download format before falling back to the motif CSV parser
        if FR3DConverter._detect_bgsu_loops(csv_path):
            return FR3DConverter._convert_bgsu_loops_csv(csv_path)

        # Detect FR3D query-search CSV candidate output.
        if FR3DConverter._detect_fr3d_search_csv(csv_path):
            return FR3DConverter._convert_fr3d_search_csv(csv_path)
        
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
    def _convert_alignment_log(file_path: str, motif_type: str,
                               apply_filters: bool, custom_pvalues: dict) -> Dict[str, List[MotifInstanceSimple]]:
        """Parse RMSX alignment-report blocks emitted by the lab dataset."""
        blocks = []
        current = None
        with open(file_path, 'r', encoding='utf-8') as handle:
            for raw_line in handle:
                line = raw_line.strip()
                match = re.search(r"Aligning\s+\S+\s+and\s+([^:]+):([^_]+(?:_[^:]+)*)", line, re.IGNORECASE)
                if match:
                    if current:
                        blocks.append(current)
                    current = {
                        'fragment_id': f"{match.group(1)}:{match.group(2).rstrip(':')}",
                        'score': 0.0,
                        'p_value': 1.0,
                    }
                    continue
                if current is None:
                    continue
                score = re.search(r"Alignment score:\s*([-+0-9.eE]+)", line, re.IGNORECASE)
                if score:
                    current['score'] = float(score.group(1))
                pvalue = re.search(r"P-value:\s*([-+0-9.eE]+)", line, re.IGNORECASE)
                if pvalue:
                    current['p_value'] = float(pvalue.group(1))
            if current:
                blocks.append(current)

        instances = []
        threshold = rmsx_pvalue_threshold(motif_type, custom_pvalues)
        for index, block in enumerate(blocks, start=1):
            pdb_id, chain, ranges = RNAMotifScanXConverter.parse_fragment_id(block['fragment_id'])
            if not ranges or (apply_filters and block['p_value'] > threshold):
                continue
            residues = [
                ('N', number, chain)
                for start, end in ranges
                for number in range(start, end + 1)
            ]
            instances.append(MotifInstanceSimple(
                motif_id=motif_type,
                instance_id=f"RMSX_{block['fragment_id'].replace(':', '_').replace('-', '_')}_{index}",
                residues=residues,
                annotation=f"Score: {block['score']}, P-value: {block['p_value']}",
                metadata={
                    'p_value': block['p_value'],
                    'alignment_score': block['score'],
                    'fragment_id': block['fragment_id'],
                    'pdb_id': pdb_id,
                    'chain': chain,
                },
            ))
        instances.sort(key=lambda item: item.metadata.get('alignment_score', 0.0), reverse=True)
        return {motif_type: instances}
    
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

        with open(file_path, 'r', encoding='utf-8') as probe:
            first_content = probe.read()
        if re.search(r"^#\s+Aligning\s+", first_content, re.MULTILINE):
            return RNAMotifScanXConverter._convert_alignment_log(
                file_path, motif_type, apply_filters, custom_pvalues
            )
        
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
                threshold = rmsx_pvalue_threshold(motif_type, custom_pvalues)
                
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
    


