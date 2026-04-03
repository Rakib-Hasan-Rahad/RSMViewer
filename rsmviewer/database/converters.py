"""
RSMViewer - Format Converters
Provides converters for different motif database formats to unified format.

Supported formats:
- RNA 3D Atlas JSON format
- Rfam Stockholm (SEED) format
- Future: CIF annotations, custom JSON schemas

Each converter takes input in its native format and produces standardized
MotifInstance and MotifType objects.

Author: CBB Lab
Version: 2.0.0
"""

from __future__ import annotations

import json
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .base_provider import MotifInstance, MotifType, ResidueSpec


class BaseConverter(ABC):
    """Abstract base class for format converters."""
    
    @abstractmethod
    def convert_file(self, file_path: Path) -> List[MotifType]:
        """
        Convert a file to list of MotifType objects.
        
        Args:
            file_path: Path to the source file
            
        Returns:
            List of MotifType objects
        """
        pass
    
    @abstractmethod
    def convert_data(self, data: Any, source_info: Dict) -> List[MotifType]:
        """
        Convert raw data to MotifType objects.
        
        Args:
            data: Raw data in native format
            source_info: Metadata about the source
            
        Returns:
            List of MotifType objects
        """
        pass


class AtlasJSONConverter(BaseConverter):
    """
    Converter for RNA 3D Motif Atlas JSON format.
    
    Atlas JSON structure:
    [
        {
            "motif_id": "HL_00317.1",
            "common_name": "...",
            "annotation": "...",
            "bp_signature": "...",
            "annotations": {"instance_id": "annotation", ...},
            "num_instances": N,
            "alignment": {
                "instance_id": {
                    "position": "PDB|Model|Chain|Nuc|ResNum|...",
                    ...
                }
            },
            "num_nucleotides": M,
            "chainbreak": []
        },
        ...
    ]
    """
    
    def convert_file(self, file_path: Path) -> List[MotifType]:
        """Convert Atlas JSON file to MotifType objects."""
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Infer motif type from filename (e.g., hl_4.5.json -> HL)
            motif_type_id = self._extract_motif_type_from_filename(file_path.name)
            
            return self.convert_data(data, {
                'type_id': motif_type_id,
                'file': str(file_path),
                'source': 'RNA 3D Motif Atlas'
            })
            
        except Exception as e:
            print(f"Error converting Atlas file {file_path}: {e}")
            return []
    
    def convert_data(self, data: Any, source_info: Dict) -> List[MotifType]:
        """Convert Atlas JSON data to MotifType objects."""
        if not isinstance(data, list):
            return []
        
        type_id = source_info.get('type_id', 'UNKNOWN')
        instances: List[MotifInstance] = []
        
        for entry in data:
            if not isinstance(entry, dict):
                continue
            
            motif_id = entry.get('motif_id', 'unknown')
            alignment = entry.get('alignment', {})
            annotations = entry.get('annotations', {})
            
            if not isinstance(alignment, dict):
                continue
            
            for instance_id, residue_map in alignment.items():
                if not isinstance(residue_map, dict):
                    continue
                
                pdb_id = self._extract_pdb_id(instance_id)
                if not pdb_id:
                    continue
                
                residues = self._parse_residue_map(residue_map)
                annotation = annotations.get(instance_id, '') if isinstance(annotations, dict) else ''
                
                instance = MotifInstance(
                    instance_id=instance_id,
                    motif_id=motif_id,
                    pdb_id=pdb_id,
                    residues=residues,
                    annotation=annotation,
                    metadata={
                        'bp_signature': entry.get('bp_signature', ''),
                        'num_nucleotides': entry.get('num_nucleotides'),
                        'common_name': entry.get('common_name', ''),
                        'chainbreak': entry.get('chainbreak', []),
                    }
                )
                instances.append(instance)
        
        if not instances:
            return []
        
        motif_type = MotifType(
            type_id=type_id,
            name=self._get_type_name(type_id),
            description=f"{type_id} motifs from RNA 3D Motif Atlas",
            instances=instances,
            source='RNA 3D Motif Atlas',
            metadata={'source_file': source_info.get('file', '')}
        )
        
        return [motif_type]
    
    def _extract_motif_type_from_filename(self, filename: str) -> str:
        """Extract motif type from filename like 'hl_4.5.json'."""
        match = re.match(r'^([a-zA-Z]+\d*)_', filename)
        if match:
            return match.group(1).upper()
        return filename.split('.')[0].upper()
    
    def _extract_pdb_id(self, instance_id: str) -> Optional[str]:
        """Extract PDB ID from Atlas instance ID."""
        # Handle residue spec style (rare)
        if '|' in instance_id:
            head = instance_id.split('|', 1)[0]
            return head.upper() if len(head) == 4 else None
        
        # Typical Atlas instance_id like HL_6SVS_002
        parts = instance_id.split('_')
        for part in parts:
            if len(part) == 4 and part.isalnum():
                return part.upper()
        return None
    
    def _parse_residue_map(self, residue_map: Dict[str, str]) -> List[ResidueSpec]:
        """Parse residue map to list of ResidueSpec."""
        result = []
        
        # Sort by position index
        sorted_items = sorted(
            residue_map.items(),
            key=lambda kv: int(kv[0]) if str(kv[0]).isdigit() else 10**9
        )
        
        for _, spec in sorted_items:
            residue = ResidueSpec.from_atlas_spec(spec)
            if residue:
                result.append(residue)
        
        return result
    
    def _get_type_name(self, type_id: str) -> str:
        """Get human-readable name for motif type."""
        names = {
            'HL': 'Hairpin Loops',
            'IL': 'Internal Loops',
            'J3': '3-Way Junctions',
            'J4': '4-Way Junctions',
            'J5': '5-Way Junctions',
            'J6': '6-Way Junctions',
            'J7': '7-Way Junctions',
        }
        return names.get(type_id, type_id)


class StockholmConverter(BaseConverter):
    """
    Converter for Rfam Stockholm (SEED) format.
    
    Stockholm format contains:
    - Metadata lines starting with #=GF (Global Feature)
    - Sequence alignments with structure annotations
    - Column annotations starting with #=GC (Global Column)
    
    Structure annotation encodes:
    - Paired positions: < and > or ( and )
    - Unpaired: . or :
    
    This converter extracts PDB mappings from sequence IDs like:
    - 3OWI_A/41-61 (PDB_Chain/start-end)
    """
    
    def convert_file(self, file_path: Path) -> List[MotifType]:
        """Convert Stockholm SEED file to MotifType objects."""
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Get motif name from parent directory or file
            motif_name = file_path.parent.name
            
            return self.convert_data(content, {
                'type_id': self._normalize_type_id(motif_name),
                'name': motif_name,
                'file': str(file_path),
                'source': 'Rfam'
            })
            
        except Exception as e:
            print(f"Error converting Stockholm file {file_path}: {e}")
            return []
    
    def convert_data(self, data: Any, source_info: Dict) -> List[MotifType]:
        """Convert Stockholm content to MotifType objects."""
        if not isinstance(data, str):
            return []
        
        type_id = source_info.get('type_id', 'UNKNOWN')
        type_name = source_info.get('name', type_id)
        
        lines = data.strip().split('\n')
        metadata = self._parse_metadata(lines)
        sequences = self._parse_sequences(lines)
        
        instances: List[MotifInstance] = []
        
        for seq_id, seq_data in sequences.items():
            pdb_info = self._parse_sequence_id(seq_id)
            if not pdb_info:
                continue
            
            pdb_id, chain, start, end = pdb_info
            
            residues = self._generate_residues(
                chain, start, end, 
                seq_data.get('sequence', ''),
                seq_data.get('structure', '')
            )
            
            instance = MotifInstance(
                instance_id=seq_id,
                motif_id=type_id,
                pdb_id=pdb_id,
                residues=residues,
                annotation=metadata.get('DE', ''),  # Description
                metadata={
                    'sequence': seq_data.get('sequence', ''),
                    'structure': seq_data.get('structure', ''),
                    'rfam_id': metadata.get('AC', ''),
                    'rfam_name': metadata.get('ID', type_name),
                }
            )
            instances.append(instance)
        
        if not instances:
            return []
        
        motif_type = MotifType(
            type_id=type_id,
            name=type_name,
            description=metadata.get('DE', f'{type_name} motif from Rfam'),
            instances=instances,
            source='Rfam',
            metadata={
                'rfam_accession': metadata.get('AC', ''),
                'source_file': source_info.get('file', ''),
                'references': self._extract_references(metadata),
            }
        )
        
        return [motif_type]
    
    def _normalize_type_id(self, name: str) -> str:
        """Normalize motif name to type ID."""
        # Replace spaces and special chars with underscores
        normalized = re.sub(r'[^a-zA-Z0-9]+', '_', name)
        return normalized.upper()
    
    def _parse_metadata(self, lines: List[str]) -> Dict[str, str]:
        """Parse Stockholm metadata lines (#=GF)."""
        metadata: Dict[str, str] = {}
        
        for line in lines:
            if line.startswith('#=GF'):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    key = parts[1]
                    value = parts[2]
                    if key in metadata:
                        metadata[key] += ' ' + value
                    else:
                        metadata[key] = value
        
        return metadata
    
    def _parse_sequences(self, lines: List[str]) -> Dict[str, Dict[str, str]]:
        """
        Parse sequence and structure annotations.
        
        Returns dict: {seq_id: {'sequence': ..., 'structure': ...}}
        """
        sequences: Dict[str, Dict[str, str]] = {}
        
        for line in lines:
            line = line.strip()
            
            # Skip empty, comments, and metadata lines
            if not line or line.startswith('#') or line.startswith('//'):
                continue
            
            # Parse sequence line: ID/range   sequence
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                sequence = parts[1]
                
                if seq_id not in sequences:
                    sequences[seq_id] = {}
                sequences[seq_id]['sequence'] = sequence
            
            # Structure annotation: #=GR ID SS structure
            if line.startswith('#=GR') and ' SS ' in line:
                match = re.match(r'#=GR\s+(\S+)\s+SS\s+(\S+)', line)
                if match:
                    seq_id = match.group(1)
                    structure = match.group(2)
                    if seq_id not in sequences:
                        sequences[seq_id] = {}
                    sequences[seq_id]['structure'] = structure
        
        return sequences
    
    def _parse_sequence_id(self, seq_id: str) -> Optional[Tuple[str, str, int, int]]:
        """
        Parse sequence ID to extract PDB info.
        
        Format: PDBID_CHAIN/start-end (e.g., 3OWI_A/41-61)
        
        Returns: (pdb_id, chain, start, end) or None
        """
        # Pattern: XXXX_Y/N-M where XXXX is PDB ID, Y is chain, N-M is range
        match = re.match(r'(\w{4})_(\w+)/(\d+)-(\d+)', seq_id)
        if match:
            pdb_id = match.group(1).upper()
            chain = match.group(2)
            start = int(match.group(3))
            end = int(match.group(4))
            return (pdb_id, chain, start, end)
        
        # Alternative: just PDB/range (no chain specified)
        match = re.match(r'(\w{4})/(\d+)-(\d+)', seq_id)
        if match:
            pdb_id = match.group(1).upper()
            start = int(match.group(2))
            end = int(match.group(3))
            # Infer default chain based on PDB and residue range
            chain = self._infer_rna_chain(pdb_id, start)
            return (pdb_id, chain, start, end)
        
        return None
    
    def _infer_rna_chain(self, pdb_id: str, start_residue: int) -> str:
        """
        Infer RNA chain for PDBs without explicit chain in SEED file.
        
        Many ribosomal structures use numeric chain IDs for RNA.
        For these structures, we default to chain '0' (the main 23S rRNA)
        since that's where most motifs are located.
        """
        # Known ribosomal structures with numeric RNA chains
        # These use chain '0' for 23S rRNA and '9' for 5S rRNA
        # Default to '0' since most motifs are in the 23S rRNA
        ribosome_pdbs = {
            '1S72', '1FFK', '1NKW', '1S1I', '1J5E', '1GIY', '1C2W', '1YHQ',
            '2GYA', '3CPW', '2AW4', '1J5A', '1C2X', '1VQ6', '1VQ8', '1VQO',
            '1VQP', '1VQ4', '1VQ5', '1VQ7', '1VQ9', '1YIT', '1YIJ', '1YHQ',
            '2OTJ', '2OTL', '3CC2', '3CC4', '3CC7', '3CCE', '3CCJ', '3CCL',
            '3CCM', '3CCQ', '3CCR', '3CCS', '3CCU', '3CD6', '3CMA', '3CME',
        }
        
        if pdb_id.upper() in ribosome_pdbs:
            # Default to chain 0 (23S rRNA) for ribosomal structures
            # The 5S rRNA (chain 9) entries typically have explicit chain in SEED
            return '0'
        
        # Default to 'A' for other structures
        return 'A'
    
    def _generate_residues(self, chain: str, start: int, end: int,
                          sequence: str, structure: str) -> List[ResidueSpec]:
        """
        Generate ResidueSpec list from range and optional sequence.
        
        Args:
            chain: Chain identifier
            start: Start residue number
            end: End residue number
            sequence: Aligned sequence (may contain gaps)
            structure: Secondary structure annotation
        
        Returns:
            List of ResidueSpec for actual (non-gap) positions
        """
        residues = []
        
        if sequence:
            # Use sequence to determine actual positions (skip gaps)
            res_num = start
            for i, char in enumerate(sequence):
                if char not in '.~-':  # Not a gap
                    nuc = char.upper() if char.isalpha() else ''
                    residues.append(ResidueSpec(
                        chain=chain,
                        residue_number=res_num,
                        nucleotide=nuc
                    ))
                    res_num += 1
        else:
            # No sequence, just use range
            for res_num in range(start, end + 1):
                residues.append(ResidueSpec(
                    chain=chain,
                    residue_number=res_num
                ))
        
        return residues
    
    def _extract_references(self, metadata: Dict[str, str]) -> List[Dict]:
        """Extract publication references from metadata."""
        refs = []
        for key in ['RM', 'RT', 'RA', 'RL']:
            if key in metadata:
                refs.append({key: metadata[key]})
        return refs


# Factory function for getting appropriate converter
def get_converter(format_type: str) -> BaseConverter:
    """
    Get appropriate converter for a format type.
    
    Args:
        format_type: Format identifier ('atlas_json', 'stockholm', etc.)
        
    Returns:
        Converter instance
    """
    converters = {
        'atlas_json': AtlasJSONConverter,
        'atlas': AtlasJSONConverter,
        'stockholm': StockholmConverter,
        'rfam': StockholmConverter,
        'seed': StockholmConverter,
    }
    
    converter_class = converters.get(format_type.lower())
    if converter_class:
        return converter_class()
    
    raise ValueError(f"Unknown format type: {format_type}")
