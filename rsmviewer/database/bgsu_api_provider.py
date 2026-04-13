"""
RSMViewer - BGSU RNA 3D Hub Hybrid Provider
Fetches motif data from BGSU RNA 3D Hub using HTML scraping (primary) 
and CSV API (fallback).

Primary Method: HTML web scraping from https://rna.bgsu.edu/rna3dhub/pdb/{PDB_ID}/motifs
  - Provides semantic annotations (K-turn, GNRA, C-loop, etc.)
  - Provides motif group IDs
  
Fallback Method: CSV API from https://rna.bgsu.edu/rna3dhub/loops/download/{PDB_ID}
  - Used if HTML parsing fails
  - Provides basic loop data only

This provider enables visualization of RNA motifs for ANY PDB structure
in the RNA 3D Hub database (~3000+ RNA structures), not just those
bundled locally.

Author: CBB Lab
Version: 1.0.0
"""

from __future__ import annotations

import re
import ssl
import urllib.request
import urllib.error
from html.parser import HTMLParser
from typing import Dict, List, Optional, Set
from collections import defaultdict

from .base_provider import (
    BaseProvider,
    DatabaseInfo,
    DatabaseSourceType,
    MotifInstance,
    MotifType,
    ResidueSpec,
)


class BGSUHTMLParser(HTMLParser):
    """HTML parser for BGSU motif annotation table"""
    
    def __init__(self):
        super().__init__()
        self.annotations = {}
        self.current_row = []
        self.current_cell = []
        self.in_tbody = False
        self.in_tr = False
        self.in_td = False
        
    def handle_starttag(self, tag, attrs):
        if tag == 'tbody':
            self.in_tbody = True
        elif tag == 'tr' and self.in_tbody:
            self.in_tr = True
            self.current_row = []
        elif tag == 'td' and self.in_tr:
            self.in_td = True
            self.current_cell = []
    
    def handle_endtag(self, tag):
        if tag == 'tbody':
            self.in_tbody = False
        elif tag == 'tr' and self.in_tr:
            self.in_tr = False
            self._process_row()
        elif tag == 'td' and self.in_td:
            self.in_td = False
            cell_text = ''.join(self.current_cell).strip()
            self.current_row.append(cell_text)
            self.current_cell = []
    
    def handle_data(self, data):
        if self.in_td:
            self.current_cell.append(data)
    
    def _process_row(self):
        """Process a completed table row.
        
        HTML cell 3 contains BOTH the annotation text and the motif group link,
        separated by <br><a>. The HTMLParser concatenates all text into one string,
        e.g. "GNRAHL_34789.4" or "No text annotationHL_75660.8".
        We extract the motif group ID via regex and strip it from the annotation.
        """
        if len(self.current_row) >= 4:
            # Extract loop ID from second cell
            loop_id_match = re.search(r'((?:HL|IL|J\d)_\S+_\d{3,4})', self.current_row[1])
            if loop_id_match:
                loop_id = loop_id_match.group(1)
                residues = self.current_row[2].strip()
                raw_col3 = self.current_row[3].strip()
                
                # Extract motif group ID from column 3 (embedded via <a> tag)
                # Pattern: HL_34789.4, IL_58355.2, J3_01004.1, etc.
                motif_group = ""
                annotation_text = raw_col3
                mg_match = re.search(r'((?:IL|HL|J\d)_\d+\.\d+)', raw_col3)
                if mg_match:
                    motif_group = mg_match.group(1)
                    # Strip the motif group ID from annotation text
                    annotation_text = raw_col3[:mg_match.start()].strip()
                
                # Fallback: check column 4 if motif group not found in column 3
                if not motif_group and len(self.current_row) >= 5:
                    mg_match4 = re.search(r'((?:IL|HL|J\d)_\d+\.\d+)', self.current_row[4])
                    if mg_match4:
                        motif_group = mg_match4.group(1)
                
                loop_type = loop_id.split('_')[0]
                
                self.annotations[loop_id] = {
                    'loop_id': loop_id,
                    'residues': residues,
                    'annotation': annotation_text if annotation_text and annotation_text != 'No text annotation' else None,
                    'motif_group': motif_group,
                    'type': loop_type
                }


class BGSUAPIProvider(BaseProvider):
    """
    Provider that fetches RNA motif data from BGSU RNA 3D Hub.
    
    Uses hybrid approach:
    1. PRIMARY: HTML scraping for semantic annotations (K-turn, GNRA, C-loop, etc.)
    2. FALLBACK: CSV API for basic loop data if HTML parsing fails
    
    Supports:
    - Semantic annotations: Kink-turn, GNRA, C-loop, Sarcin-Ricin, UAA/GAN, etc.
    - Loop types (fallback): HL, IL, J3, J4, J5, J6, J7, J8
    """
    
    # Base URLs
    HTML_BASE_URL = "https://rna.bgsu.edu/rna3dhub/pdb"
    API_BASE_URL = "https://rna.bgsu.edu/rna3dhub/loops/download"
    
    # Timeout for requests (seconds)
    REQUEST_TIMEOUT = 30
    
    # Loop type names for fallback (when no semantic annotation)
    LOOP_TYPE_NAMES = {
        'HL': 'Hairpin Loop (HL)',
        'IL': 'Internal Loop (IL)',
        'J3': '3-way Junction (J3)',
        'J4': '4-way Junction (J4)',
        'J5': '5-way Junction (J5)',
        'J6': '6-way Junction (J6)',
        'J7': '7-way Junction (J7)',
        'J8': '8-way Junction (J8)',
    }
    
    # Semantic annotation patterns to recognize
    SEMANTIC_PATTERNS = {
        'Kink-turn': ['Kink-turn', 'kink turn'],
        'C-loop': ['C-loop', 'mini C-loop'],
        'GNRA': ['GNRA'],
        'Sarcin-Ricin': ['Sarcin', 'Sarcin-Ricin'],
        'E-loop': ['E-loop'],
        'UAA/GAN': ['UAA/GAN'],
        'Triple sheared': ['Triple sheared'],
        'Major groove platform': ['Major groove platform'],
        'Minor groove platform': ['Minor groove platform'],
        'Tetraloop': ['tetraloop', 'Tetraloop'],
        'Bulged': ['bulged', 'Bulged'],
        'UNCG': ['UNCG'],
        'T-loop': ['T-loop'],
        'Pseudoknot': ['Pseudoknot', 'pseudoknot'],
    }
    
    def __init__(self, cache_manager=None):
        """
        Initialize the BGSU hybrid provider.
        
        Args:
            cache_manager: Optional cache manager for storing responses
        """
        super().__init__("api://bgsu.rna3dhub")
        
        self._info = DatabaseInfo(
            id="bgsu_api",
            name="BGSU RNA 3D Hub (Online)",
            version="3.0 Hybrid",
            description="Live data from BGSU RNA 3D Hub with semantic annotations - supports ~3000+ RNA structures",
            source_type=DatabaseSourceType.API,
        )
        self.cache_manager = cache_manager
        self._available_motif_types = []  # Will be populated dynamically
        self._fetched_pdbs: Set[str] = set()
        self._motif_cache: Dict[str, Dict[str, List[MotifInstance]]] = {}
        self._annotation_cache: Dict[str, Dict] = {}  # Store raw annotations
    
    @property
    def info(self) -> DatabaseInfo:
        """Get database metadata."""
        return self._info
    
    def initialize(self) -> bool:
        """Initialize the provider (always succeeds for API provider)."""
        self._initialized = True
        return True
    
    def get_available_motif_types(self) -> List[str]:
        """Get list of motif types supported by BGSU RNA 3D Hub."""
        return self._available_motif_types
    
    def get_available_pdb_ids(self) -> List[str]:
        """
        Get list of PDB IDs that have been successfully fetched.
        
        Note: This only returns PDBs we've already fetched. The actual
        BGSU database has ~3000+ structures, but we don't enumerate them.
        """
        return list(self._fetched_pdbs)
    
    def get_motif_type(self, motif_type_id: str) -> Optional[MotifType]:
        """
        Get information about a motif type.
        
        Args:
            motif_type_id: Motif type ID (e.g., 'Kink-turn', 'GNRA', 'HL', 'IL')
            
        Returns:
            MotifType object or None
        """
        # Check if it's a semantic annotation category
        if motif_type_id in self.SEMANTIC_PATTERNS:
            return MotifType(
                type_id=motif_type_id,
                name=motif_type_id,
                description=f"{motif_type_id} - Semantic annotation from BGSU RNA 3D Hub",
            )
        
        # Check if it's a loop type
        if motif_type_id in self.LOOP_TYPE_NAMES:
            return MotifType(
                type_id=motif_type_id,
                name=self.LOOP_TYPE_NAMES[motif_type_id],
                description=f"{self.LOOP_TYPE_NAMES[motif_type_id]} from BGSU RNA 3D Hub",
            )
        
        # For any other category found dynamically
        return MotifType(
            type_id=motif_type_id,
            name=motif_type_id,
            description=f"{motif_type_id} from BGSU RNA 3D Hub",
        )
    
    def get_motifs_for_pdb(self, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """
        Fetch motifs for a PDB structure using hybrid approach.
        
        Strategy:
        1. Fetch CSV API for detailed residue data
        2. Fetch HTML for semantic annotations
        3. Merge annotations into CSV data
        
        Args:
            pdb_id: PDB structure ID
            
        Returns:
            Dict mapping motif type IDs to lists of MotifInstances
        """
        pdb_id = pdb_id.strip().upper()
        
        # Check internal cache first
        if pdb_id in self._motif_cache:
            return self._motif_cache[pdb_id]
        
        # Check file cache if available
        if self.cache_manager:
            cached = self.cache_manager.get_cached_motifs(pdb_id, "bgsu_api")
            if cached is not None:
                self._fetched_pdbs.add(pdb_id)
                self._motif_cache[pdb_id] = cached
                # Update available motif types
                for motif_type in cached.keys():
                    if motif_type not in self._available_motif_types:
                        self._available_motif_types.append(motif_type)
                return cached
        
        # Fetch CSV first (for residue data)
        try:
            csv_data = self._fetch_from_api(pdb_id)
            if csv_data is None:
                return {}
            
            # Try to fetch HTML annotations
            annotations = None
            try:
                annotations = self._fetch_html_annotations(pdb_id)
                if annotations:
                    print(f"Successfully fetched {len(annotations)} annotations from HTML for {pdb_id}")
            except Exception as e:
                print(f"HTML annotation fetch failed for {pdb_id}: {e}")
            
            # Parse CSV with annotations (if available)
            if annotations:
                motifs = self._parse_csv_with_annotations(csv_data, annotations, pdb_id)
            else:
                import logging as _logging
                _lg = _logging.getLogger(__name__)
                _lg.warning(
                    f"Could not retrieve semantic motif annotations from "
                    f"BGSU RNA 3D Hub HTML page for {pdb_id}. "
                    f"Using CSV loop data with generic names instead."
                )
                print(
                    f"\u26a0 Could not retrieve semantic annotations from "
                    f"BGSU RNA 3D Hub for {pdb_id}. Using CSV loop data instead."
                )
                motifs = self._parse_csv_response(csv_data, pdb_id)
            
            # Cache the results
            if self.cache_manager and motifs:
                self.cache_manager.cache_motifs(pdb_id, "bgsu_api", motifs)
            
            if motifs:
                self._fetched_pdbs.add(pdb_id)
                self._motif_cache[pdb_id] = motifs
                
                # Update available motif types
                for motif_type in motifs.keys():
                    if motif_type not in self._available_motif_types:
                        self._available_motif_types.append(motif_type)
            
            return motifs
            
        except Exception as e:
            print(f"Error fetching motifs for {pdb_id} from BGSU: {e}")
            return {}
    
    def _fetch_html_annotations(self, pdb_id: str) -> Optional[Dict[str, Dict]]:
        """
        Fetch semantic annotations from BGSU HTML page.
        
        Args:
            pdb_id: PDB ID to fetch
            
        Returns:
            Dict mapping loop IDs to annotation data, or None if failed
        """
        url = f"{self.HTML_BASE_URL}/{pdb_id}/motifs"
        
        try:
            # Create request with proper headers
            request = urllib.request.Request(
                url,
                headers={
                    'User-Agent': 'RSMViewer/3.0',
                    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
                }
            )
            
            # Create SSL context that doesn't verify certificates
            ssl_context = ssl.create_default_context()
            ssl_context.check_hostname = False
            ssl_context.verify_mode = ssl.CERT_NONE
            
            with urllib.request.urlopen(request, timeout=self.REQUEST_TIMEOUT, context=ssl_context) as response:
                if response.status == 200:
                    html_content = response.read().decode('utf-8')
                    
                    # Parse HTML to extract annotations
                    parser = BGSUHTMLParser()
                    parser.feed(html_content)
                    
                    if parser.annotations:
                        self._annotation_cache[pdb_id] = parser.annotations
                        return parser.annotations
                    else:
                        print(f"No annotations found in HTML for {pdb_id}")
                        return None
                else:
                    print(f"BGSU HTML page returned status {response.status} for {pdb_id}")
                    return None
                    
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(f"PDB {pdb_id} not found on BGSU HTML page")
            else:
                print(f"HTTP error fetching HTML for {pdb_id}: {e.code} {e.reason}")
            return None
            
        except urllib.error.URLError as e:
            print(f"Network error fetching HTML for {pdb_id}: {e.reason}")
            return None
            
        except Exception as e:
            print(f"Error parsing HTML for {pdb_id}: {e}")
            return None
    
    def _parse_csv_with_annotations(self, csv_data: str, annotations: Dict[str, Dict], pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """
        Parse BGSU CSV response and enrich with HTML annotations.
        
        Strategy:
        1. Parse CSV for accurate residue data
        2. Match loop IDs with HTML annotations
        3. Categorize by semantic annotation or fallback to loop type
        
        Args:
            csv_data: Raw CSV string from API
            annotations: Dict from _fetch_html_annotations
            pdb_id: PDB ID for validation
            
        Returns:
            Dict mapping motif categories to lists of MotifInstances
        """
        result: Dict[str, List[MotifInstance]] = defaultdict(list)
        
        # Pattern to match CSV entries: "LOOP_ID","RESIDUES"
        pattern = r'"([^"]+)","([^"]+)"'
        
        for match in re.finditer(pattern, csv_data):
            loop_id = match.group(1).strip()
            residue_specs = match.group(2).strip()
            
            # Parse loop ID to get motif type
            # Format: {TYPE}_{PDB}_{NUMBER} e.g., HL_4V9F_001
            parts = loop_id.split('_')
            if len(parts) < 2:
                continue
            
            loop_type = parts[0]  # HL, IL, J3, etc.
            
            # Get annotation from HTML if available
            annotation = None
            motif_group = ""
            if loop_id in annotations:
                annotation = annotations[loop_id].get('annotation')
                motif_group = annotations[loop_id].get('motif_group', '')
            
            # Determine category (semantic annotation or loop type fallback)
            category = self._categorize_motif(annotation, loop_type)
            
            # Parse residues from CSV (accurate data)
            try:
                residues = self._parse_residue_specs(residue_specs, pdb_id)
                if not residues:
                    continue
            except Exception as e:
                print(f"Error parsing residues for {loop_id}: {e}")
                continue
            
            # Create MotifInstance
            instance = MotifInstance(
                instance_id=loop_id,
                motif_id=category,
                pdb_id=pdb_id,
                residues=residues,
                annotation=annotation if annotation else category,
                metadata={
                    'source': 'bgsu_hybrid',
                    'motif_group': motif_group,
                    'loop_type': loop_type,
                },
            )
            
            result[category].append(instance)
        
        # Convert defaultdict to regular dict
        return dict(result)
    
    def _parse_with_annotations(self, annotations: Dict[str, Dict], pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """
        DEPRECATED: Parse annotations into MotifInstance objects with smart categorization.
        
        This method is not used anymore because HTML doesn't provide detailed residue information.
        Use _parse_csv_with_annotations instead.
        """
        raise NotImplementedError("Use _parse_csv_with_annotations instead")
    
    def _categorize_motif(self, annotation: Optional[str], loop_type: str) -> str:
        """
        Determine motif category from semantic annotation or loop type.
        
        Args:
            annotation: Semantic annotation from HTML (e.g., "Kink-turn")
            loop_type: Loop type prefix (e.g., "HL", "IL")
            
        Returns:
            Category string for grouping
        """
        # Priority 1: Match semantic annotation patterns
        if annotation:
            annotation_lower = annotation.lower()
            
            for category, patterns in self.SEMANTIC_PATTERNS.items():
                for pattern in patterns:
                    if pattern.lower() in annotation_lower:
                        return category
            
            # If has annotation but no pattern match, use annotation as-is
            return annotation
        
        # Priority 2: Fallback to loop type
        return self.LOOP_TYPE_NAMES.get(loop_type, loop_type)
    
    def _fetch_from_api(self, pdb_id: str) -> Optional[str]:
        """
        Fetch raw CSV data from BGSU API (fallback method).
        
        Args:
            pdb_id: PDB ID to fetch
            
        Returns:
            Raw CSV string or None if failed
        """
        url = f"{self.API_BASE_URL}/{pdb_id}"
        
        try:
            # Create request with proper headers
            request = urllib.request.Request(
                url,
                headers={
                    'User-Agent': 'RSMViewer/2.0',
                    'Accept': 'text/csv, text/plain, */*',
                }
            )
            
            # Create SSL context that doesn't verify certificates
            # This handles macOS certificate issues
            ssl_context = ssl.create_default_context()
            ssl_context.check_hostname = False
            ssl_context.verify_mode = ssl.CERT_NONE
            
            with urllib.request.urlopen(request, timeout=self.REQUEST_TIMEOUT, context=ssl_context) as response:
                if response.status == 200:
                    return response.read().decode('utf-8')
                else:
                    print(f"BGSU API returned status {response.status} for {pdb_id}")
                    return None
                    
        except urllib.error.HTTPError as e:
            if e.code == 404:
                # PDB not found in RNA 3D Hub - this is normal for non-RNA structures
                print(f"PDB {pdb_id} not found in RNA 3D Hub database")
            else:
                print(f"HTTP error fetching {pdb_id}: {e.code} {e.reason}")
            return None
            
        except urllib.error.URLError as e:
            print(f"Network error fetching {pdb_id}: {e.reason}")
            return None
            
        except Exception as e:
            print(f"Error fetching {pdb_id} from BGSU API: {e}")
            return None
    
    def _parse_csv_response(self, csv_data: str, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """
        Parse BGSU CSV response into MotifInstance objects (fallback method).
        
        CSV Format:
            "HL_4V9F_001","4V9F|1|0|U|55,4V9F|1|0|G|56,..."
            
        Args:
            csv_data: Raw CSV string from API
            pdb_id: PDB ID for validation
            
        Returns:
            Dict mapping motif categories to lists of MotifInstances
        """
        result: Dict[str, List[MotifInstance]] = defaultdict(list)
        
        # Pattern to match CSV entries: "LOOP_ID","RESIDUES"
        # Handle various whitespace and newline formats
        pattern = r'"([^"]+)","([^"]+)"'
        
        for match in re.finditer(pattern, csv_data):
            loop_id = match.group(1).strip()
            residue_specs = match.group(2).strip()
            
            # Parse loop ID to get motif type
            # Format: {TYPE}_{PDB}_{NUMBER} e.g., HL_4V9F_001
            parts = loop_id.split('_')
            if len(parts) < 2:
                continue
            
            loop_type = parts[0]  # HL, IL, J3, etc.
            
            # Get category name using loop type fallback
            category = self.LOOP_TYPE_NAMES.get(loop_type, loop_type)
            
            # Parse residues
            try:
                residues = self._parse_residue_specs(residue_specs, pdb_id)
                if not residues:
                    continue
            except Exception as e:
                print(f"Error parsing residues for {loop_id}: {e}")
                continue
            
            # Create MotifInstance
            instance = MotifInstance(
                instance_id=loop_id,
                motif_id=category,
                pdb_id=pdb_id,
                residues=residues,
                annotation=category,
                metadata={'source': 'bgsu_csv', 'loop_type': loop_type},
            )
            
            result[category].append(instance)
        
        # Convert defaultdict to regular dict
        return dict(result)
    
    def _parse_residue_specs(self, specs_str: str, pdb_id: str = "") -> List[ResidueSpec]:
        """
        Parse comma-separated residue specifications.
        
        Format: PDB|Model|Chain|Nucleotide|ResNum
        Example: 4V9F|1|0|U|55,4V9F|1|0|G|56
        
        Args:
            specs_str: Comma-separated residue specifications
            pdb_id: PDB ID for validation (optional)
            
        Returns:
            List of ResidueSpec objects
        """
        residues = []
        
        for spec in specs_str.split(','):
            spec = spec.strip()
            if not spec:
                continue
            
            parts = spec.split('|')
            if len(parts) < 5:
                continue
            
            try:
                # PDB|Model|Chain|Nucleotide|ResNum
                model = int(parts[1]) if parts[1].isdigit() else 1
                chain = parts[2]
                nucleotide = parts[3]
                res_num = int(parts[4])
                
                residues.append(ResidueSpec(
                    chain=chain,
                    residue_number=res_num,
                    nucleotide=nucleotide,
                    model=model,
                ))
            except (ValueError, IndexError):
                continue
        
        return residues
    
    def has_pdb(self, pdb_id: str) -> bool:
        """
        Check if a PDB has motif data in BGSU.
        
        This actually fetches data to check, but results are cached.
        
        Args:
            pdb_id: PDB ID to check
            
        Returns:
            True if PDB has motif data
        """
        motifs = self.get_motifs_for_pdb(pdb_id)
        return len(motifs) > 0
    
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
        pdb_id = pdb_id.upper()
        
        # Check cache first
        if pdb_id in self._motif_cache:
            motifs = self._motif_cache[pdb_id]
        else:
            motifs = self.get_motifs_for_pdb(pdb_id)
            self._motif_cache[pdb_id] = motifs
        
        # Find the specific instance
        if motif_type not in motifs:
            return []
        
        for instance in motifs[motif_type]:
            if instance.instance_id == instance_id:
                return instance.residues
        
        return []
    
    def get_motif_instances_for_pdb(self, pdb_id: str, motif_type_id: str) -> List[MotifInstance]:
        """
        Get all instances of a specific motif type in a PDB.
        
        Args:
            pdb_id: PDB ID
            motif_type_id: Motif type (e.g., 'HL', 'IL')
            
        Returns:
            List of MotifInstance objects
        """
        all_motifs = self.get_motifs_for_pdb(pdb_id)
        return all_motifs.get(motif_type_id, [])
