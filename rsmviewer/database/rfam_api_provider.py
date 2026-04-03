"""
RSMViewer - Rfam API Provider
Fetches motif data from the Rfam REST API.

API Endpoint: https://rfam.org/family/{RF_ID}/structures?content-type=application/json
Returns JSON data with PDB structure mappings for RNA families.

This provider enables visualization of named RNA motifs (GNRA, K-turn, T-loop, etc.)
for structures that have Rfam annotations.

Author: Structural Biology Lab
Version: 2.0.0
"""

from __future__ import annotations

import json
import ssl
import urllib.request
import urllib.error
from typing import Dict, List, Optional, Set

from .base_provider import (
    BaseProvider,
    DatabaseInfo,
    DatabaseSourceType,
    MotifInstance,
    MotifType,
    ResidueSpec,
)


class RfamAPIProvider(BaseProvider):
    """
    Provider that fetches RNA motif data from Rfam API.
    
    Supports named RNA motifs like:
    - GNRA tetraloop
    - UNCG tetraloop
    - K-turn
    - T-loop
    - C-loop
    - U-turn
    - And more...
    
    The Rfam database focuses on RNA families identified by sequence/structure
    conservation, mapped to PDB structures.
    """
    
    # Base URL for Rfam API
    API_BASE_URL = "https://rfam.org"
    
    # Timeout for API requests (seconds)
    REQUEST_TIMEOUT = 30
    
    # Mapping of Rfam motif IDs to readable names
    # These are the main structural motifs in Rfam
    MOTIF_IDS = {
        'RM00008': {'name': 'GNRA tetraloop', 'short': 'GNRA'},
        'RM00029': {'name': 'UNCG tetraloop', 'short': 'UNCG'},
        'RM00010': {'name': 'Kink-turn', 'short': 'K-turn'},
        'RM00024': {'name': 'T-loop', 'short': 'T-loop'},
        'RM00003': {'name': 'C-loop', 'short': 'C-loop'},
        'RM00030': {'name': 'U-turn', 'short': 'U-turn'},
        'RM00021': {'name': 'Tandem GA/AG', 'short': 'tandem-GA'},
        'RM00028': {'name': 'UMAC tetraloop', 'short': 'UMAC'},
        'RM00007': {'name': 'Splicing Domain V', 'short': 'Domain-V'},
        'RM00022': {'name': 'Rho terminator 1', 'short': 'Terminator1'},
        'RM00023': {'name': 'Rho terminator 2', 'short': 'Terminator2'},
        'RM00005': {'name': 'CsrA/RsmA binding', 'short': 'CsrA_binding'},
    }
    
    def __init__(self, cache_manager=None):
        """
        Initialize the Rfam API provider.
        
        Args:
            cache_manager: Optional cache manager for storing API responses
        """
        # Create a fake path for the base class (API doesn't use local files)
        super().__init__("api://rfam.org")
        
        self._info = DatabaseInfo(
            id="rfam_api",
            name="Rfam (Online)",
            version="API",
            description="Live data from Rfam database - named RNA motifs",
            source_type=DatabaseSourceType.API,
        )
        self.cache_manager = cache_manager
        self._fetched_pdbs: Set[str] = set()
        self._motif_pdb_cache: Dict[str, Dict] = {}  # Cache of motif->PDB mappings
        self._pdb_motif_cache: Dict[str, Dict[str, List[MotifInstance]]] = {}  # pdb -> motifs
    
    @property
    def info(self) -> DatabaseInfo:
        """Get database metadata."""
        return self._info
    
    def initialize(self) -> bool:
        """Initialize the provider (always succeeds for API provider)."""
        self._initialized = True
        return True
    
    def get_available_motif_types(self) -> List[str]:
        """Get list of motif type IDs supported."""
        return [info['short'] for info in self.MOTIF_IDS.values()]
    
    def get_available_pdb_ids(self) -> List[str]:
        """Get list of PDB IDs that have been successfully fetched."""
        return list(self._fetched_pdbs)
    
    def get_motif_type(self, motif_type_id: str) -> Optional[MotifType]:
        """
        Get information about a motif type.
        
        Args:
            motif_type_id: Motif type ID (e.g., 'GNRA', 'K-turn')
            
        Returns:
            MotifType object or None
        """
        # Find by short name
        for rm_id, info in self.MOTIF_IDS.items():
            if info['short'] == motif_type_id:
                return MotifType(
                    type_id=motif_type_id,
                    name=info['name'],
                    description=f"Rfam motif: {info['name']}",
                    metadata={'rfam_id': rm_id}
                )
        return None
    
    def get_motifs_for_pdb(self, pdb_id: str) -> Dict[str, List[MotifInstance]]:
        """
        Get all Rfam motifs for a PDB structure.
        
        This queries each known Rfam motif family to find matches.
        Results are cached.
        
        Args:
            pdb_id: PDB structure ID
            
        Returns:
            Dict mapping motif type IDs to lists of MotifInstances
        """
        pdb_id = pdb_id.strip().upper()
        
        # Check internal cache
        if pdb_id in self._pdb_motif_cache:
            return self._pdb_motif_cache[pdb_id]
        
        # Check file cache
        if self.cache_manager:
            cached = self.cache_manager.get_cached_motifs(pdb_id, "rfam_api")
            if cached is not None:
                self._fetched_pdbs.add(pdb_id)
                self._pdb_motif_cache[pdb_id] = cached
                return cached
        
        result: Dict[str, List[MotifInstance]] = {}
        
        # Query each motif family
        for rm_id, info in self.MOTIF_IDS.items():
            try:
                instances = self._get_motif_instances_for_pdb(pdb_id, rm_id, info)
                if instances:
                    result[info['short']] = instances
            except Exception as e:
                # Silently continue - many motif types won't have matches for a given PDB
                # This is expected behavior, not an error
                continue
        
        # Cache the results
        if self.cache_manager and result:
            self.cache_manager.cache_motifs(pdb_id, "rfam_api", result)
        
        if result:
            self._fetched_pdbs.add(pdb_id)
            self._pdb_motif_cache[pdb_id] = result
        
        return result
    
    def _get_motif_instances_for_pdb(
        self, pdb_id: str, rfam_motif_id: str, motif_info: Dict
    ) -> List[MotifInstance]:
        """
        Get instances of a specific Rfam motif in a PDB.
        
        Args:
            pdb_id: PDB ID
            rfam_motif_id: Rfam motif ID (e.g., 'RM00008')
            motif_info: Motif metadata dict
            
        Returns:
            List of MotifInstance objects
        """
        # First, get the PDB mappings for this motif family
        pdb_mappings = self._get_pdb_mappings_for_motif(rfam_motif_id)
        
        if pdb_id not in pdb_mappings:
            return []
        
        instances = []
        mapping_data = pdb_mappings[pdb_id]
        
        # Create MotifInstance from mapping data
        for idx, mapping in enumerate(mapping_data):
            instance_id = f"{motif_info['short']}_{pdb_id}_{idx + 1:03d}"
            
            # Parse residue information from mapping
            residues = self._parse_rfam_residues(mapping, pdb_id)
            
            instances.append(MotifInstance(
                instance_id=instance_id,
                motif_id=motif_info['short'],
                pdb_id=pdb_id,
                residues=residues,
                annotation=motif_info['name'],
                metadata={
                    'source': 'rfam_api',
                    'rfam_id': rfam_motif_id,
                    'raw_mapping': mapping
                }
            ))
        
        return instances
    
    def _get_pdb_mappings_for_motif(self, rfam_motif_id: str) -> Dict[str, List]:
        """
        Get PDB structure mappings for a Rfam motif.
        
        Uses Rfam /motif/{id} or /family endpoint.
        
        Args:
            rfam_motif_id: Rfam motif ID
            
        Returns:
            Dict mapping PDB IDs to lists of mapping data
        """
        # Check in-memory cache
        if rfam_motif_id in self._motif_pdb_cache:
            return self._motif_pdb_cache[rfam_motif_id]
        
        url = f"{self.API_BASE_URL}/motif/{rfam_motif_id}?content-type=application/json"
        
        try:
            request = urllib.request.Request(
                url,
                headers={
                    'User-Agent': 'RSMViewer/2.0',
                    'Accept': 'application/json',
                }
            )
            
            # Create SSL context that doesn't verify certificates
            # This handles macOS certificate issues
            ssl_context = ssl.create_default_context()
            ssl_context.check_hostname = False
            ssl_context.verify_mode = ssl.CERT_NONE
            
            with urllib.request.urlopen(request, timeout=self.REQUEST_TIMEOUT, context=ssl_context) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode('utf-8'))
                    
                    # Parse the response to extract PDB mappings
                    pdb_mappings = self._parse_rfam_motif_response(data)
                    
                    # Cache the result
                    self._motif_pdb_cache[rfam_motif_id] = pdb_mappings
                    return pdb_mappings
                    
        except urllib.error.HTTPError as e:
            # 404 is expected for motifs not in Rfam - silently return empty
            pass
        except urllib.error.URLError as e:
            # Network errors - silently return empty
            pass
        except Exception as e:
            # Other errors - silently return empty
            pass
        
        return {}
    
    def _parse_rfam_motif_response(self, data: Dict) -> Dict[str, List]:
        """
        Parse Rfam API motif response to extract PDB mappings.
        
        The Rfam API returns structure mappings in various formats.
        We extract PDB IDs and residue ranges.
        
        Args:
            data: JSON response from Rfam API
            
        Returns:
            Dict mapping PDB IDs to lists of residue mapping data
        """
        pdb_mappings: Dict[str, List] = {}
        
        # Try different response formats
        # Rfam may return 'pdb' or 'structures' field
        structures = data.get('structures', data.get('pdb', []))
        
        if isinstance(structures, list):
            for struct in structures:
                if isinstance(struct, dict):
                    pdb_id = struct.get('pdb_id', struct.get('pdb', '')).upper()
                    if pdb_id and len(pdb_id) == 4:
                        if pdb_id not in pdb_mappings:
                            pdb_mappings[pdb_id] = []
                        pdb_mappings[pdb_id].append(struct)
        
        return pdb_mappings
    
    def _parse_rfam_residues(self, mapping: Dict, pdb_id: str) -> List[ResidueSpec]:
        """
        Parse residue information from Rfam mapping data.
        
        Rfam provides chain and residue range information.
        
        Args:
            mapping: Mapping dict from Rfam
            pdb_id: PDB ID for context
            
        Returns:
            List of ResidueSpec objects
        """
        residues = []
        
        chain = mapping.get('chain', mapping.get('auth_asym_id', 'A'))
        seq_start = mapping.get('seq_start', mapping.get('pdb_start', 1))
        seq_end = mapping.get('seq_end', mapping.get('pdb_end', seq_start))
        
        try:
            start = int(seq_start)
            end = int(seq_end)
            
            for res_num in range(start, end + 1):
                residues.append(ResidueSpec(
                    chain=str(chain),
                    residue_number=res_num,
                    nucleotide='',  # Rfam doesn't always provide nucleotide type
                ))
        except (ValueError, TypeError):
            pass
        
        return residues
    
    def has_pdb(self, pdb_id: str) -> bool:
        """
        Check if a PDB has any Rfam motif annotations.
        
        Args:
            pdb_id: PDB ID to check
            
        Returns:
            True if PDB has Rfam motif data
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
        if pdb_id in self._pdb_motif_cache:
            motifs = self._pdb_motif_cache[pdb_id]
        else:
            motifs = self.get_motifs_for_pdb(pdb_id)
            self._pdb_motif_cache[pdb_id] = motifs
        
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
            motif_type_id: Motif type short name (e.g., 'GNRA', 'K-turn')
            
        Returns:
            List of MotifInstance objects
        """
        all_motifs = self.get_motifs_for_pdb(pdb_id)
        return all_motifs.get(motif_type_id, [])
