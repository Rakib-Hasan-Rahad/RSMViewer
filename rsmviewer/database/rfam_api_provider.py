"""
RSMViewer - Rfam API Provider
Fetches motif data from the Rfam REST API.

Downloads SEED alignments in Stockholm format from:
  https://rfam.org/motif/{RM_ID}/alignment?acc={RM_ID}&alnType=seed&nseLabels=0&format=stockholm

The Stockholm data is parsed with the same StockholmConverter used by the
local Rfam provider (Source 2), producing identical residue-level annotations.

Author: Structural Biology Lab
Version: 1.0.0
"""

from __future__ import annotations

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
from .converters import StockholmConverter


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
    # Complete set of all 34 Rfam structural motifs
    MOTIF_IDS = {
        'RM00001': {'name': 'ANYA', 'short': 'ANYA'},
        'RM00002': {'name': 'AUF1 binding', 'short': 'AUF1_binding'},
        'RM00003': {'name': 'C-loop', 'short': 'C-loop'},
        'RM00004': {'name': 'CRC binding', 'short': 'CRC_binding'},
        'RM00005': {'name': 'CsrA/RsmA binding', 'short': 'CsrA_binding'},
        'RM00006': {'name': 'CUYG tetraloop', 'short': 'CUYG'},
        'RM00007': {'name': 'Splicing Domain V', 'short': 'Domain-V'},
        'RM00008': {'name': 'GNRA tetraloop', 'short': 'GNRA'},
        'RM00009': {'name': 'HuR binding', 'short': 'HuR_binding'},
        'RM00010': {'name': 'Kink-turn', 'short': 'K-turn'},
        'RM00011': {'name': 'Kink-turn type 2', 'short': 'K-turn-2'},
        'RM00012': {'name': 'Pseudo-kink-turn', 'short': 'pK-turn'},
        'RM00013': {'name': 'RBS B. subtilis', 'short': 'RBS_B_subtilis'},
        'RM00014': {'name': 'RBS E. coli', 'short': 'RBS_E_coli'},
        'RM00015': {'name': 'RBS H. pylori', 'short': 'RBS_H_pylori'},
        'RM00016': {'name': 'Right angle type 2', 'short': 'right_angle-2'},
        'RM00017': {'name': 'Right angle type 3', 'short': 'right_angle-3'},
        'RM00018': {'name': 'Sarcin-ricin type 1', 'short': 'sarcin-ricin-1'},
        'RM00019': {'name': 'Sarcin-ricin type 2', 'short': 'sarcin-ricin-2'},
        'RM00020': {'name': 'SRP S domain', 'short': 'SRP_S_domain'},
        'RM00021': {'name': 'Tandem GA/AG', 'short': 'tandem-GA'},
        'RM00022': {'name': 'Rho terminator 1', 'short': 'Terminator1'},
        'RM00023': {'name': 'Rho terminator 2', 'short': 'Terminator2'},
        'RM00024': {'name': 'T-loop', 'short': 'T-loop'},
        'RM00025': {'name': 'TRIT', 'short': 'TRIT'},
        'RM00026': {'name': 'Twist-up', 'short': 'twist_up'},
        'RM00027': {'name': 'UAA/GAN', 'short': 'UAA_GAN'},
        'RM00028': {'name': 'UMAC tetraloop', 'short': 'UMAC'},
        'RM00029': {'name': 'UNCG tetraloop', 'short': 'UNCG'},
        'RM00030': {'name': 'U-turn', 'short': 'U-turn'},
        'RM00031': {'name': 'VapC target', 'short': 'vapC_target'},
        'RM00032': {'name': 'Docking elbow', 'short': 'docking_elbow'},
        'RM00033': {'name': 'VTS1 binding', 'short': 'VTS1_binding'},
        'RM00034': {'name': 'Roquin binding', 'short': 'Roquin_binding'},
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
        self._converter = StockholmConverter()
        self._fetched_pdbs: Set[str] = set()
        # motif RM ID -> {PDB_ID: [MotifInstance, ...]}
        self._motif_instances_cache: Dict[str, Dict[str, List[MotifInstance]]] = {}
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
        
        Downloads the motif SEED alignment in Stockholm format and parses
        it with `StockholmConverter`, then filters for the requested PDB.
        
        Args:
            pdb_id: PDB ID
            rfam_motif_id: Rfam motif ID (e.g., 'RM00008')
            motif_info: Motif metadata dict
            
        Returns:
            List of MotifInstance objects
        """
        pdb_instances = self._get_parsed_instances_for_motif(rfam_motif_id, motif_info)
        return pdb_instances.get(pdb_id, [])
    
    def _get_parsed_instances_for_motif(
        self, rfam_motif_id: str, motif_info: Dict
    ) -> Dict[str, List[MotifInstance]]:
        """
        Download and parse the SEED alignment for a single Rfam motif.
        
        Returns a dict mapping PDB IDs to their MotifInstance lists.
        Results are cached per motif ID.
        """
        if rfam_motif_id in self._motif_instances_cache:
            return self._motif_instances_cache[rfam_motif_id]
        
        url = (
            f"{self.API_BASE_URL}/motif/{rfam_motif_id}/alignment"
            f"?acc={rfam_motif_id}&alnType=seed&nseLabels=0&format=stockholm"
        )
        
        try:
            request = urllib.request.Request(
                url,
                headers={
                    'User-Agent': 'RSMViewer/2.0',
                    'Accept': 'text/plain',
                }
            )
            
            ssl_context = ssl.create_default_context()
            ssl_context.check_hostname = False
            ssl_context.verify_mode = ssl.CERT_NONE
            
            with urllib.request.urlopen(request, timeout=self.REQUEST_TIMEOUT, context=ssl_context) as response:
                if response.status == 200:
                    body = response.read().decode('utf-8', errors='replace')
                    if not body or not body.startswith('# STOCKHOLM'):
                        self._motif_instances_cache[rfam_motif_id] = {}
                        return {}
                    
                    motif_types = self._converter.convert_data(body, {
                        'type_id': motif_info['short'],
                        'name': motif_info['name'],
                        'file': url,
                        'source': 'Rfam API',
                    })
                    
                    # Group instances by PDB ID
                    pdb_instances: Dict[str, List[MotifInstance]] = {}
                    for mt in motif_types:
                        for inst in mt.instances:
                            pdb_instances.setdefault(inst.pdb_id, []).append(inst)
                    
                    self._motif_instances_cache[rfam_motif_id] = pdb_instances
                    return pdb_instances
                    
        except (urllib.error.HTTPError, urllib.error.URLError):
            pass
        except Exception:
            pass
        
        self._motif_instances_cache[rfam_motif_id] = {}
        return {}

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
