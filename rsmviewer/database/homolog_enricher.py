"""
Homolog Enricher - Replace generic motif names with specific annotations

Uses the BGSU Non-Redundant (NR) representative set to look up the
representative structure for a query PDB, then fetches the representative's
semantic annotations from the BGSU HTML scraper. Generic motif names
(HL, IL, J3-J7, etc.) are replaced with specific names (Kink-turn, GNRA,
C-loop, etc.) when the same motif GROUP (structural class) has a specific
annotation in the representative structure.

Pipeline:
  1. Query PDB (e.g., 1S72) -> NR list -> Representative (e.g., 4V9F)
  2. Fetch 4V9F annotations via BGSU hybrid provider
  3. Build motif_group -> annotation lookup from representative
  4. For each generic motif in 1S72:
     - Check its motif_group (e.g., HL_75660.8)
     - Look up that group in the representative's lookup
     - If found with specific name: replace generic with specific
  5. Return enriched motif dataset (original residues, new names)

Author: CBB Lab
Version: 2.0.0
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Set, Tuple

from .base_provider import MotifInstance, ResidueSpec

logger = logging.getLogger(__name__)


# Names that are considered "generic" and candidates for enrichment.
# These are the fallback names assigned when no semantic annotation exists.
GENERIC_NAMES: Set[str] = {
    # Short codes
    'HL', 'IL',
    'J3', 'J4', 'J5', 'J6', 'J7', 'J8',
    # Descriptive names from BGSU provider
    'Hairpin Loop (HL)',
    'Internal Loop (IL)',
    '3-way Junction (J3)',
    '4-way Junction (J4)',
    '5-way Junction (J5)',
    '6-way Junction (J6)',
    '7-way Junction (J7)',
    '8-way Junction (J8)',
    # Other generic names
    'Hairpin Loop', 'Internal Loop', 'Hairpin', 'Internal',
    'Junction', 'Unknown', 'Bulge',
}


def _is_generic_name(name: str) -> bool:
    """Check if a motif name is generic (candidate for enrichment)."""
    return name.strip() in GENERIC_NAMES


def _jaccard_similarity(set_a: Set, set_b: Set) -> float:
    """
    Calculate Jaccard similarity between two sets.
    
    J(A, B) = |A intersection B| / |A union B|
    
    Args:
        set_a: First set
        set_b: Second set
        
    Returns:
        Jaccard similarity (0.0 - 1.0), or 0.0 if both empty
    """
    if not set_a or not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


class HomologEnricher:
    """
    Enriches generic motif names using representative structure annotations.
    
    Takes motifs with generic names (HL, IL, etc.) and tries to find
    the corresponding specific annotation (Kink-turn, GNRA, etc.) from
    the representative structure in the NR equivalence class.
    
    Primary matching strategy: BGSU motif_group IDs.
    Each loop in the BGSU database is assigned to a motif group (structural
    class, e.g., HL_34789.4). If the same motif group has a specific
    annotation in the representative PDB, that annotation is used.
    
    Fallback: Jaccard residue similarity (only for self-representative PDBs).
    """
    
    def __init__(
        self,
        representative_loader,
        bgsu_provider,
        jaccard_threshold: float = 0.60,
    ):
        """
        Initialize the HomologEnricher.
        
        Args:
            representative_loader: RepresentativeSetLoader instance
            bgsu_provider: BGSUAPIProvider instance (for fetching rep annotations)
            jaccard_threshold: Minimum Jaccard similarity for residue-based fallback
        """
        self.rep_loader = representative_loader
        self.bgsu_provider = bgsu_provider
        self.threshold = jaccard_threshold
        # Cache: {rep_pdb: {motif_category: [MotifInstance, ...]}}
        self._rep_annotation_cache: Dict[str, Dict[str, List[MotifInstance]]] = {}
        # Cache: {rep_pdb: {motif_group_id: specific_category_name}}
        self._rep_group_lookup: Dict[str, Dict[str, str]] = {}
    
    def enrich(
        self,
        pdb_id: str,
        motifs: Dict[str, List[MotifInstance]],
    ) -> Dict[str, List[MotifInstance]]:
        """
        Enrich generic motif names using representative homolog annotations.
        
        For each motif with a generic name, looks up its motif_group in
        the representative's group->annotation mapping. If the same group
        has a specific name in the representative, that name is used.
        
        Args:
            pdb_id: The query PDB ID (e.g., '1S72')
            motifs: Dict mapping motif_type -> [MotifInstance, ...]
        
        Returns:
            New motifs dict with generic names replaced where possible.
            Original residues are preserved; only names change.
        """
        pdb_id = pdb_id.upper()
        
        # Get all chain-to-representative mappings for this PDB
        chain_reps = self.rep_loader.get_all_representatives(pdb_id)
        
        if not chain_reps:
            logger.info(f"[ENRICHER] {pdb_id} not in NR list, skipping enrichment")
            return motifs
        
        logger.info(f"[ENRICHER] {pdb_id} has {len(chain_reps)} chain mappings")
        
        # Collect unique representative PDBs we need to fetch
        rep_pdbs_needed = {rep_pdb for (rep_pdb, _) in chain_reps.values()}
        
        # Fetch and build group lookups for each representative PDB
        for rep_pdb in rep_pdbs_needed:
            self._build_group_lookup(rep_pdb)
        
        # Check if any representative has useful lookups
        has_lookup = any(
            rep_pdb in self._rep_group_lookup and self._rep_group_lookup[rep_pdb]
            for rep_pdb in rep_pdbs_needed
        )
        if not has_lookup:
            logger.info(f"[ENRICHER] No representative group lookups available")
            return motifs
        
        # Determine if this is a self-representative case (same PDB)
        is_self_rep = all(
            rep_pdb == pdb_id for rep_pdb, _ in chain_reps.values()
        )
        
        # Now enrich: scan each motif, check if generic, try to match
        enriched: Dict[str, List[MotifInstance]] = {}
        stats = {'total': 0, 'generic': 0, 'enriched': 0, 'kept_generic': 0}
        
        for motif_type, instances in motifs.items():
            for instance in instances:
                stats['total'] += 1
                
                # Check if this motif type is generic
                if not _is_generic_name(motif_type):
                    enriched.setdefault(motif_type, []).append(instance)
                    continue
                
                stats['generic'] += 1
                
                # Try to find a specific name via motif_group
                new_name = self._find_homolog_name(
                    instance, pdb_id, chain_reps, is_self_rep
                )
                
                if new_name and new_name != motif_type:
                    stats['enriched'] += 1
                    enriched_instance = MotifInstance(
                        instance_id=instance.instance_id,
                        motif_id=new_name,
                        pdb_id=instance.pdb_id,
                        residues=instance.residues,
                        annotation=new_name,
                        metadata={
                            **instance.metadata,
                            'enriched_from': motif_type,
                            'enrichment_source': 'homolog',
                        },
                    )
                    enriched.setdefault(new_name, []).append(enriched_instance)
                else:
                    stats['kept_generic'] += 1
                    enriched.setdefault(motif_type, []).append(instance)
        
        logger.info(
            f"[ENRICHER] Results for {pdb_id}: "
            f"{stats['total']} total, {stats['generic']} generic, "
            f"{stats['enriched']} enriched, {stats['kept_generic']} kept generic"
        )
        
        return enriched
    
    def _build_group_lookup(self, rep_pdb: str) -> None:
        """
        Build the motif_group -> specific_category lookup for a representative PDB.
        
        Scans all motifs from the representative and maps each motif_group ID
        to its category name, but ONLY for non-generic categories.
        
        Args:
            rep_pdb: Representative PDB ID
        """
        rep_pdb = rep_pdb.upper()
        
        if rep_pdb in self._rep_group_lookup:
            return
        
        # Fetch representative motifs
        annotations = self._get_rep_annotations(rep_pdb)
        if not annotations:
            self._rep_group_lookup[rep_pdb] = {}
            return
        
        # Build group -> category mapping (non-generic only)
        lookup: Dict[str, str] = {}
        for category, instances in annotations.items():
            if _is_generic_name(category):
                continue  # Skip generic categories
            for inst in instances:
                mg = inst.metadata.get('motif_group', '')
                if mg:
                    lookup[mg] = category
        
        self._rep_group_lookup[rep_pdb] = lookup
        logger.info(
            f"[ENRICHER] Built group lookup for {rep_pdb}: "
            f"{len(lookup)} specific motif groups"
        )
    
    def _get_rep_annotations(
        self, rep_pdb: str
    ) -> Optional[Dict[str, List[MotifInstance]]]:
        """
        Get annotations for a representative PDB, using cache if available.
        
        Always uses the BGSU HTML scraper (even if the original source is
        different) because only BGSU provides semantic annotations.
        
        Args:
            rep_pdb: Representative PDB ID
            
        Returns:
            Dict mapping motif_type -> [MotifInstance, ...], or None
        """
        rep_pdb = rep_pdb.upper()
        
        # Check cache
        if rep_pdb in self._rep_annotation_cache:
            return self._rep_annotation_cache[rep_pdb]
        
        # Fetch from BGSU provider
        try:
            annotations = self.bgsu_provider.get_motifs_for_pdb(rep_pdb)
            if annotations:
                self._rep_annotation_cache[rep_pdb] = annotations
                return annotations
        except Exception as e:
            logger.error(f"[ENRICHER] Error fetching representative {rep_pdb}: {e}")
        
        return None
    
    def _find_homolog_name(
        self,
        instance: MotifInstance,
        pdb_id: str,
        chain_reps: Dict[str, Tuple[str, str]],
        is_self_rep: bool,
    ) -> Optional[str]:
        """
        Find a specific motif name from the representative for a generic motif.
        
        Primary strategy: Use motif_group ID to look up the specific annotation
        in the representative's group->category mapping.
        
        Fallback (self-representative only): Jaccard residue number comparison.
        
        Args:
            instance: The generic motif instance to enrich
            pdb_id: The query PDB ID
            chain_reps: Chain-to-representative mappings
            is_self_rep: Whether the PDB is its own representative
            
        Returns:
            Specific motif name, or None if no match found
        """
        # --- Strategy 1: motif_group matching ---
        motif_group = instance.metadata.get('motif_group', '')
        if motif_group:
            # Check each representative PDB's lookup
            for (rep_pdb, _) in chain_reps.values():
                lookup = self._rep_group_lookup.get(rep_pdb, {})
                if motif_group in lookup:
                    specific_name = lookup[motif_group]
                    logger.debug(
                        f"[ENRICHER] Group match: {instance.instance_id} "
                        f"group={motif_group} -> '{specific_name}'"
                    )
                    return specific_name
        
        # --- Strategy 2: Jaccard fallback (self-representative only) ---
        if is_self_rep:
            return self._jaccard_fallback(instance, pdb_id, chain_reps)
        
        return None
    
    def _jaccard_fallback(
        self,
        instance: MotifInstance,
        pdb_id: str,
        chain_reps: Dict[str, Tuple[str, str]],
    ) -> Optional[str]:
        """
        Fallback matching using Jaccard residue similarity.
        
        Only works for self-representative PDBs where residue numbers
        are identical between member and representative.
        
        Args:
            instance: The generic motif instance
            pdb_id: The PDB ID
            chain_reps: Chain-to-representative mappings
            
        Returns:
            Specific motif name, or None
        """
        motif_chains = instance.get_chains()
        best_name = None
        best_score = 0.0
        
        for mem_chain in motif_chains:
            if mem_chain not in chain_reps:
                continue
            rep_pdb, rep_chain = chain_reps[mem_chain]
            
            annotations = self._rep_annotation_cache.get(rep_pdb)
            if not annotations:
                continue
            
            member_residues = {
                r.residue_number for r in instance.residues
                if r.chain == mem_chain
            }
            if not member_residues:
                continue
            
            for rep_type, rep_instances in annotations.items():
                if _is_generic_name(rep_type):
                    continue
                for rep_inst in rep_instances:
                    rep_residues = {
                        r.residue_number for r in rep_inst.residues
                        if r.chain == rep_chain
                    }
                    if not rep_residues:
                        continue
                    score = _jaccard_similarity(member_residues, rep_residues)
                    if score >= self.threshold and score > best_score:
                        best_score = score
                        best_name = rep_type
        
        return best_name
