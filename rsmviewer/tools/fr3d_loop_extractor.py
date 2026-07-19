#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fr3d_loop_extractor.py

Extracts RNA structural loops (Hairpin Loops HL, Internal Loops IL) from a
CIF/PDB file by running the FR3D pairwise interaction annotator and applying
the flankSS loop-boundary algorithm described in:

    Petrov et al. "Automated classification of RNA 3D motifs and the RNA 3D
    Motif Atlas." RNA 19:1327-1340 (2013). doi:10.1261/rna.039438.113

Pipeline:
  1. Run NA_pairwise_interactions.py  (FR3D annotation)
  2. Parse output: collect canonical cWW basepairs
  3. Apply flankSS logic to identify loop boundaries
  4. Classify loops: HL (hairpin), IL (internal loop)
  5. Output {PDB_ID}_fr3d_loops.csv  -- same format as
     https://rna.bgsu.edu/rna3dhub/loops/download/{PDB_ID}

Usage (standalone):
    python fr3d_loop_extractor.py --config fr3d_pipeline_config.json --pdb 1S72
    python fr3d_loop_extractor.py --config fr3d_pipeline_config.json --pdb 1S72 --cif /path/1s72.cif

Called by RSMViewer gui.py when:
    rmv_db 5 /path/to/fr3d_pipeline_config.json
    rmv_fetch 1S72          (or rmv_fr3d run 1S72)
"""

import argparse
import json
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# Canonical Watson-Crick and G-U wobble base pairs that form helical stems.
# (Petrov et al. 2013: "flanking base pairs" == canonical cWW pairs)
CANONICAL_PAIRS = frozenset([
    ('A', 'U'), ('U', 'A'),
    ('G', 'C'), ('C', 'G'),
    ('G', 'U'), ('U', 'G'),
])

# FR3D interaction labels that count as canonical cWW.
# Includes true (cWW), alternate geometry (acWW), and mixed-case variants.
def _is_canonical_cww(interaction: str) -> bool:
    """Return True iff interaction is a true (non-near) canonical cWW pair."""
    s = interaction.strip()
    # Reject "near" pairs (prefix 'n')
    if s.startswith('n') or s.startswith('N'):
        return False
    # Accept cWW, acWW, cWw, cwW, acWw, acwW  (case-insensitive after stripping 'a')
    core = s.lstrip('a').lstrip('A')
    return core.lower().startswith('cww')


# ─────────────────────────────────────────────────────────────────────────────
# Unit ID helpers
# ─────────────────────────────────────────────────────────────────────────────

def _parse_uid(uid: str):
    """Parse BGSU unit ID: pdbid|model|chain|base|seqnum[|ins|alt]
    Returns (pdb_id, model, chain, base, seqnum:int) or None."""
    parts = uid.split('|')
    if len(parts) < 5:
        return None
    try:
        return parts[0].upper(), parts[1], parts[2], parts[3].upper(), int(parts[4])
    except (ValueError, IndexError):
        return None


def _chain_key(uid: str):
    """Return (pdb_id, model, chain) tuple from a unit ID, or None."""
    p = _parse_uid(uid)
    if p is None:
        return None
    return (p[0], p[1], p[2])


def _seqnum(uid: str) -> int:
    """Return sequence number from unit ID, or -1 on failure."""
    p = _parse_uid(uid)
    return p[4] if p is not None else -1


def _base(uid: str) -> str:
    """Return single-letter base from unit ID."""
    p = _parse_uid(uid)
    return p[3] if p is not None else ''


# ─────────────────────────────────────────────────────────────────────────────
# Parsing FR3D output
# ─────────────────────────────────────────────────────────────────────────────

def load_basepairs(basepairs_file: str):
    """Parse FR3D _basepairs.txt output (tab-separated).

    File format (from FR3D README):
        uid1 \\t interaction \\t uid2 \\t crossing_number

    Returns list of (uid1, interaction, uid2, crossing:int).
    """
    pairs = []
    with open(basepairs_file, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            uid1 = parts[0].strip()
            interaction = parts[1].strip()
            uid2 = parts[2].strip()
            try:
                crossing = int(parts[3]) if len(parts) > 3 else 0
            except ValueError:
                crossing = 0
            pairs.append((uid1, interaction, uid2, crossing))
    return pairs


def extract_canonical_cww_pairs(pairs):
    """Filter pairs to retain only TRUE NESTED canonical cWW pairs.

    Two criteria from Petrov et al. 2013 (RNA 3D Motif Atlas paper):
    1. True (non-near) cWW interaction — 'near' pairs (prefix 'n') are excluded.
       The paper states: "The closing pairs are required to make cWW pairs;
       'near' base pairs defined by FR3D are not allowed."
    2. Crossing number = 0 (nested, not pseudoknot). FR3D's 4th column records
       how many nested pairs the interaction crosses.  Only nested pairs
       (crossing=0) form the regular secondary structure helix framework used
       to identify loop boundaries.  Pseudoknot pairs (crossing>0) cross the
       secondary structure and must NOT be used as loop boundaries.

    Returns list of (uid1, uid2).
    """
    canonical = []
    for uid1, interaction, uid2, crossing in pairs:
        if not _is_canonical_cww(interaction):
            continue
        # Exclude pseudoknot base pairs — these cross the secondary structure
        # and do not constitute regular helix boundaries.
        if crossing != 0:
            continue
        b1 = _base(uid1)
        b2 = _base(uid2)
        if (b1, b2) not in CANONICAL_PAIRS:
            continue
        canonical.append((uid1, uid2))
    return canonical


# ─────────────────────────────────────────────────────────────────────────────
# Secondary structure graph
# ─────────────────────────────────────────────────────────────────────────────

def build_secondary_structure(canonical_pairs, all_uids):
    """Build the secondary structure from canonical cWW pairs.

    Returns:
        paired   : dict {uid -> partner_uid}
        chain_map: dict {(pdb,model,chain) -> sorted list of (seqnum, uid)}
    """
    paired = {}
    for uid1, uid2 in canonical_pairs:
        paired[uid1] = uid2
        paired[uid2] = uid1

    # Group nucleotides by chain
    chain_map = defaultdict(list)
    for uid in all_uids:
        ck = _chain_key(uid)
        if ck is None:
            continue
        sn = _seqnum(uid)
        if sn < 0:
            continue
        chain_map[ck].append((sn, uid))

    # Sort each chain by sequence position
    for ck in chain_map:
        chain_map[ck].sort()

    return paired, dict(chain_map)


# ─────────────────────────────────────────────────────────────────────────────
# flankSS helpers
# ─────────────────────────────────────────────────────────────────────────────

def _has_nested_cww(chain_seqnum_list, paired, lo, hi, chain_key):
    """Return True iff any canonical cWW pair is nested strictly between
    positions lo and hi (exclusive) on chain_key.

    "Nested" means: a pair (k, l) where lo < k < l < hi and both k, l
    are on chain_key.
    """
    uid_at = {sn: uid for sn, uid in chain_seqnum_list}
    for sn, uid in chain_seqnum_list:
        if sn <= lo or sn >= hi:
            continue
        partner = paired.get(uid)
        if partner is None:
            continue
        ck_partner = _chain_key(partner)
        sn_partner = _seqnum(partner)
        if ck_partner == chain_key and lo < sn_partner < hi:
            return True
    return False


# ─────────────────────────────────────────────────────────────────────────────
# Hairpin loop extraction
# ─────────────────────────────────────────────────────────────────────────────

def extract_hairpin_loops(chain_map, paired):
    """Extract hairpin loops (HL) using the flankSS definition.

    Hairpin: two nucleotides X (pos i) and Y (pos j, i < j) where
        - X is paired with Y  (same chain)
        - No canonical cWW pair is nested between positions i and j
        - At least one nucleotide exists between i and j

    Closing pair (X, Y) is included in the loop residue list, as in the
    Motif Atlas convention.

    Returns list of ('HL', [uid_x, ...internal..., uid_y]).
    """
    loops = []
    seen = set()

    for chain_key, nts in chain_map.items():
        uid_at = {sn: uid for sn, uid in nts}
        positions = sorted(uid_at)

        for sn_x, uid_x in nts:
            partner = paired.get(uid_x)
            if partner is None:
                continue
            ck_p = _chain_key(partner)
            sn_p = _seqnum(partner)
            if ck_p != chain_key or sn_p <= sn_x:
                continue  # not same chain, or partner comes before X

            between = [s for s in positions if sn_x < s < sn_p]
            if not between:
                continue  # no residues between → degenerate, skip

            if _has_nested_cww(nts, paired, sn_x, sn_p, chain_key):
                continue  # nested pair → this is a helix, not a hairpin

            key = (chain_key, sn_x, sn_p)
            if key in seen:
                continue
            seen.add(key)

            loop_uids = [uid_x] + [uid_at[s] for s in between] + [partner]
            loops.append(('HL', loop_uids))

    return loops


# ─────────────────────────────────────────────────────────────────────────────
# Internal loop extraction
# ─────────────────────────────────────────────────────────────────────────────

def extract_internal_loops(chain_map, paired):
    """Extract internal loops (IL) using the flankSS definition.

    Internal loop: two consecutive canonical cWW pairs P1=(A,B) and P2=(C,D)
    forming the closing base pairs of an internal loop where:
        - A and C are on strand 1 (same chain, A < C)
        - B and D are on strand 2 (same or different chain, D < B)
        - No canonical cWW pair is nested between A and C on strand 1
        - No canonical cWW pair is nested between D and B on strand 2
        - At least one unpaired residue on at least one strand

    Includes closing pairs (A, C, D, B) in the loop residue list.

    Returns list of ('IL', [uid_A, ...strand1..., uid_C, uid_D, ...strand2..., uid_B]).
    """
    loops = []
    seen = set()

    # For each chain, collect the 5'-end nucleotide of each canonical cWW pair
    # i.e. the member of the pair that has the smaller sequence number on its chain.
    cww_on_chain = defaultdict(list)  # chain_key -> [(seqnum, uid_5prime, partner_uid)]
    for chain_key, nts in chain_map.items():
        for sn, uid in nts:
            partner = paired.get(uid)
            if partner is None:
                continue
            ck_p = _chain_key(partner)
            sn_p = _seqnum(partner)
            # Record the pair once: the member with smaller seqnum on its chain
            # For cross-chain pairs: uid is on chain_key, partner is on another chain
            # For same-chain pairs: only record when sn < sn_p
            if ck_p == chain_key and sn_p < sn:
                continue  # will be recorded from the other side
            cww_on_chain[chain_key].append((sn, uid, partner))

    for chain_key, cw_list in cww_on_chain.items():
        cw_list_sorted = sorted(cw_list)  # by seqnum on strand 1
        uid_at_1 = {sn: uid for sn, uid in chain_map.get(chain_key, [])}

        for i in range(len(cw_list_sorted) - 1):
            sn_a, uid_a, uid_b = cw_list_sorted[i]
            sn_c, uid_c, uid_d = cw_list_sorted[i + 1]

            # Strand 2: B and D must be on the same chain, with D < B
            ck_b = _chain_key(uid_b)
            ck_d = _chain_key(uid_d)
            if ck_b != ck_d:
                continue
            sn_b = _seqnum(uid_b)
            sn_d = _seqnum(uid_d)
            if sn_d >= sn_b:
                continue  # D must come before B on strand 2

            # Strand 1: no canonical cWW pair nested between A and C
            nts_s1 = chain_map.get(chain_key, [])
            if _has_nested_cww(nts_s1, paired, sn_a, sn_c, chain_key):
                continue

            # Strand 2: no canonical cWW pair nested between D and B
            nts_s2 = chain_map.get(ck_b, [])
            if _has_nested_cww(nts_s2, paired, sn_d, sn_b, ck_b):
                continue

            # Must have at least one unpaired residue on at least one strand
            between_ac = [s for s in sorted(uid_at_1) if sn_a < s < sn_c]
            uid_at_2 = {sn: uid for sn, uid in nts_s2}
            between_db = [s for s in sorted(uid_at_2) if sn_d < s < sn_b]

            if not between_ac and not between_db:
                continue  # adjacent closing pairs with no loop residues → skip

            key = frozenset([uid_a, uid_b, uid_c, uid_d])
            if key in seen:
                continue
            seen.add(key)

            # Build ordered residue list:
            #   strand 1: A, residues(A→C), C
            #   strand 2: D, residues(D→B), B
            loop_uids = (
                [uid_a]
                + [uid_at_1[s] for s in between_ac]
                + [uid_c, uid_d]
                + [uid_at_2[s] for s in between_db]
                + [uid_b]
            )
            loops.append(('IL', loop_uids))

    return loops


def extract_junction_loops(chain_map, paired, max_order: int = 9):
    """Extract junction loops (J3–J9) for same-chain RNA.

    Algorithm (from BGSU FR3DMotifs/extractLoops.m, aFlankJunctionSearch):
    BGSU searches for a CYCLIC pattern of n closing pairs (J_n junction
    has n arms = n closing pairs).  For a single-chain linear RNA, the
    canonical topology is:

        Outer pair P0=(A,B): the closing pair of the helix arm leading
        INTO the junction (A < B).
        n-1 immediate inner pairs P1=(C1,D1)…Pn-1=(Cn-1,Dn-1) at the
        same nesting level inside P0 (A < C1 < D1 < … < Cn-1 < Dn-1 < B),
        where none of the inner pairs is nested inside another.

        FlankSS conditions between consecutive endpoints:
          flankSS(A, C1)        — junction strand 0
          flankSS(Di, Ci+1)    — junction strand i (i=1..n-2)
          flankSS(Dn-1, B)     — junction strand n-1

        This yields a J_n junction (n = number of immediate inner pairs + 1).

    Limitation: multi-chain junctions (arms on different chains) are not
    handled here; use `rmv_fr3d run <PDB>` (BGSU download) for complete
    junction coverage including multi-chain cases.

    Returns list of (junction_type, [uid_list]) e.g. ('J3', [...]).
    """
    loops = []
    seen: set = set()

    for chain_key, nts in chain_map.items():
        uid_at = {sn: uid for sn, uid in nts}
        positions = sorted(uid_at)

        # Collect same-chain canonical cWW pairs, 5' member first
        same_chain_pairs = []
        for sn_a, uid_a in nts:
            partner = paired.get(uid_a)
            if partner is None:
                continue
            ck_p = _chain_key(partner)
            sn_p = _seqnum(partner)
            if ck_p == chain_key and sn_p > sn_a:
                same_chain_pairs.append((sn_a, sn_p, uid_a, partner))
        same_chain_pairs.sort()

        for sn_a, sn_b, uid_a, uid_b in same_chain_pairs:

            # All canonical cWW pairs strictly inside (sn_a, sn_b)
            inner = [
                (sn_c, sn_d, uid_c, uid_d)
                for sn_c, sn_d, uid_c, uid_d in same_chain_pairs
                if sn_a < sn_c and sn_d < sn_b
            ]

            # Keep only IMMEDIATE inner pairs (not nested inside any other inner pair)
            immediate = []
            for pair in sorted(inner):
                sn_c, sn_d, uid_c, uid_d = pair
                if not any(sn_x < sn_c and sn_d < sn_y
                           for sn_x, sn_y, _, _ in immediate):
                    immediate.append(pair)

            n = len(immediate)
            if n < 2 or n + 1 > max_order:
                continue  # J2 = IL (handled elsewhere); above max_order: skip

            # Verify flankSS between consecutive pair endpoints
            # (no canonical cWW nested between consecutive boundary positions)
            if _has_nested_cww(nts, paired, sn_a, immediate[0][0], chain_key):
                continue
            valid = True
            for k in range(n - 1):
                if _has_nested_cww(nts, paired,
                                   immediate[k][1], immediate[k + 1][0],
                                   chain_key):
                    valid = False
                    break
            if not valid:
                continue
            if _has_nested_cww(nts, paired, immediate[-1][1], sn_b, chain_key):
                continue

            # Deduplication key: frozenset of all closing-pair unit IDs
            pair_uids = frozenset(
                [uid_a, uid_b]
                + [u for p in immediate for u in (p[2], p[3])]
            )
            if pair_uids in seen:
                continue
            seen.add(pair_uids)

            # Build residue list (BGSU convention: closing pairs + SS strands)
            loop_uids = [uid_a]
            for sn in positions:          # strand between A and first inner 5'
                if sn_a < sn < immediate[0][0]:
                    loop_uids.append(uid_at[sn])
            for k, (sn_c, sn_d, uid_c, uid_d) in enumerate(immediate):
                loop_uids.extend([uid_c, uid_d])
                if k < n - 1:            # strand between Di and next Ci+1
                    for sn in positions:
                        if sn_d < sn < immediate[k + 1][0]:
                            loop_uids.append(uid_at[sn])
            for sn in positions:          # strand between last inner 3' and B
                if immediate[-1][1] < sn < sn_b:
                    loop_uids.append(uid_at[sn])
            loop_uids.append(uid_b)

            loops.append((f"J{n + 1}", loop_uids))

    return loops


# ─────────────────────────────────────────────────────────────────────────────
# CSV output (BGSU format)
# ─────────────────────────────────────────────────────────────────────────────

def write_bgsu_csv(hl_loops, il_loops, output_file: str, pdb_id: str,
                   loop_types=None, min_loop_size=3, max_loop_size=500,
                   chains=None, jn_loops=None) -> int:
    """Write loops in BGSU loops/download CSV format, applying structural filters.

    Filters applied (from Petrov et al. 2013 quality criteria):
        loop_types     : only write listed types (default all)
        min_loop_size  : skip loops with fewer than this many residues
        max_loop_size  : skip loops with more than this many residues
        chains         : if non-empty, only keep loops whose residues are all
                         on one of the listed chain IDs

    Returns total number of loop instances written.
    """
    if loop_types is None:
        loop_types = ['HL', 'IL', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9']
    allowed_types = set(t.upper() for t in loop_types)
    chain_filter = set(chains) if chains else set()

    lines = []
    counters = {'HL': 0, 'IL': 0}

    all_loops = [('HL', hl_loops), ('IL', il_loops)]
    # Add junction loops grouped by type (J3, J4, …)
    if jn_loops:
        from collections import defaultdict
        jn_by_type = defaultdict(list)
        for lt, uids in jn_loops:
            jn_by_type[lt].append((lt, uids))
        for jt in sorted(jn_by_type.keys()):
            all_loops.append((jt, jn_by_type[jt]))
    for loop_type, loop_list in all_loops:
        if loop_type not in allowed_types:
            continue
        idx = 0
        for _lt, uids in loop_list:
            n = len(uids)
            if n < min_loop_size or n > max_loop_size:
                continue
            if chain_filter:
                loop_chains = set()
                for uid in uids:
                    p = _parse_uid(uid)
                    if p:
                        loop_chains.add(p[2])
                if not loop_chains.issubset(chain_filter):
                    continue
            idx += 1
            loop_id = f"{loop_type}_{pdb_id.upper()}_{idx:03d}"
            lines.append(f'"{loop_id}","{",".join(uids)}"')
        counters[loop_type] = idx

    with open(output_file, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines))
        if lines:
            fh.write('\n')

    return len(lines)


# ─────────────────────────────────────────────────────────────────────────────
# FR3D annotation runner
# ─────────────────────────────────────────────────────────────────────────────

def run_fr3d_annotation(fr3d_python_dir: str, python_exe: str,
                         cif_file: str, output_dir: str, pdb_id: str):
    """Run NA_pairwise_interactions.py to produce _basepairs.txt.

    Returns path to the basepairs file, or None on failure.
    """
    # Locate classifier script
    candidates = [
        os.path.join(fr3d_python_dir, 'fr3d', 'classifiers', 'NA_pairwise_interactions.py'),
        os.path.join(fr3d_python_dir, 'NA_pairwise_interactions.py'),
    ]
    script = next((c for c in candidates if os.path.isfile(c)), None)
    if script is None:
        print(f"[fr3d_loop_extractor] ERROR: NA_pairwise_interactions.py not found in {fr3d_python_dir}")
        return None

    cif_dir = os.path.dirname(os.path.abspath(cif_file))
    pdb_upper = pdb_id.upper()

    cmd = [
        python_exe, script,
        '--input', cif_dir,
        '--output', output_dir,
        '-c', 'basepair',
        pdb_upper,
    ]

    print(f"[fr3d_loop_extractor] Running: {' '.join(cmd)}")

    env = os.environ.copy()
    env['PYTHONPATH'] = fr3d_python_dir + os.pathsep + env.get('PYTHONPATH', '')

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600, env=env
        )
        if result.returncode != 0:
            print(f"[fr3d_loop_extractor] FR3D failed (exit {result.returncode})")
            if result.stderr:
                print(result.stderr[-3000:])
            return None
    except subprocess.TimeoutExpired:
        print("[fr3d_loop_extractor] ERROR: FR3D annotation timed out (>600 s)")
        return None
    except Exception as exc:
        print(f"[fr3d_loop_extractor] ERROR: {exc}")
        return None

    # Find output file — FR3D writes _basepair.txt (no trailing 's') on some versions
    for name in [
        f'{pdb_upper}_basepair.txt', f'{pdb_id.lower()}_basepair.txt',
        f'{pdb_upper}_basepairs.txt', f'{pdb_id.lower()}_basepairs.txt',
    ]:
        path = os.path.join(output_dir, name)
        if os.path.isfile(path):
            print(f"[fr3d_loop_extractor] Basepairs file: {path}")
            return path

    print(f"[fr3d_loop_extractor] ERROR: Expected basepairs file not found in {output_dir}")
    return None


# ─────────────────────────────────────────────────────────────────────────────
# CIF download helper
# ─────────────────────────────────────────────────────────────────────────────

def _download_cif(pdb_id: str, dest_dir: str) -> str:
    """Download CIF from RCSB. Returns local path or empty string on failure."""
    import ssl
    import urllib.request
    pdb_lower = pdb_id.lower()
    url = f'https://files.rcsb.org/download/{pdb_lower}.cif.gz'
    dest = os.path.join(dest_dir, f'{pdb_lower}.cif.gz')
    print(f"[fr3d_loop_extractor] Downloading CIF: {url}")
    # Try with verified SSL first; fall back to unverified (common on macOS)
    contexts = [ssl.create_default_context(), ssl._create_unverified_context()]
    for ctx in contexts:
        try:
            with urllib.request.urlopen(url, context=ctx, timeout=60) as r:
                data = r.read()
            with open(dest, 'wb') as fh:
                fh.write(data)
            return dest
        except Exception as exc:
            last_exc = exc
            continue
    print(f"[fr3d_loop_extractor] Download failed: {last_exc}")
    return ''


# ─────────────────────────────────────────────────────────────────────────────
# High-level API (called by RSMViewer gui.py)
# ─────────────────────────────────────────────────────────────────────────────

def _load_all_rna_nucleotides(cif_file: str, pdb_id: str, fr3d_dir: str) -> set:
    """Load ALL RNA nucleotide unit IDs from the CIF file using FR3D's structure loader.

    The Motif Atlas pipeline (Petrov et al. 2013) extracts loops from the full
    RNA chain, including residues with no basepair interactions.  This function
    ensures the chain maps are complete so that loop residue lists are
    continuous (no gaps from unpaired residues).

    Returns a set of BGSU-format unit IDs, or empty set on failure.
    """
    if fr3d_dir and fr3d_dir not in sys.path:
        sys.path.insert(0, fr3d_dir)
    try:
        from fr3d.classifiers.NA_pairwise_interactions import load_structure  # type: ignore
        structure, messages = load_structure(cif_file, pdb_id)
        if structure is None:
            return set()
        all_uids = set()
        for nt in structure.residues(type=["RNA linking", "DNA linking"]):
            uid = nt.unit_id()
            if uid:
                all_uids.add(uid)
        print(f"[fr3d_loop_extractor] Full nucleotide list: {len(all_uids)} residues from structure")
        return all_uids
    except Exception as exc:
        print(f"[fr3d_loop_extractor] Warning: could not load full structure ({exc}); "
              f"loop residue lists may be discontinuous")
        return set()


def run_pipeline(config: dict, pdb_id: str, cif_file: str = '') -> str:
    """Run the full FR3D loop extraction pipeline.

    Args:
        config   : Parsed fr3d_pipeline_config.json dict
        pdb_id   : PDB ID (e.g. '1S72')
        cif_file : Optional explicit path to CIF file

    Returns:
        Path to the output CSV file on success, empty string on failure.
    """
    pdb_id = pdb_id.strip().upper()
    fr3d_dir = os.path.expanduser(str(config.get('fr3d_python_dir', '') or ''))
    python_exe = str(config.get('python_executable', '') or sys.executable).strip() or sys.executable
    output_dir = os.path.expanduser(str(config.get('output_dir', '.') or '.'))
    cif_input_dir = os.path.expanduser(str(config.get('cif_input_dir', '') or ''))
    auto_dl = bool(config.get('auto_download_cif', True))
    force_fresh = bool(config.get('force_fresh', False))

    if not fr3d_dir or not os.path.isdir(fr3d_dir):
        print(f"[fr3d_loop_extractor] ERROR: fr3d_python_dir not found: {fr3d_dir}")
        return ''

    os.makedirs(output_dir, exist_ok=True)

    # Resolve CIF file
    if not cif_file or not os.path.isfile(cif_file):
        pdb_lower = pdb_id.lower()
        search_dirs = [d for d in [cif_input_dir, output_dir, '.'] if d]
        for sd in search_dirs:
            for name in [f'{pdb_lower}.cif', f'{pdb_id}.cif',
                         f'{pdb_lower}.cif.gz', f'{pdb_id}.cif.gz']:
                candidate = os.path.join(sd, name)
                if os.path.isfile(candidate):
                    cif_file = candidate
                    break
            if cif_file:
                break

    if not cif_file or not os.path.isfile(cif_file):
        if auto_dl:
            cif_file = _download_cif(pdb_id, output_dir)
        if not cif_file or not os.path.isfile(cif_file):
            print(f"[fr3d_loop_extractor] ERROR: CIF file not found for {pdb_id}")
            return ''

    print(f"[fr3d_loop_extractor] CIF: {cif_file}")

    # Run FR3D annotation.
    # When force_fresh is enabled, previous basepair outputs are removed so
    # NA_pairwise_interactions.py runs every time.
    bpf = ''
    basepair_names = [
        f'{pdb_id}_basepair.txt', f'{pdb_id.lower()}_basepair.txt',
        f'{pdb_id}_basepairs.txt', f'{pdb_id.lower()}_basepairs.txt',
    ]

    if force_fresh:
        for _name in basepair_names:
            _candidate = os.path.join(output_dir, _name)
            if os.path.isfile(_candidate):
                try:
                    os.remove(_candidate)
                    print(f"[fr3d_loop_extractor] Removed cached basepairs file: {_candidate}")
                except OSError as exc:
                    print(f"[fr3d_loop_extractor] Warning: could not remove {_candidate}: {exc}")
    else:
        for _name in basepair_names:
            _candidate = os.path.join(output_dir, _name)
            if os.path.isfile(_candidate):
                bpf = _candidate
                print(f"[fr3d_loop_extractor] Using existing basepairs file: {bpf}")
                break

    if not bpf:
        bpf = run_fr3d_annotation(fr3d_dir, python_exe, cif_file, output_dir, pdb_id)
        if bpf is None:
            return ''

    # Parse annotations
    print("[fr3d_loop_extractor] Parsing pairwise interactions...")
    all_pairs = load_basepairs(bpf)
    print(f"  {len(all_pairs)} total interactions")

    canonical = extract_canonical_cww_pairs(all_pairs)
    print(f"  {len(canonical)} canonical cWW pairs")

    # Load ALL nucleotides from the CIF structure so chain maps are complete.
    # The Motif Atlas pipeline uses the full chain (all residues, not just
    # those in basepairs) so that loop residue lists are continuous with no
    # gaps from unpaired residues.
    all_uids = _load_all_rna_nucleotides(cif_file, pdb_id, fr3d_dir)
    if not all_uids:
        # Fallback: only nucleotides seen in basepairs (produces gapped loops)
        print("[fr3d_loop_extractor] Fallback: using nucleotides from basepairs only "
              "(loops may have discontinuous residue ranges)")
        for u1, _i, u2, _c in all_pairs:
            all_uids.add(u1)
            all_uids.add(u2)
    print(f"  {len(all_uids)} unique nucleotides")

    # Build secondary structure
    paired, chain_map = build_secondary_structure(canonical, all_uids)
    print(f"  {len(chain_map)} RNA chains")

    # Extract loops
    print("[fr3d_loop_extractor] Extracting loops...")
    hl_loops = extract_hairpin_loops(chain_map, paired)
    il_loops = extract_internal_loops(chain_map, paired)
    requested_types = config.get('loop_types') or ['HL', 'IL']
    jn_types = [t for t in requested_types if t.startswith('J')]
    jn_loops = []
    if jn_types:
        max_jn = max((int(t[1:]) for t in jn_types if t[1:].isdigit()), default=9)
        jn_loops = extract_junction_loops(chain_map, paired, max_order=max_jn)
    print(f"  Hairpin loops  (HL): {len(hl_loops)}")
    print(f"  Internal loops (IL): {len(il_loops)}")
    if jn_loops:
        from collections import Counter
        jn_counts = Counter(lt for lt, _ in jn_loops)
        for jt in sorted(jn_counts):
            print(f"  {jt} junction loops   : {jn_counts[jt]}")

    # Write CSV with structural filters from config
    out_csv = os.path.join(output_dir, f'{pdb_id}_fr3d_loops.csv')
    total = write_bgsu_csv(
        hl_loops, il_loops, out_csv, pdb_id,
        loop_types=requested_types,
        min_loop_size=int(config.get('min_loop_size', 3)),
        max_loop_size=int(config.get('max_loop_size', 500)),
        chains=config.get('chains') or [],
        jn_loops=jn_loops,
    )
    print(f"[fr3d_loop_extractor] Wrote {total} loop instances → {out_csv}")

    return out_csv if total > 0 else ''


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Extract HL/IL loops from CIF using FR3D (Motif Atlas pipeline)'
    )
    parser.add_argument('--config', required=True,
                        help='Path to fr3d_pipeline_config.json')
    parser.add_argument('--pdb', required=True,
                        help='PDB ID (e.g. 1S72)')
    parser.add_argument('--cif', default='',
                        help='Explicit path to CIF file (optional)')
    args = parser.parse_args()

    with open(args.config, 'r', encoding='utf-8') as fh:
        config = json.load(fh)

    result = run_pipeline(config, args.pdb, args.cif)
    if not result:
        print("[fr3d_loop_extractor] Pipeline failed.")
        return 1
    print(f"[fr3d_loop_extractor] Done: {result}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
