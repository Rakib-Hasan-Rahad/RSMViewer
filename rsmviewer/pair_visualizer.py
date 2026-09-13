"""
Base-Pair Visualizer Module for RSMViewer.

Provides the rmv_pair command for visualizing RNA base pairs with
Leontis-Westhof edge labels (Watson-Crick, Hoogsteen, Sugar).

Input format: pdbid_chain1_resnum1_chain2_resnum2
Example:      3j6y_2S_681_2S_696

The command will:
1. Fetch the PDB structure (if not already loaded)
2. Create a pair object with the two selected residues
3. Auto-detect nucleotide identities from the structure
4. Show sticks for the pair, color each residue distinctly
5. Place W, H, S edge labels at the centroid of defining atoms
6. Zoom to the pair
"""

from pymol import cmd


# ------------------------------------------------------------------ #
# Leontis-Westhof edge atom definitions per nucleotide
# ------------------------------------------------------------------ #

EDGE_ATOMS = {
    "A": {
        "W": ["N1", "C2", "N6"],
        "H": ["C5", "C6", "N7", "C8"],
        "S": ["N3", "C2"],
    },
    "G": {
        "W": ["O6", "N1", "N2"],
        "H": ["C5", "N7", "C8"],
        "S": ["N3", "C2", "N2"],
    },
    "C": {
        "W": ["N3", "O2", "N4"],
        "H": ["C5", "C6"],
        "S": ["N1", "C2", "O2"],
    },
    "U": {
        "W": ["O4", "N3", "O2"],
        "H": ["C5", "C6"],
        "S": ["N1", "C2", "O2"],
    },
}

# Default colors for residue 1 and residue 2
COLOR_RES1 = "yellow"
COLOR_RES2 = "cyan"


# Map 3-letter residue names (from PDB) to 1-letter codes for edge lookup
RESNAME_MAP = {
    "A": "A", "G": "G", "C": "C", "U": "U",
    "ADE": "A", "GUA": "G", "CYT": "C", "URA": "U",
    "DA": "A", "DG": "G", "DC": "C", "DT": "U",
    "RA": "A", "RG": "G", "RC": "C", "RU": "U",
}


# ------------------------------------------------------------------ #
# Input parser
# ------------------------------------------------------------------ #

def parse_pair_descriptor(descriptor: str) -> dict:
    """
    Parse a base-pair descriptor string into its components.

    Supported formats (5.field preferred, 7.field legacy):
        pdbid_chain1_resnum1_chain2_resnum2
        pdbid_chain1_resnum1_chain2_resnum2_resname1_resname2

    For the 5-field format, residue names are auto-detected from the
    PDB structure after fetching.

    Returns:
        dict with keys: pdb_id, chain1, resnum1, chain2, resnum2,
                        resname1, resname2  (resnames may be None)
    Raises:
        ValueError on malformed input.
    """
    parts = descriptor.strip().split("_")

    if len(parts) < 5:
        raise ValueError(
            f"Descriptor must have at least 5 underscore-separated fields.\n"
            f"Format: pdbid_chain1_resnum1_chain2_resnum2\n"
            f"Got: {descriptor}"
        )

    # First token is PDB ID
    pdb_id = parts[0].lower()

    # Check if the last two tokens are valid residue names (legacy 7-field)
    resname1 = None
    resname2 = None
    tail = parts[-2:]
    if (len(parts) >= 7
            and tail[0].upper() in EDGE_ATOMS
            and tail[1].upper() in EDGE_ATOMS):
        resname1 = tail[0].upper()
        resname2 = tail[1].upper()
        middle = parts[1:-2]
    else:
        middle = parts[1:]

    if len(middle) != 4:
        raise ValueError(
            f"Expected exactly 4 fields after PDB ID "
            f"(chain1, resnum1, chain2, resnum2).\n"
            f"Got {len(middle)} field(s): {middle}\n"
            f"Format: pdbid_chain1_resnum1_chain2_resnum2"
        )

    chain1 = middle[0]
    resnum1 = middle[1]
    chain2 = middle[2]
    resnum2 = middle[3]

    return {
        "pdb_id": pdb_id,
        "chain1": chain1,
        "resnum1": resnum1,
        "chain2": chain2,
        "resnum2": resnum2,
        "resname1": resname1,
        "resname2": resname2,
    }


# ------------------------------------------------------------------ #
# Edge label placement
# ------------------------------------------------------------------ #

def _place_edge_label(pair_obj: str, chain: str, resi: str,
                      resname: str, label_prefix: str):
    """
    Place W, H, S pseudoatom labels on a single residue.

    Args:
        pair_obj:      Name of the PyMOL pair object.
        chain:         Chain identifier.
        resi:          Residue number (string).
        resname:       One-letter nucleotide name (A, G, C, U).
        label_prefix:  Unique prefix for pseudoatom names.
    """
    edges = EDGE_ATOMS.get(resname)
    if edges is None:
        print(f"  [rmv_pair] Warning: no edge definitions for '{resname}', "
              f"skipping labels for resi {resi}")
        return

    for edge_name, atom_names in edges.items():
        atom_list = "+".join(atom_names)
        sel = (f"{pair_obj} and chain {chain} and resi {resi} "
               f"and name {atom_list}")

        # Check that the selection is non-empty
        if cmd.count_atoms(sel) == 0:
            print(f"  [rmv_pair] Warning: no atoms found for {edge_name} "
                  f"edge of {resname}{resi} (chain {chain})")
            continue

        try:
            xyz = cmd.centerofmass(sel)
        except Exception:
            coords = cmd.get_coords(sel)
            if coords is None or len(coords) == 0:
                continue
            xyz = [sum(c[i] for c in coords) / len(coords) for i in range(3)]

        pseudo_name = f"{label_prefix}_{edge_name}"
        cmd.pseudoatom(pseudo_name, pos=list(xyz))
        cmd.hide("everything", pseudo_name)
        cmd.label(pseudo_name, f'"{edge_name}"')


def _get_resname_from_structure(obj_name: str, chain: str, resi: str) -> str:
    """
    Query PyMOL for the residue name at a given chain/resi and return
    the 1-letter nucleotide code (A, G, C, U) or None.
    """
    stored_names = []
    cmd.iterate(
        f"model {obj_name} and chain {chain} and resi {resi} and name C1'",
        "stored_names.append(resn)",
        space={"stored_names": stored_names},
    )
    if not stored_names:
        # Fallback: try any atom in that residue
        cmd.iterate(
            f"first (model {obj_name} and chain {chain} and resi {resi})",
            "stored_names.append(resn)",
            space={"stored_names": stored_names},
        )
    if stored_names:
        raw = stored_names[0].strip().upper()
        return RESNAME_MAP.get(raw, None)
    return None


# ------------------------------------------------------------------ #
# Main visualization routine
# ------------------------------------------------------------------ #

def visualize_pair(descriptor: str, label_size: int = 24,
                   label_color: str = "black",
                   color1: str = COLOR_RES1,
                   color2: str = COLOR_RES2,
                   zoom_buffer: float = 8.0,
                   bg_style: str = "cartoon"):
    """
    Visualize an RNA base pair from a descriptor string.

    Args:
        descriptor: e.g. "3j6y_2S_681_2S_696"
        label_size: font size for edge labels (default 24)
        label_color: color name for labels (default "black")
        color1: color for residue 1 (default "yellow")
        color2: color for residue 2 (default "cyan")
        zoom_buffer: angstrom buffer around the pair (default 8)
        bg_style: how to show the background structure (default "cartoon",
                  use "none" to hide completely)
    """
    try:
        info = parse_pair_descriptor(descriptor)
    except ValueError as e:
        print(f"\n  [rmv_pair] Error: {e}\n")
        return

    pdb_id = info["pdb_id"]
    chain1 = info["chain1"]
    resnum1 = info["resnum1"]
    chain2 = info["chain2"]
    resnum2 = info["resnum2"]
    resname1 = info["resname1"]
    resname2 = info["resname2"]

    # Step 1: Fetch structure if not already loaded
    existing = cmd.get_object_list()
    existing_lower = [o.lower() for o in existing]
    if pdb_id.lower() in existing_lower:
        obj_idx = existing_lower.index(pdb_id.lower())
        pdb_obj = existing[obj_idx]
        print(f"  [rmv_pair] Structure {pdb_obj} already loaded.")
    else:
        print(f"  [rmv_pair] Fetching {pdb_id}...")
        try:
            cmd.fetch(pdb_id, pdb_id, async_=0)
            pdb_obj = pdb_id
        except Exception as e:
            print(f"  [rmv_pair] Error fetching {pdb_id}: {e}")
            return

    # Step 1b: Auto-detect residue names if not provided
    if resname1 is None:
        resname1 = _get_resname_from_structure(pdb_obj, chain1, resnum1)
        if resname1 is None:
            print(f"  [rmv_pair] Error: could not detect residue name for "
                  f"chain {chain1} resi {resnum1} in {pdb_obj}")
            return
    if resname2 is None:
        resname2 = _get_resname_from_structure(pdb_obj, chain2, resnum2)
        if resname2 is None:
            print(f"  [rmv_pair] Error: could not detect residue name for "
                  f"chain {chain2} resi {resnum2} in {pdb_obj}")
            return

    print(f"\n  [rmv_pair] Visualizing base pair:")
    print(f"    PDB:      {pdb_id}")
    print(f"    Residue1: {resname1} (chain {chain1}, resi {resnum1})")
    print(f"    Residue2: {resname2} (chain {chain2}, resi {resnum2})")

    # Step 2: Build selection strings
    # Use the actual PyMOL object name (pdb_obj) in selections
    if chain1 == chain2:
        resi_combined = f"{resnum1}+{resnum2}"
        create_sel = (f"model {pdb_obj} and chain {chain1} "
                      f"and resi {resi_combined}")
    else:
        create_sel = (f"(model {pdb_obj} and chain {chain1} and resi {resnum1}) or "
                      f"(model {pdb_obj} and chain {chain2} and resi {resnum2})")

    # Sanitize the object name
    safe_chain1 = chain1.replace(" ", "")
    safe_chain2 = chain2.replace(" ", "")
    if chain1 == chain2:
        pair_obj = f"pair_{pdb_id}_{safe_chain1}_{resnum1}_{resnum2}"
    else:
        pair_obj = f"pair_{pdb_id}_{safe_chain1}_{resnum1}_{safe_chain2}_{resnum2}"

    # Step 3: Create pair object
    if cmd.count_atoms(create_sel) == 0:
        print(f"  [rmv_pair] Error: no atoms found for selection:\n"
              f"    {create_sel}")
        print(f"  [rmv_pair] Check chain IDs and residue numbers.")
        return

    cmd.create(pair_obj, create_sel)
    print(f"  [rmv_pair] Created object: {pair_obj}")

    # Step 4: Remove solvent and style the background structure
    cmd.remove("solvent")
    cmd.hide("everything", pdb_obj)
    if bg_style and bg_style.lower() != "none":
        cmd.show(bg_style, pdb_obj)

    # Step 5: Show sticks and color the pair
    cmd.show("sticks", pair_obj)

    color_sel1 = f"{pair_obj} and chain {chain1} and resi {resnum1}"
    color_sel2 = f"{pair_obj} and chain {chain2} and resi {resnum2}"

    cmd.color(color1, color_sel1)
    cmd.color(color2, color_sel2)

    # Step 6: Place edge labels
    label_prefix1 = f"lbl_{pdb_id}_{safe_chain1}_{resnum1}"
    label_prefix2 = f"lbl_{pdb_id}_{safe_chain2}_{resnum2}"

    _place_edge_label(pair_obj, chain1, resnum1, resname1, label_prefix1)
    _place_edge_label(pair_obj, chain2, resnum2, resname2, label_prefix2)

    # Step 7: Label styling
    cmd.set("label_size", label_size)
    cmd.set("label_color", label_color)

    # Step 8: Zoom to the pair
    cmd.zoom(pair_obj, zoom_buffer)

    print(f"  [rmv_pair] Done. Edge labels (W/H/S) placed.\n")


# ------------------------------------------------------------------ #
# Batch visualization
# ------------------------------------------------------------------ #

def visualize_pairs_from_file(filepath: str, **kwargs):
    """
    Read a file with one pair descriptor per line and visualize all.

    Blank lines and lines starting with '#' are skipped.

    Args:
        filepath: path to a text file with descriptors
        **kwargs: passed through to visualize_pair()
    """
    from pathlib import Path

    fpath = Path(filepath)
    if not fpath.exists():
        print(f"\n  [rmv_pair] Error: file not found: {filepath}\n")
        return

    lines = fpath.read_text().strip().splitlines()
    descriptors = [
        line.strip() for line in lines
        if line.strip() and not line.strip().startswith("#")
    ]

    if not descriptors:
        print(f"\n  [rmv_pair] No descriptors found in {filepath}\n")
        return

    print(f"\n  [rmv_pair] Processing {len(descriptors)} pair(s) "
          f"from {filepath}...\n")

    for desc in descriptors:
        visualize_pair(desc, **kwargs)

    print(f"\n  [rmv_pair] Batch complete: {len(descriptors)} pair(s) "
          f"processed.\n")


# ------------------------------------------------------------------ #
# PyMOL command registration
# ------------------------------------------------------------------ #

def register_pair_commands():
    """Register rmv_pair and rmv_pair_batch commands in PyMOL."""

    def _rmv_pair(descriptor="", label_size="24", label_color="black",
                  color1=COLOR_RES1, color2=COLOR_RES2,
                  zoom_buffer="8", bg_style="cartoon"):
        """
        Visualize an RNA base pair with Leontis-Westhof edge labels.

        Usage:
            rmv_pair 3j6y_2S_681_2S_696
            rmv_pair 1s72_A_100_B_200, label_size=30
            rmv_pair 5lzs_5_1582_5_1610, color1=red, color2=blue

        Format: pdbid_chain1_resnum1_chain2_resnum2
        """
        if not descriptor.strip():
            print("\n  [rmv_pair] Usage:")
            print("    rmv_pair <descriptor>")
            print("")
            print("  Format: pdbid_chain1_resnum1_chain2_resnum2")
            print("  Example: rmv_pair 3j6y_2S_681_2S_696")
            print("")
            print("  Residue names (A/G/C/U) are auto-detected from the structure.")
            print("")
            print("  Options:")
            print("    label_size  = font size for W/H/S labels (default: 24)")
            print("    label_color = label color (default: black)")
            print("    color1      = color for residue 1 (default: yellow)")
            print("    color2      = color for residue 2 (default: cyan)")
            print("    zoom_buffer = zoom buffer in angstroms (default: 8)")
            print("    bg_style    = background style: cartoon|none (default: cartoon)")
            print("")
            return

        visualize_pair(
            descriptor.strip(),
            label_size=int(label_size),
            label_color=label_color,
            color1=color1,
            color2=color2,
            zoom_buffer=float(zoom_buffer),
            bg_style=bg_style,
        )

    def _rmv_pair_batch(filepath="", **kwargs):
        """
        Visualize multiple RNA base pairs from a file.

        Usage:
            rmv_pair_batch /path/to/pairs.txt

        File format: one descriptor per line.
        Lines starting with # are comments.
        """
        if not filepath.strip():
            print("\n  [rmv_pair_batch] Usage:")
            print("    rmv_pair_batch /path/to/pairs.txt")
            print("")
            print("  File format: one descriptor per line.")
            print("  Example line: 3j6y_2S_681_2S_696")
            print("  Lines starting with # are treated as comments.")
            print("")
            return

        visualize_pairs_from_file(filepath.strip(), **kwargs)

    cmd.extend("rmv_pair", _rmv_pair)
    cmd.extend("rmv_pair_batch", _rmv_pair_batch)
