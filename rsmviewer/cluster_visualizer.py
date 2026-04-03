"""
RSMViewer - Cluster Visualizer Module

Provides rmv_cluster and rmv_cluster_clear commands for loading and
visualizing structural motif clusters from CSV-based cluster analysis.

Each cluster is a set of PDB entries (with residue ranges).  The
visualizer fetches each PDB, creates per-entry motif objects with
distinct colors, and prints a summary table.
"""

from pymol import cmd

# Palette of 15 distinguishable colours (PyMOL named colours)
CLUSTER_COLORS = [
    "red", "blue", "green", "orange", "magenta",
    "cyan", "yellow", "salmon", "lime", "slate",
    "hotpink", "teal", "olive", "violet", "deepteal",
]


def _get_gui():
    """Retrieve the GUI singleton (holds cluster_provider reference)."""
    from . import gui as gui_module
    return getattr(gui_module, "gui", None)


def _object_name(cluster_name, entry):
    """Build a unique PyMOL object name for one cluster entry.

    Format: cluster_<CLUSTER>_<PDB>_<CHAIN>[_<INDEX>]
    Index suffix is only added when duplicate PDB+chain exists (index > 0).
    """
    name = f"cluster_{cluster_name}_{entry.pdb_id}_{entry.chain}"
    if entry.index > 0:
        name += f"_{entry.index}"
    return name


def visualize_cluster(cluster_name):
    """Fetch PDBs, create per-entry objects, colour them, and zoom."""
    gui = _get_gui()
    if gui is None or not hasattr(gui, "cluster_provider") or gui.cluster_provider is None:
        print("\n  No cluster data loaded. Run 'rmv_db 8' first.\n")
        return

    provider = gui.cluster_provider
    entries = provider.get_cluster(cluster_name)
    if entries is None:
        available = ", ".join(provider.get_cluster_names())
        print(f"\n  Cluster '{cluster_name}' not found.")
        print(f"  Available: {available}\n")
        return

    # Determine unique PDB IDs to fetch
    pdb_ids = sorted({e.pdb_id for e in entries})

    # Fetch structures (skip those already loaded)
    existing = set(cmd.get_object_list())
    for pdb_id in pdb_ids:
        if pdb_id.upper() not in {n.upper() for n in existing}:
            print(f"  Fetching {pdb_id} …")
            cmd.fetch(pdb_id, type="cif")
            cmd.remove(f"{pdb_id} and solvent")

    # Create per-entry selection objects
    created = []
    for idx, entry in enumerate(entries):
        obj_name = _object_name(cluster_name, entry)
        color = CLUSTER_COLORS[idx % len(CLUSTER_COLORS)]

        # Build residue selection: chain X and (resi 10-15 or resi 20-25 …)
        resi_parts = " or ".join(
            f"resi {s}-{e}" for s, e in entry.regions
        )
        sele = f"{entry.pdb_id} and chain {entry.chain} and ({resi_parts})"

        cmd.create(obj_name, sele)
        cmd.color(color, obj_name)
        cmd.show("cartoon", obj_name)
        created.append((obj_name, entry, color))

    # Print summary table
    print(f"\n  ┌{'─'*70}┐")
    print(f"  │  Cluster: {cluster_name:<57} │")
    print(f"  ├{'─'*70}┤")
    print(f"  │  {'#':<4} {'Object':<34} {'PDB':<6} {'Chain':<6} {'Regions':<16} │")
    print(f"  ├{'─'*70}┤")
    for i, (obj_name, entry, color) in enumerate(created, 1):
        regions_str = ", ".join(entry.region_strings())
        print(f"  │  {i:<4} {obj_name:<34} {entry.pdb_id:<6} "
              f"{entry.chain:<6} {regions_str:<16} │")
    print(f"  └{'─'*70}┘")
    print(f"  {len(created)} object(s) created  |  color palette cycles every {len(CLUSTER_COLORS)}")
    print()

    # Zoom to the first object
    if created:
        cmd.zoom(created[0][0], buffer=10)


def clear_cluster_objects():
    """Remove all PyMOL objects whose name starts with 'cluster_'."""
    existing = cmd.get_object_list()
    removed = [n for n in existing if n.startswith("cluster_")]
    for name in removed:
        cmd.delete(name)
    print(f"  Removed {len(removed)} cluster object(s).")


def list_clusters():
    """Print available clusters from the loaded provider."""
    gui = _get_gui()
    if gui is None or not hasattr(gui, "cluster_provider") or gui.cluster_provider is None:
        print("\n  No cluster data loaded. Run 'rmv_db 8' first.\n")
        return

    provider = gui.cluster_provider
    names = provider.get_cluster_names()
    print(f"\n  Loaded clusters ({len(names)}):")
    for name in names:
        entries = provider.get_cluster(name)
        print(f"    {name:<20} ({len(entries)} entries)")
    print(f"\n  Use: rmv_cluster <NAME>  to visualize a cluster\n")


# ------------------------------------------------------------------ #
# Command registration
# ------------------------------------------------------------------ #

def register_cluster_commands():
    """Register rmv_cluster and rmv_cluster_clear commands in PyMOL."""

    def _rmv_cluster(name=""):
        """
        Load and visualize a structural motif cluster.

        Usage:
            rmv_cluster              List available clusters
            rmv_cluster ML1_2        Visualize cluster ML1_2
        """
        name = name.strip()
        if not name:
            list_clusters()
        else:
            visualize_cluster(name)

    def _rmv_cluster_clear():
        """Remove all cluster objects from the PyMOL session."""
        clear_cluster_objects()

    cmd.extend("rmv_cluster", _rmv_cluster)
    cmd.extend("rmv_cluster_clear", _rmv_cluster_clear)
