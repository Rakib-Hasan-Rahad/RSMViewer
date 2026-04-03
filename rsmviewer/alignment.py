"""
RSMViewer — Medoid-Based Superimposition Module

Provides rmv_super (sequence-independent, cmd.super) and rmv_align
(sequence-dependent, cmd.align) commands.  Both commands share the
same pipeline:

1. Resolve the user's motif name via alias table + fuzzy matching.
2. Collect individual instance objects (auto-create if needed).
3. Build an N×N pairwise RMSD matrix from temporary copies.
4. Identify the medoid (instance with minimum average RMSD).
5. Superimpose every instance onto the medoid.
6. Colour each instance uniquely (medoid = green).
7. Print a formatted report.
"""

from pymol import cmd


# ------------------------------------------------------------------ #
#  Motif-name aliasing — maps common user input variants to the
#  canonical uppercase key that appears in loaded_motifs.
# ------------------------------------------------------------------ #

MOTIF_ALIASES = {
    # K-turns
    'KTURN':            'K-TURN',
    'K_TURN':           'K-TURN',
    'KINK-TURN':        'K-TURN',
    'KINK_TURN':        'K-TURN',
    # C-loops
    'CLOOP':            'C-LOOP',
    'C_LOOP':           'C-LOOP',
    # Sarcin-ricin
    'SARCIN':           'SARCIN-RICIN',
    'SARCINRICIN':      'SARCIN-RICIN',
    'SARCIN_RICIN':     'SARCIN-RICIN',
    # Reverse K-turn
    'REVERSE-KTURN':    'REVERSE-K-TURN',
    'REVERSE_KTURN':    'REVERSE-K-TURN',
    'REVERSEKTURN':     'REVERSE-K-TURN',
    'REVERSE_K_TURN':   'REVERSE-K-TURN',
    # E-loops
    'ELOOP':            'E-LOOP',
    'E_LOOP':           'E-LOOP',
    # T-loops
    'TLOOP':            'T-LOOP',
    'T_LOOP':           'T-LOOP',
}

# Instance colours (cycled when N > len).  Medoid always gets MEDOID_COLOR.
SUPER_COLORS = [
    "red", "blue", "orange", "magenta", "cyan",
    "yellow", "salmon", "lime", "slate", "hotpink",
    "teal", "olive", "violet", "deepteal", "white",
]
MEDOID_COLOR = "green"


# ------------------------------------------------------------------ #
#  GUI access helper
# ------------------------------------------------------------------ #

def _get_gui():
    """Retrieve the GUI singleton."""
    from . import gui as gui_module
    return getattr(gui_module, "gui", None)


# ------------------------------------------------------------------ #
#  Motif name resolution
# ------------------------------------------------------------------ #

def resolve_motif_name(user_input, loaded_motifs):
    """Return list of loaded_motifs keys that match *user_input*.

    Resolution order:
    1. Exact match (case-insensitive).
    2. Canonical alias → check if alias target exists.
    3. Reverse alias → find keys that alias TO the exact match.
    4. Substring match (user input is a substring of a key, or vice-versa).

    Steps 1-3 are combined so that both K-TURN and KINK-TURN are
    returned when they co-exist in loaded_motifs and the user types
    either name.

    Returns an empty list if nothing matches.
    """
    up = user_input.upper().strip()
    keys = list(loaded_motifs.keys())
    results = set()

    # 1. Exact match
    if up in loaded_motifs:
        results.add(up)

    # 2. Forward alias: user typed an alias → canonical target
    canonical = MOTIF_ALIASES.get(up)
    if canonical and canonical in loaded_motifs:
        results.add(canonical)

    # 3. Reverse alias: find keys that alias TO the user input *or*
    #    to the canonical form found in step 2.  E.g. user types KTURN →
    #    canonical K-TURN → KINK-TURN aliases to K-TURN → add KINK-TURN.
    check_targets = {up}
    if canonical:
        check_targets.add(canonical)
    for alias, target in MOTIF_ALIASES.items():
        if target in check_targets and alias in loaded_motifs:
            results.add(alias)
        if alias in check_targets and target in loaded_motifs:
            results.add(target)

    if results:
        return sorted(results)

    # 4. Substring fallback — collect all keys where one is a substring of the other
    up_norm = up.replace('-', '').replace('_', '')
    for k in keys:
        k_norm = k.replace('-', '').replace('_', '')
        if up_norm in k_norm or k_norm in up_norm:
            results.add(k)

    return sorted(results)


# ------------------------------------------------------------------ #
#  Batch object creation  (suppress GUI updates)
# ------------------------------------------------------------------ #

def _batch_create_instance_objects(motif_type, motif_details, structure_name,
                                   source_suffix, pdb_id, indices=None,
                                   padding=0):
    """Create individual PyMOL objects for requested instances (or all).

    Returns list of created/existing object names.  Only names whose
    underlying PyMOL object *actually exists* are returned, preventing
    phantom entries from silent selector failures.

    Per-instance ``_source_suffix``, ``_pdb_id``, and ``_structure_name``
    fields in each detail dict override the function-level defaults,
    enabling cross-PDB and multi-source workflows.
    """
    from .utils.selectors import sanitize_pymol_name
    from .utils.parser import SelectionParser
    from . import colors

    existing_objects = set(cmd.get_object_list())
    created = []

    # Determine which instances to create (1-indexed)
    if indices:
        target_indices = set(indices)
    else:
        target_indices = None  # means "all"

    for i, detail in enumerate(motif_details, 1):
        if target_indices is not None and i not in target_indices:
            continue

        # Per-instance overrides (set by fetch_motif_data_action tagging)
        inst_suffix = detail.get('_source_suffix', source_suffix)
        inst_pdb = detail.get('_pdb_id', pdb_id)
        inst_struct = detail.get('_structure_name', structure_name)

        pdb_tag = f"_{inst_pdb}" if inst_pdb else ""
        obj_name = sanitize_pymol_name(
            f"{motif_type}_{i}{pdb_tag}{inst_suffix}")

        if obj_name in existing_objects:
            created.append(obj_name)
            continue

        residues = detail.get('residues', [])
        if not residues:
            continue

        chain_residues = {}
        for res in residues:
            if isinstance(res, tuple) and len(res) >= 3:
                _nuc, resi, chain = res[0], res[1], res[2]
                chain_residues.setdefault(chain, []).append(resi)

        selections = []
        for chain, resi_list in chain_residues.items():
            if padding > 0:
                expanded = set()
                for r in resi_list:
                    for offset in range(-padding, padding + 1):
                        expanded.add(r + offset)
                resi_list = sorted(expanded)
            sel = SelectionParser.create_selection_string(
                chain, sorted(resi_list))
            if sel:
                selections.append(f"({sel})")
        if not selections:
            continue

        combined_sel = " or ".join(selections)
        instance_sel = f"(model {inst_struct}) and ({combined_sel})"
        try:
            cmd.create(obj_name, instance_sel)
            # Verify the object was actually created (cmd.create can fail
            # silently — it prints a Selector-Error but doesn't raise).
            if obj_name not in cmd.get_object_list():
                continue
            cmd.show('cartoon', obj_name)
            cmd.set('cartoon_nucleic_acid_mode', 4, obj_name, quiet=1)
            cmd.set('cartoon_tube_radius', 0.4, obj_name, quiet=1)
            colors.set_motif_color_in_pymol(cmd, obj_name, motif_type)
            created.append(obj_name)
        except Exception:
            pass
        if i % 10 == 0:
            print(f"  Creating objects... {i}/{len(motif_details)}")

    cmd.rebuild()
    cmd.refresh()

    return created


# ------------------------------------------------------------------ #
#  Collect target objects for a motif type
# ------------------------------------------------------------------ #

def _collect_motif_objects(motif_type, indices=None, pdb_src_tags=None, padding=0):
    """Gather (or create) individual instance objects for *motif_type*.

    *pdb_src_tags* is an optional list of ``"PDB_SN"`` strings (e.g.
    ``["1S72_S7", "4V9F_S3"]``).  When given, only instances whose
    per-instance ``_pdb_id`` + ``_source_suffix`` match one of these
    tags are created/returned.

    Returns ``(objects, error_msg)`` — *error_msg* is ``None`` on success.
    """
    gui = _get_gui()
    if gui is None:
        return [], "GUI not initialised."

    loaded_motifs = gui.viz_manager.motif_loader.get_loaded_motifs()
    if not loaded_motifs:
        return [], "No motifs loaded. Use rmv_load_motif first."

    # Resolve aliases
    resolved_keys = resolve_motif_name(motif_type, loaded_motifs)
    if not resolved_keys:
        avail = ", ".join(sorted(loaded_motifs.keys()))
        return [], (f"Motif type '{motif_type}' not found.\n"
                    f"  Available: {avail}")

    if len(resolved_keys) > 1:
        print(f"  Alias '{motif_type}' matched multiple keys: "
              f"{', '.join(resolved_keys)}")
        print(f"  Collecting instances from all of them.")

    # Build the set of allowed PDB_SRC tags (normalised to uppercase)
    tag_set = None
    if pdb_src_tags:
        tag_set = set()
        for t in pdb_src_tags:
            # Normalise: "1S72_S7" → ("1S72", "_S7")
            upper_t = t.upper().replace(' ', '')
            if '_S' in upper_t:
                idx = upper_t.rfind('_S')
                pdb_part = upper_t[:idx]
                src_part = upper_t[idx:]      # includes "_S"
                tag_set.add((pdb_part, src_part))
            else:
                # No source suffix — match any source for this PDB
                tag_set.add((upper_t, None))

    all_objects = []

    for key in resolved_keys:
        info = loaded_motifs[key]
        motif_details = info.get('motif_details', [])
        structure_name = info.get('structure_name', '')
        source_suffix = info.get('source_suffix', '')
        pid = info.get('pdb_id', '')

        if not motif_details:
            continue

        # If tag filtering is active, pre-filter motif_details to only
        # matching instances, then create objects from just those.
        if tag_set is not None:
            filtered_details = []
            filtered_indices = []
            for i, d in enumerate(motif_details, 1):
                d_pdb = d.get('_pdb_id', pid).upper()
                d_suffix = d.get('_source_suffix', source_suffix)
                # Check if this instance matches any requested tag
                for req_pdb, req_suffix in tag_set:
                    if d_pdb == req_pdb and (req_suffix is None or d_suffix == req_suffix):
                        filtered_details.append(d)
                        filtered_indices.append(i)
                        break
            if not filtered_details:
                continue
            # If specific instance numbers were also requested, narrow further
            if indices:
                idx_set_inner = set(indices)
                filtered_indices = [i for i in filtered_indices if i in idx_set_inner]
                if not filtered_indices:
                    continue
            # Create only the matching instances (using their original indices)
            objs = _batch_create_instance_objects(
                key, motif_details, structure_name, source_suffix, pid,
                indices=filtered_indices, padding=padding)
        else:
            # Create all (or indexed subset)
            objs = _batch_create_instance_objects(
                key, motif_details, structure_name, source_suffix, pid,
                indices=indices, padding=padding)

        all_objects.extend(objs)

    if not all_objects:
        return [], f"No instance objects could be created for '{motif_type}'."

    # --- apply index filter (safety net for any code path) ---
    if indices:
        idx_set = set(indices)
        filtered = [o for o in all_objects if _extract_index(o) in idx_set]
        if not filtered:
            return [], (f"Instances {sorted(indices)} not found for "
                        f"'{motif_type}'.  Available: "
                        f"{sorted(_extract_index(o) for o in all_objects)}")
        all_objects = filtered

    return sorted(all_objects), None


def _extract_pdb(obj_name):
    """Best-effort extraction of PDB ID from an object name like
    K_TURN_3_1S72_S7  →  1S72."""
    parts = obj_name.split('_')
    # Walk backwards; source suffix is last (e.g. S7), PDB is second-to-last
    if len(parts) >= 3:
        candidate = parts[-2]
        if len(candidate) == 4 and candidate[0].isdigit():
            return candidate.upper()
    # Fallback: search all parts for a 4-char PDB-like token
    for p in reversed(parts):
        if len(p) == 4 and p[0].isdigit() and p.isalnum():
            return p.upper()
    return ""


def _extract_index(obj_name):
    """Extract the 1-based instance number from an object name like
    K_TURN_3_1S72_S7  →  3."""
    parts = obj_name.split('_')
    # Instance index is the part after the motif-type prefix, before PDB tag.
    # Motif types can have underscores (K_TURN), so we scan for the first
    # purely-numeric part that isn't the PDB ID or source suffix.
    for i, p in enumerate(parts):
        if p.isdigit():
            return int(p)
    return -1


def _extract_source(obj_name):
    """Extract the source ID from an object name like
    K_TURN_3_1S72_S7  →  7."""
    parts = obj_name.split('_')
    # Source suffix is the last part, format S<N>
    if parts:
        last = parts[-1]
        if last.startswith('S') and last[1:].isdigit():
            return int(last[1:])
    return -1


# ------------------------------------------------------------------ #
#  Pairwise RMSD matrix
# ------------------------------------------------------------------ #

def compute_pairwise_rmsd(objects, method="super"):
    """Compute NxN RMSD matrix using temporary copies.

    Returns ``(matrix, skipped)`` where *skipped* is a list of
    ``(obj_i, obj_j)`` pairs that failed.
    """
    n = len(objects)
    matrix = [[0.0] * n for _ in range(n)]
    skipped = []

    align_fn = cmd.super if method == "super" else cmd.align

    # Create temp reference copies (prefixed with _ to hide from panel)
    temps = []
    for i, obj in enumerate(objects):
        tmp = f"_medoid_ref_{i}"
        cmd.create(tmp, obj)
        temps.append(tmp)

    total = n * (n - 1) // 2
    done = 0

    for i in range(n):
        for j in range(i + 1, n):
            work = "_medoid_work"
            cmd.create(work, temps[i])
            try:
                result = align_fn(work, temps[j])
                rmsd = result[0]
                matrix[i][j] = rmsd
                matrix[j][i] = rmsd
            except Exception:
                matrix[i][j] = float('inf')
                matrix[j][i] = float('inf')
                skipped.append((objects[i], objects[j]))
            cmd.delete(work)
            done += 1
            if done % 50 == 0:
                print(f"  Computing RMSD... {done}/{total}")

    # Clean up temps
    for tmp in temps:
        cmd.delete(tmp)

    return matrix, skipped


# ------------------------------------------------------------------ #
#  Medoid selection
# ------------------------------------------------------------------ #

def find_medoid(matrix):
    """Return ``(medoid_index, avg_rmsd_list)``.

    The medoid is the row whose average RMSD (excluding ∞ pairs) is
    smallest.  Ties go to the first occurrence.
    """
    n = len(matrix)
    avg_rmsd = []
    for i in range(n):
        finite = [matrix[i][j] for j in range(n)
                  if j != i and matrix[i][j] != float('inf')]
        avg_rmsd.append(sum(finite) / len(finite) if finite else float('inf'))
    idx = avg_rmsd.index(min(avg_rmsd))
    return idx, avg_rmsd


# ------------------------------------------------------------------ #
#  Superimpose all onto medoid
# ------------------------------------------------------------------ #

def superimpose_onto_medoid(objects, medoid_idx, method="super"):
    """Move every non-medoid object onto the medoid in place.

    Returns list of ``(obj_name, rmsd, success)`` tuples.
    """
    medoid_obj = objects[medoid_idx]
    align_fn = cmd.super if method == "super" else cmd.align
    results = []
    for i, obj in enumerate(objects):
        if i == medoid_idx:
            results.append((obj, 0.0, True))
            continue
        try:
            result = align_fn(obj, medoid_obj)
            results.append((obj, result[0], True))
        except Exception:
            results.append((obj, float('inf'), False))
    return results


# ------------------------------------------------------------------ #
#  Colouring
# ------------------------------------------------------------------ #

def color_superimposed(objects, medoid_idx):
    """Assign a unique colour to each object; medoid gets green."""
    color_idx = 0
    for i, obj in enumerate(objects):
        if i == medoid_idx:
            cmd.color(MEDOID_COLOR, obj)
        else:
            cmd.color(SUPER_COLORS[color_idx % len(SUPER_COLORS)], obj)
            color_idx += 1


# ------------------------------------------------------------------ #
#  Report printer
# ------------------------------------------------------------------ #

def print_medoid_report(method, motif_type, objects, medoid_idx,
                        avg_rmsd_list, super_results, skipped):
    """Print a formatted medoid superimposition report."""
    n = len(objects)
    medoid_obj = objects[medoid_idx]

    box_w = 64
    print(f"\n  ┌{'─' * box_w}┐")
    print(f"  │  {method} — {motif_type}  ({n} instances)"
          f"{' ' * max(0, box_w - len(method) - len(motif_type) - len(str(n)) - 17)}│")
    print(f"  ├{'─' * box_w}┤")
    print(f"  │  Medoid:  {medoid_obj}"
          f"  (avg RMSD: {avg_rmsd_list[medoid_idx]:.3f} Å)"
          f"{' ' * max(0, box_w - len(medoid_obj) - 30 - len(f'{avg_rmsd_list[medoid_idx]:.3f}'))}│")
    print(f"  ├{'─' * box_w}┤")

    header = f"  │  {'#':<4} {'Instance':<30} {'Color':<10} {'RMSD to medoid':<16} │"
    print(header)

    color_idx = 0
    for i, (obj, rmsd, ok) in enumerate(super_results):
        num = i + 1
        if i == medoid_idx:
            c = MEDOID_COLOR.upper()
            rmsd_str = "— (medoid)"
        else:
            c = SUPER_COLORS[color_idx % len(SUPER_COLORS)]
            color_idx += 1
            rmsd_str = f"{rmsd:.3f} Å" if ok else "FAILED"
        line = f"  │  {num:<4} {obj:<30} {c:<10} {rmsd_str:<16} │"
        print(line)

    print(f"  ├{'─' * box_w}┤")

    # Overall average RMSD (exclude medoid-self and failed pairs)
    finite = [r for _, r, ok in super_results if ok and r > 0]
    overall = sum(finite) / len(finite) if finite else 0.0
    print(f"  │  Overall avg RMSD: {overall:.3f} Å"
          f"{' ' * max(0, box_w - 24 - len(f'{overall:.3f}'))}│")
    print(f"  │  Skipped pairs: {len(skipped)}"
          f"{' ' * max(0, box_w - 20 - len(str(len(skipped))))}│")
    print(f"  └{'─' * box_w}┘")

    if skipped:
        print(f"\n  ⚠  {len(skipped)} pair(s) skipped (superimposition failed):")
        for a, b in skipped:
            print(f"     {a} ↔ {b}")
    print()


# ------------------------------------------------------------------ #
#  CSV export helper
# ------------------------------------------------------------------ #

def _save_matrix_csv(path, objects, matrix):
    """Write the pairwise RMSD matrix to a CSV file."""
    import csv
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([""] + objects)
        for i, obj in enumerate(objects):
            writer.writerow([obj] + [f"{matrix[i][j]:.4f}" for j in range(len(objects))])
    print(f"  RMSD matrix saved to {path}")


def _get_current_pdb_src_tag():
    """Return the current PDB+source tag (e.g. ``"4V9F_S3"``) or ``None``.

    Used as the default filter when ``rmv_super`` is invoked without
    explicit PDB_SRC tags, so that only the most-recently-loaded data
    is used.
    """
    gui = _get_gui()
    if gui is None:
        return None
    pdb = getattr(gui, 'loaded_pdb_id', None)
    if not pdb:
        return None
    suffix = gui._get_source_suffix()
    if not suffix:
        return None
    return f"{pdb.upper()}{suffix}"


# ------------------------------------------------------------------ #
#  Main pipeline
# ------------------------------------------------------------------ #

def _run_medoid_pipeline(motif_type, pdb_src_tags, method, indices=None, padding=0):
    """Shared pipeline for rmv_super / rmv_align.

    *pdb_src_tags* is a list of ``"PDB_SN"`` strings (e.g.
    ``["1S72_S7", "4V9F_S3"]``) or ``None`` for all loaded instances.
    """

    # ---------- collect objects ----------
    objects, err = _collect_motif_objects(motif_type, indices=indices,
                                         pdb_src_tags=pdb_src_tags,
                                         padding=padding)
    if err:
        print(f"\n  {err}\n")
        return

    n = len(objects)
    if n < 2:
        print(f"\n  Need ≥ 2 instances for medoid analysis.  Found {n}.\n")
        return

    # Large-set warning
    pairs = n * (n - 1) // 2
    if n > 100:
        print(f"\n  {n} instances → {pairs} pairwise comparisons.  This may take a while.")
        if n > 200:
            print("  Tip: narrow with specific PDB_SRC tags.")

    # ---------- enable instance objects, hide parent structure & combined objects ----------
    all_pymol_objects = cmd.get_object_list()
    for obj in all_pymol_objects:
        upper = obj.upper()
        if "_ALL_" in upper:
            cmd.disable(obj)
        elif len(obj) == 4 and obj[0].isdigit() and obj.isalnum():
            cmd.disable(obj)

    for obj in objects:
        cmd.enable(obj)

    method_label = "rmv_super" if method == "super" else "rmv_align"
    tag_desc = ", ".join(pdb_src_tags) if pdb_src_tags else "all loaded"
    print(f"\n  [{method_label}] {motif_type}: {n} instances ({tag_desc}), {pairs} pairs")

    # ---------- pairwise RMSD ----------
    print(f"  Computing pairwise RMSD matrix ({method})...")
    matrix, skipped = compute_pairwise_rmsd(objects, method=method)

    # ---------- medoid ----------
    medoid_idx, avg_rmsd = find_medoid(matrix)

    # ---------- superimpose ----------
    print(f"  Superimposing onto medoid: {objects[medoid_idx]}")
    super_results = superimpose_onto_medoid(objects, medoid_idx, method=method)

    # ---------- colour ----------
    color_superimposed(objects, medoid_idx)

    # ---------- zoom ----------
    sel = " or ".join(objects)
    cmd.zoom(sel)

    # ---------- report ----------
    print_medoid_report(method_label, motif_type, objects, medoid_idx,
                        avg_rmsd, super_results, skipped)


# ------------------------------------------------------------------ #
#  PDB_SRC tag validation with typo detection
# ------------------------------------------------------------------ #

def _get_loaded_tags():
    """Return a set of available ``PDB_SN`` tag strings.

    Derives tags from *both* ``gui.loaded_sources`` **and** the
    per-instance ``_pdb_id`` / ``_source_suffix`` metadata in
    ``loaded_motifs``, so tags are always accurate even if the
    registry was populated before the tracking code was added.
    """
    gui = _get_gui()
    if gui is None:
        return set()
    tags = set()
    # 1. From the explicit registry
    for pdb_id, suffix in gui.loaded_sources:
        tags.add(f"{pdb_id}{suffix}")    # e.g. "1S72_S7"
    # 2. Scan loaded_motifs per-instance metadata
    loaded_motifs = gui.viz_manager.motif_loader.get_loaded_motifs()
    if loaded_motifs:
        for info in loaded_motifs.values():
            default_pdb = info.get('pdb_id', '').upper()
            default_sfx = info.get('source_suffix', '')
            for detail in info.get('motif_details', []):
                d_pdb = detail.get('_pdb_id', default_pdb).upper()
                d_sfx = detail.get('_source_suffix', default_sfx)
                if d_pdb and d_sfx:
                    tags.add(f"{d_pdb}{d_sfx}")
    return tags


def _suggest_tag(bad_tag, available_tags):
    """Find the closest matching tag for a mis-typed input.

    Checks:
    1. Case normalisation  (1s72_s7 → 1S72_S7)
    2. Missing underscore  (1S72S7  → 1S72_S7)
    3. Space instead of _  (1S72 S7 → 1S72_S7)
    4. Lowercase 's'       (1S72_s7 → 1S72_S7)
    5. PDB exists but wrong source
    6. Source exists but wrong PDB
    """
    upper = bad_tag.upper().replace(' ', '_')

    # Direct case-fix match
    if upper in available_tags:
        return upper, "case"

    # Missing underscore before S  (e.g. 1S72S3 → 1S72_S3)
    import re
    m = re.match(r'^([A-Z0-9]{4})(S\d+)$', upper)
    if m:
        fixed = f"{m.group(1)}_{m.group(2)}"
        if fixed in available_tags:
            return fixed, "missing_underscore"

    # Try adding underscore at position 4 if length is right
    if len(upper) >= 6 and '_' not in upper:
        candidate = upper[:4] + '_' + upper[4:]
        if candidate in available_tags:
            return candidate, "missing_underscore"

    # Extract PDB and source from the bad tag
    bad_pdb, bad_src = None, None
    if '_S' in upper:
        idx = upper.rfind('_S')
        bad_pdb = upper[:idx]
        bad_src = upper[idx:]
    elif len(upper) == 4:
        bad_pdb = upper

    # PDB exists but wrong source
    if bad_pdb:
        matching_pdb = [t for t in available_tags if t.startswith(bad_pdb + '_')]
        if matching_pdb:
            return matching_pdb, "wrong_source"

    # Source exists but wrong PDB
    if bad_src:
        matching_src = [t for t in available_tags if t.endswith(bad_src)]
        if matching_src:
            return matching_src, "wrong_pdb"

    return None, None


def _validate_pdb_src_tags(tags):
    """Validate user-supplied PDB_SRC tags.

    Returns ``(valid_tags, error_printed)`` where *error_printed* is
    True if an error was shown and the pipeline should abort.
    """
    available = _get_loaded_tags()

    if not available:
        print("\n  ⚠  No PDB+source combinations loaded yet.")
        print("  Load data first:")
        print("    rmv_fetch 1S72")
        print("    rmv_db 7")
        print("    rmv_load_motif\n")
        return [], True

    valid = []
    has_error = False

    for raw_tag in tags:
        upper = raw_tag.upper().replace(' ', '_')

        # Exact match
        if upper in available:
            valid.append(upper)
            continue

        # Try to suggest a fix
        suggestion, kind = _suggest_tag(raw_tag, available)

        if kind == "case":
            print(f"  ⚠  '{raw_tag}' → auto-corrected to '{suggestion}'")
            valid.append(suggestion)
            continue

        if kind == "missing_underscore":
            print(f"\n  ✗  '{raw_tag}' is not valid. Did you mean '{suggestion}'?")
            print(f"     Use underscore between PDB ID and source: {suggestion}")
            has_error = True
            continue

        if kind == "wrong_source" and isinstance(suggestion, list):
            print(f"\n  ✗  '{raw_tag}' was not loaded.")
            pdb_part = upper.split('_S')[0] if '_S' in upper else upper
            print(f"     PDB {pdb_part} was loaded with these sources:")
            for s in sorted(suggestion):
                print(f"       {s}")
            has_error = True
            continue

        if kind == "wrong_pdb" and isinstance(suggestion, list):
            print(f"\n  ✗  '{raw_tag}' was not loaded.")
            src_part = '_S' + upper.split('_S')[-1] if '_S' in upper else ''
            print(f"     Source {src_part} was loaded for these PDBs:")
            for s in sorted(suggestion):
                print(f"       {s}")
            has_error = True
            continue

        # No suggestion at all
        print(f"\n  ✗  '{raw_tag}' was not loaded and does not match any known data.")
        has_error = True

    if has_error:
        avail_sorted = sorted(available)
        print(f"\n  Currently loaded PDB+source tags:")
        for t in avail_sorted:
            print(f"    {t}")
        print(f"\n  Correct format:  rmv_super MOTIF_TYPE, PDB_SRC1, PDB_SRC2, ...")
        print(f"  Example:         rmv_super KINK-TURN, 1S72_S7, 4V9F_S3\n")
        return [], True

    return valid, False


# ------------------------------------------------------------------ #
#  Argument parser
# ------------------------------------------------------------------ #

def _parse_super_args(args, kwargs):
    """Parse rmv_super / rmv_align arguments.

    Syntax::

        rmv_super MOTIF_TYPE                        All instances
        rmv_super MOTIF_TYPE, PDB_SRC1, PDB_SRC2    Cross-PDB / multi-source
        rmv_super MOTIF_TYPE 1,3,5                   Specific instance numbers

    PyMOL splits on commas first, then spaces within each piece.
    Examples of what arrives::

        "rmv_super KINK-TURN, 1S72_S7, 4V9F_S3"
        → args = ("KINK-TURN",  "1S72_S7",  "4V9F_S3")

        "rmv_super K-TURN"
        → args = ("K-TURN",)

        "rmv_super K-TURN 1, 3, 5"
        → args = ("K-TURN 1",  "3",  "5")

    Returns ``(motif_type, pdb_src_tags, indices, padding)``.
    """
    if not args:
        return "", [], [], 0

    # Flatten: split every positional arg on whitespace
    raw_parts = []
    for a in args:
        s = str(a).strip()
        if s:
            raw_parts.extend(s.split())
    if not raw_parts:
        return "", [], [], 0

    # Classify each token
    motif_parts = []
    pdb_src_tags = []
    instance_nums = []
    padding = 0
    save_matrix = str(kwargs.get("save_matrix", "")).strip()

    for p in raw_parts:
        upper = p.upper().replace(' ', '')
        if upper.startswith("SAVE_MATRIX="):
            save_matrix = p.split("=", 1)[1].strip()
        elif upper.startswith("PADDING="):
            try:
                padding = int(p.split("=", 1)[1])
            except ValueError:
                pass
        elif _looks_like_pdb_src_tag(upper):
            pdb_src_tags.append(upper)
        elif p.isdigit():
            instance_nums.append(int(p))
        else:
            motif_parts.append(p)

    motif_type = " ".join(motif_parts) if motif_parts else ""

    return motif_type, pdb_src_tags, instance_nums, padding


def _looks_like_pdb_src_tag(token):
    """Return True if *token* looks like a PDB_SRC tag (e.g. 1S72_S7).

    Recognised patterns:
    - ``1S72_S7``   (PDB + underscore + S + digit(s))
    - ``1S72S7``    (missing underscore — still recognised for typo help)
    - ``1S72``      (bare PDB ID — only if 4-char alnum starting with digit)
    """
    import re
    # Canonical: 4-char PDB + _S + digits
    if re.match(r'^[0-9][A-Z0-9]{3}_S\d+$', token):
        return True
    # Missing underscore variant: 4-char PDB + S + digits
    if re.match(r'^[0-9][A-Z0-9]{3}S\d+$', token):
        return True
    # Bare PDB ID (only recognised if it has 4 alphanumeric chars starting with digit
    # AND there's no motif name that looks that way)
    if re.match(r'^[0-9][A-Z0-9]{3}$', token):
        return True
    return False


# ------------------------------------------------------------------ #
#  Command registration
# ------------------------------------------------------------------ #

def register_alignment_commands():
    """Register rmv_super and rmv_align commands in PyMOL."""

    def _print_usage(method_label):
        """Print usage help for rmv_super / rmv_align."""
        print(f"\n  [{method_label}] Medoid-based structural superimposition")
        print(f"  Usage:")
        print(f"    {method_label} <MOTIF_TYPE>                              All loaded instances")
        print(f"    {method_label} <MOTIF_TYPE>, <PDB_SRC1>, <PDB_SRC2>      Cross-PDB / multi-source")
        print(f"    {method_label} <MOTIF_TYPE> 1,3,5                        Specific instance numbers")
        print(f"\n  Examples:")
        print(f"    {method_label} K-TURN                           All K-TURN from all loaded data")
        print(f"    {method_label} KINK-TURN, 1S72_S7, 4V9F_S3      Only 1S72 source 7 + 4V9F source 3")
        print(f"    {method_label} K-TURN, 1S72_S7, 1S72_S3         Same PDB, compare 2 sources")
        print(f"    {method_label} K-TURN 2,5,8                     Instance numbers 2, 5, 8 only")
        # Show currently loaded tags if any
        tags = sorted(_get_loaded_tags())
        if tags:
            print(f"\n  Currently loaded PDB_SRC tags:")
            for t in tags:
                print(f"    {t}")
        print()

    def _rmv_super(*args, **kwargs):
        """
        Medoid-based superimposition (sequence-independent, uses cmd.super).

        Usage:
            rmv_super <MOTIF_TYPE>                              All instances
            rmv_super <MOTIF_TYPE>, <PDB_SRC1>, <PDB_SRC2>      Cross-PDB / multi-source
            rmv_super <MOTIF_TYPE> 1,3,5                        Specific instances
        """
        motif_type, pdb_src_tags, indices, padding = _parse_super_args(args, kwargs)
        if not motif_type:
            _print_usage("rmv_super")
            return

        # Validate PDB_SRC tags if provided
        if pdb_src_tags:
            valid_tags, had_error = _validate_pdb_src_tags(pdb_src_tags)
            if had_error:
                return
            pdb_src_tags = valid_tags
        else:
            # No explicit tags → default to current PDB+source
            default_tag = _get_current_pdb_src_tag()
            if default_tag:
                pdb_src_tags = [default_tag]

        _run_medoid_pipeline(
            motif_type,
            pdb_src_tags=pdb_src_tags or None,
            method="super",
            indices=indices or None,
            padding=padding)

    def _rmv_align(*args, **kwargs):
        """
        Medoid-based superimposition (sequence-dependent, uses cmd.align).

        Usage:
            rmv_align <MOTIF_TYPE>                              All instances
            rmv_align <MOTIF_TYPE>, <PDB_SRC1>, <PDB_SRC2>      Cross-PDB / multi-source
            rmv_align <MOTIF_TYPE> 1,3,5                        Specific instances
        """
        motif_type, pdb_src_tags, indices, padding = _parse_super_args(args, kwargs)
        if not motif_type:
            _print_usage("rmv_align")
            return

        # Validate PDB_SRC tags if provided
        if pdb_src_tags:
            valid_tags, had_error = _validate_pdb_src_tags(pdb_src_tags)
            if had_error:
                return
            pdb_src_tags = valid_tags
        else:
            # No explicit tags → default to current PDB+source
            default_tag = _get_current_pdb_src_tag()
            if default_tag:
                pdb_src_tags = [default_tag]

        _run_medoid_pipeline(
            motif_type,
            pdb_src_tags=pdb_src_tags or None,
            method="align",
            indices=indices or None,
            padding=padding)

    cmd.extend("rmv_super", _rmv_super)
    cmd.extend("rmv_align", _rmv_align)
