"""
RSMViewer - GUI Module
Provides PyMOL GUI interface for the plugin with multi-database support.

This module provides:
- MotifVisualizerGUI: Main GUI class for the plugin
- PyMOL command registration (24 commands)
- Database selection and switching functionality
- Multi-source loading (kept separate; merging deferred to rmv_select/rmv_combine_groups)

Author: CBB LAB @Rakib Hasan Rahad
Version: 2.0.0
Last Updated: 10 September 2026
"""

from typing import List, Optional, Dict, Tuple
import glob
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
import urllib.request
from pymol import cmd
from .loader import VisualizationManager
from .utils import get_logger
from . import colors
from .database import get_registry
from .database.config import SOURCE_ID_MAP
from .database.consolidated_table import ConsolidatedAnnotationTable, DEFAULT_JACCARD_THRESHOLD
from .database.query_parser import QuerySyntaxError, parse_selection_query
from pathlib import Path
import platform
from .tools.runtime_layout import discover_rmsx_runtime_dir, get_runtime_platform_dir, validate_fr3d_checkout


 # ---------------------------------------------------------------------------
 # Semantic normalisation constants (shared between Atlas & BGSU pipelines)
 # ---------------------------------------------------------------------------
_SEMANTIC_PATTERNS = {
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

_LOOP_TYPE_NAMES = {
    'HL': 'Hairpin Loop (HL)',
    'IL': 'Internal Loop (IL)',
    'J3': '3-way Junction (J3)',
    'J4': '4-way Junction (J4)',
    'J5': '5-way Junction (J5)',
    'J6': '6-way Junction (J6)',
    'J7': '7-way Junction (J7)',
    'J8': '8-way Junction (J8)',
}

_GENERIC_TYPES = set(_LOOP_TYPE_NAMES.keys())

# Compact source names for table headers - SOURCE_ID_MAP's full names (e.g.
# "BGSU RNA 3D Hub", "RNAMotifScanX (RMSX)") make headers hard to read.
_SHORT_SOURCE_NAMES = {
    1: 'ATLAS', 2: 'RFAM', 3: 'BGSU', 4: 'RFAM-API',
    5: 'FR3D', 6: 'RMS', 7: 'RMSX', 8: 'NoBIAS',
}


def _short_source_name(sid: int) -> str:
    return _SHORT_SOURCE_NAMES.get(sid, SOURCE_ID_MAP.get(sid, {}).get('name', f'S{sid}'))


def _categorize_by_annotation(annotation, loop_type):
    """Determine semantic category from annotation, falling back to loop type name."""
    if annotation:
        annotation_lower = annotation.lower()
        for category, patterns in _SEMANTIC_PATTERNS.items():
            for pattern in patterns:
                if pattern.lower() in annotation_lower:
                    return category
        # Has annotation but no pattern match – use raw annotation as category
        return annotation
    # No annotation – descriptive fallback
    return _LOOP_TYPE_NAMES.get(loop_type, loop_type)


def _is_generic_equivalent(generic: str, final: str) -> bool:
    """True when *final* is just the same generic type reworded (e.g. 'IL' vs
    'Internal Loop (IL)'), not a genuine semantic refinement worth storing
    as a second hierarchy level.
    """
    g = generic.strip().upper()
    f = final.strip().upper()
    if g == f:
        return True
    return _LOOP_TYPE_NAMES.get(g, '').upper() == f


def _normalize_motif_groups(available_motifs):
    """Re-categorise generic HL / IL / J3… keys to semantic types using
    per-instance ``annotation`` (the same logic BGSU API uses).

    Non-generic keys (already semantic) pass through unchanged.
    Returns a *new* dict; the original is not mutated.
    """
    normalized = {}
    for motif_type, instances in available_motifs.items():
        key_upper = motif_type.upper()
        if key_upper not in _GENERIC_TYPES:
            normalized.setdefault(motif_type, []).extend(instances)
            continue
        for inst in instances:
            ann = getattr(inst, 'annotation', '') or ''
            category = _categorize_by_annotation(ann, key_upper)
            # Remember the generic parent type so the hierarchy cache can
            # store both levels (e.g. IL -> K-TURN) instead of only the
            # final, more specific category.
            if getattr(inst, 'metadata', None) is None:
                inst.metadata = {}
            inst.metadata.setdefault('_generic_type', key_upper)
            normalized.setdefault(category, []).append(inst)
    return normalized


def _parse_db_expression(tokens: List[str]) -> Optional[List[Tuple[str, int]]]:
    """Parse tokens like ['db','3','and','db','2'] into [('AND',3),('AND',2)].

    Left-to-right evaluation, no operator precedence beyond order; 'not'
    applies only to the 'db N' immediately following it. Returns None if
    the tokens don't match the db-expression grammar at all.
    """
    ops: List[Tuple[str, int]] = []
    i = 0
    op = 'AND'
    while i < len(tokens):
        t = tokens[i].lower()
        if t in ('and', 'or'):
            op = t.upper()
            i += 1
            continue
        if t == 'not':
            op = 'NOT'
            i += 1
            continue
        if t == 'db' and i + 1 < len(tokens) and tokens[i + 1].isdigit():
            ops.append((op, int(tokens[i + 1])))
            op = 'AND'
            i += 2
            continue
        return None
    return ops or None


def _split_motif_and_db_expression(parts: List[str]) -> Tuple[str, Optional[List[Tuple[str, int]]]]:
    """Split ['K-TURN','db','3','and','db','2'] into ('K-TURN', [ops...]).

    Also handles a leading 'not' immediately before the expression, e.g.
    ['SR','not','db','3','and','db','7'] -> ('SR', [('NOT',3),('AND',7)]) -
    anchoring on the first 'db' alone would otherwise strand that 'not' in
    the motif filter text and silently drop it from the expression.

    Returns (motif_filter, None) when no db-expression keyword is present
    at all, so callers fall back to the existing plain motif-name behaviour.
    """
    lowered = [p.lower() for p in parts]
    keywords = {'db', 'and', 'or', 'not'}
    idx = next((i for i, t in enumerate(lowered) if t in keywords), None)
    if idx is None:
        return ' '.join(parts), None
    ops = _parse_db_expression(parts[idx:])
    if ops is None:
        return ' '.join(parts), None
    return ' '.join(parts[:idx]).strip(), ops


_ALIAS_NAME_RE = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')


def _extract_trailing_alias(parts: List[str]) -> Tuple[List[str], Optional[str], Optional[str]]:
    """Strip a trailing 'as ALIAS' clause (Applications 5/6 benchmarking).

    Returns (remaining_parts, alias, error_message). *alias* is None when
    no 'as' token is present. *error_message* is set (and alias is None)
    when 'as' is present but not followed by exactly one valid name.
    """
    lowered = [p.lower() for p in parts]
    if 'as' not in lowered:
        return parts, None, None
    idx = lowered.index('as')
    remaining = parts[:idx]
    alias_tokens = parts[idx + 1:]
    if len(alias_tokens) != 1 or not _ALIAS_NAME_RE.match(alias_tokens[0]):
        return parts, None, (
            "Usage: ... as <ALIAS_NAME>  (single word: letters, digits, "
            "underscore only, e.g. 'as novel_SR')"
        )
    return remaining, alias_tokens[0], None


def _normalize_motif_alias_text(text: str) -> str:
    """Normalize a motif name/filter for alias-tolerant comparison, so
    'K-TURN', 'Kink-Turn' and 'KINK_TURN' all reduce to the same token
    (reuses alignment.py's canonical alias table).
    """
    norm = ''.join(ch for ch in text.upper() if ch.isalnum())
    try:
        from .alignment import MOTIF_ALIASES
    except Exception:
        MOTIF_ALIASES = {}
    for alias_key, canonical in MOTIF_ALIASES.items():
        alias_norm = ''.join(ch for ch in alias_key.upper() if ch.isalnum())
        if norm == alias_norm:
            return ''.join(ch for ch in canonical.upper() if ch.isalnum())
    return norm


def _format_db_expression(ops: List[Tuple[str, int]]) -> str:
    """Render parsed db-expression ops back to human-readable text for storage."""
    parts = []
    for i, (op, sid) in enumerate(ops):
        if i == 0:
            parts.append(f"not db {sid}" if op == 'NOT' else f"db {sid}")
        elif op == 'AND':
            parts.append(f"and db {sid}")
        elif op == 'OR':
            parts.append(f"or db {sid}")
        else:
            parts.append(f"and not db {sid}")
    return ' '.join(parts)


def _partition_into_known_types(all_parts: List[str], loaded_motifs: dict, resolve_fn) -> Optional[List[str]]:
    """Try to split a flat token list into 2+ distinct known motif types
    (Application 3: 'rmv_show K-TURN, Sarcin-Ricin, HL' arrives as flat
    tokens ['K-TURN','Sarcin-Ricin','HL'] with no comma marker left).

    Greedily consumes the LONGEST matching run of tokens at each position
    (so multi-word names like '4-WAY JUNCTION (J4)' are still recognised
    inside a list). Returns None unless the whole token list is consumed
    and at least two distinct types are found - anything else (a single
    multi-word type, or instance-number tokens) falls through unchanged
    to the existing single-type logic. *resolve_fn* is the GUI's
    ``_resolve_loaded_motif_type`` bound method (alias/substring lookup).
    """
    if len(all_parts) < 2:
        return None
    resolved: List[str] = []
    i = 0
    n = len(all_parts)
    while i < n:
        match = None
        for j in range(n, i, -1):
            candidate = " ".join(all_parts[i:j])
            resolved_type = resolve_fn(candidate, loaded_motifs)
            if resolved_type in loaded_motifs:
                match = (j, resolved_type)
                break
        if not match:
            return None
        j, resolved_type = match
        resolved.append(resolved_type)
        i = j
    if len(resolved) < 2:
        return None
    return resolved


def _format_residue_key_ranges(r_key: str) -> str:
    """Render a 'chain:num;chain:num;...' residue key as 'A:10-12,B:5'."""
    by_chain: Dict[str, List[int]] = {}
    for token in r_key.split(';'):
        if not token:
            continue
        chain, _, num_str = token.partition(':')
        try:
            by_chain.setdefault(chain, []).append(int(num_str))
        except ValueError:
            continue
    parts = []
    for chain in sorted(by_chain):
        nums = sorted(by_chain[chain])
        ranges = []
        start = prev = nums[0]
        for n in nums[1:]:
            if n == prev + 1:
                prev = n
                continue
            ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
            start = prev = n
        ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
        parts.append(f"{chain}:{','.join(ranges)}")
    return "; ".join(parts)


def _evaluate_db_expression(ops: List[Tuple[str, int]], present_ids: set) -> bool:
    """Left-to-right AND/OR/NOT evaluation of a parsed db-expression.

    No seed value like True/False is used for the very first term: seeding
    with True made a leading 'OR' always evaluate to True regardless of
    whether that source was actually present (True or X == True), silently
    matching every row. The first term instead just becomes the result
    (negated if it's NOT); only later terms combine with AND/OR/NOT.
    """
    result = None
    for op, sid in ops:
        present = sid in present_ids
        if result is None:
            result = (not present) if op == 'NOT' else present
            continue
        if op == 'AND':
            result = result and present
        elif op == 'OR':
            result = result or present
        elif op == 'NOT':
            result = result and (not present)
    return True if result is None else result


def _parse_residue_key(r_key: str):
    """Parse a 'chain:num;chain:num;...' hierarchy key into a set of (chain, num)."""
    pairs = set()
    for token in r_key.split(';'):
        if not token:
            continue
        chain, _, num_str = token.partition(':')
        try:
            pairs.add((chain, int(num_str)))
        except ValueError:
            continue
    return pairs


def _cluster_hierarchy_entries(hierarchy: Dict[str, List[Dict]], jaccard_threshold: float) -> List[Dict]:
    """Group hierarchy residue-keys into instances by residue set, not by name.

    Two residue-keys are the same physical instance when they are an exact
    match, one is a subset of the other, or their Jaccard similarity is at
    or above ``jaccard_threshold`` — the same residue-overlap rule used during
    multi-source annotation merging.
    """
    from .database.residue_merger import _jaccard

    keys = list(hierarchy.keys())
    parsed = {k: _parse_residue_key(k) for k in keys}
    parent = {k: k for k in keys}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(len(keys)):
        set_i = parsed[keys[i]]
        if not set_i:
            continue
        for j in range(i + 1, len(keys)):
            set_j = parsed[keys[j]]
            if not set_j:
                continue
            if (set_i == set_j or set_i.issubset(set_j) or set_j.issubset(set_i)
                    or _jaccard(set_i, set_j) >= jaccard_threshold):
                union(keys[i], keys[j])

    groups: Dict[str, List[str]] = {}
    for k in keys:
        groups.setdefault(find(k), []).append(k)

    clusters = []
    for member_keys in groups.values():
        # source_key -> {hierarchy_level: motif_label}, so a single source's
        # own generic->specific chain (e.g. BGSU IL -> K-TURN) is kept as
        # separate ordered levels instead of collapsing to one label.
        levels_by_source: Dict[str, Dict[int, str]] = {}
        custom_ids = []
        combined = set()
        for mk in member_keys:
            combined |= parsed[mk]
            for entry in hierarchy[mk]:
                levels = levels_by_source.setdefault(entry['source_key'], {})
                levels.setdefault(entry.get('hierarchy_level', 1), entry['motif_label'])
                custom_ids.append(entry['custom_id'])
        labels_by_source: Dict[str, List[str]] = {
            sid: [levels[lvl] for lvl in sorted(levels)]
            for sid, levels in levels_by_source.items()
        }
        clusters.append({
            'custom_id': sorted(custom_ids)[0] if custom_ids else '',
            'residues': combined,
            'labels_by_source': labels_by_source,
        })
    return clusters


_INSTANCE_RANGE_RE = re.compile(r'^(\d+)-(\d+)$')


def _expand_instance_token(token: str) -> Optional[List[int]]:
    """Expand an instance-selector token: '5' -> [5], '1-5' -> [1..5], '1,3,5' -> [1,3,5]."""
    token = token.strip()
    if not token:
        return None
    nums: List[int] = []
    for piece in token.split(','):
        piece = piece.strip()
        if not piece:
            continue
        if piece.isdigit():
            nums.append(int(piece))
            continue
        m = _INSTANCE_RANGE_RE.match(piece)
        if not m:
            return None
        lo, hi = int(m.group(1)), int(m.group(2))
        if lo > hi:
            lo, hi = hi, lo
        nums.extend(range(lo, hi + 1))
    return nums or None


class MotifVisualizerGUI:
    """PyMOL GUI for RNA motif visualization with multi-database support."""

    def command_error(self, message: str) -> None:
        """Report a command/input error with the central command reference."""
        self.logger.error(message)
        self.logger.info("Run rmv_help for command syntax and available options.")

    def __init__(self):
        """Initialize GUI components."""
        self.logger = get_logger()

        # Get path to motif database
        plugin_dir = Path(__file__).parent
        distribution_dir = plugin_dir.parent
        self.database_dir = plugin_dir / 'motif_database'

        # Initialize visualization manager
        self.viz_manager = VisualizationManager(cmd, str(self.database_dir))

        # Track UI state
        self.motif_visibility = {}

        # Track currently loaded PDB
        self.loaded_pdb = None
        self.loaded_pdb_id = None
        self.loaded_structures: Dict[str, str] = {}
        self.annotation_tables: Dict[str, ConsolidatedAnnotationTable] = {}
        self.query_groups: Dict[str, dict] = {}

        # Track current source mode
        self.current_source_mode = None
        self.current_source_names: List[str] = []
        self.current_user_tool = None
        self.current_local_source = None      # 'atlas', 'rfam', or None (for both)
        self.current_web_source = None         # 'bgsu', 'rfam', or None (for auto)
        self.combined_source_ids = []          # List of source IDs for combining
        self.current_source_id = None          # Numeric source ID (1-8) for object tagging

        # Per-PDB source memory: pdb_id -> {'mode', 'tool', 'jaccard'} raw
        # rmv_db args, so switching back to an already-loaded PDB (multi-PDB
        # session) auto-restores the source it was last loaded from.
        self.pdb_source_state: Dict[str, dict] = {}

        # Track filtering state for RMS, RMSX, and NoBIAS (for user annotations)
        self.user_rms_filtering_enabled = True   # Default: filters ON
        self.user_rmsx_filtering_enabled = True  # Default: filters ON
        self.user_nobias_filtering_enabled = True  # Default: filters ON

        # RMSX thresholds are loaded from config/rmsx_config.json only.
        self.user_rmsx_custom_pvalues = {}

        # Track chain ID convention: 1 = auth_asym_id (default), 0 = label_asym_id
        self.cif_use_auth = 1
        self.auth_to_label_map = {}   # auth_asym_id -> label_asym_id chain mapping

        # Track custom user annotation data paths per source (user-source data paths)
        self.user_data_paths = {}  # Dict: {source_id_int: path_str}

        # Jaccard similarity threshold for annotation merging (0.0–1.0)
        self.jaccard_threshold = DEFAULT_JACCARD_THRESHOLD

        # Store within-source dedup stats for display in source attribution report
        # Dict: {source_id: (before_count, after_count)}
        self.dedup_stats = {}

        # Track loaded PDB+source combinations for cross-PDB superimposition
        # Each entry is a tuple (pdb_id_upper, source_suffix) e.g. ("1S72", "_S7")
        self.loaded_sources = set()

        # FR3D (Source 5): minimal wrapper around an EXTERNAL, user-installed
        # official BGSU fr3d-python checkout. Registered per-session via
        # `rmv_db FR3D` (bundled default config) or `rmv_fr3d register <config>`.
        # vendored, patched, or checked out on the user's behalf.
        fr3d_dir = plugin_dir / 'database' / 'user_annotations' / 'fr3d'
        self.default_fr3d_config = distribution_dir / 'config' / 'fr3d_config.json'
        self.fr3d_cache_path = str((distribution_dir / 'external' / 'fr3d' / 'fr3d_cache').resolve())
        self.fr3d_output_dir = str((distribution_dir / 'output' / 'fr3d_runs').resolve())
        self._fr3d_reset_state()

        # RNAMotifScanX wrapper runtime settings (persisted in user_annotations/RNAMotifScanX)
        rmsx_dir = plugin_dir / 'database' / 'user_annotations' / 'RNAMotifScanX'
        self.rmsx_config_file = rmsx_dir / '.rmsx_wrapper.json'
        self.rmsx_executable_path = ''
        self.rmsx_output_path = str((distribution_dir / 'output' / 'rmsx_results').resolve())
        self.rmsx_working_dir = ''
        self.rmsx_args_template = ''
        self.rmsx_auto_run_on_fetch = False
        # RMSX runtime config (integrated source-7 runtime)
        self.rmsx_pipeline_config_path = str((distribution_dir / 'config' / 'rmsx_config.json').resolve())
        self.rmsx_pipeline_config = {}
        self.rmsx_runtime_dir = discover_rmsx_runtime_dir(
            str((distribution_dir / 'external' / 'rmsx').resolve()),
            [str((distribution_dir / 'external' / 'rmsx').resolve())],
        )
        if not self.rmsx_runtime_dir:
            self.rmsx_runtime_dir = str((distribution_dir / 'external' / 'rmsx').resolve())
        self.rmsx_setup_attempted = False
        self._load_rmsx_wrapper_config()

    def _normalize_filter_motif_name(self, motif_name: str) -> str:
        """Normalize user-entered motif names for RMS/RMSX/NoBIAS P-value maps.

        Accepts common aliases and consensus suffix variants, and returns a
        canonical uppercase family key used by converters.
        """
        key = str(motif_name or '').strip().upper()
        if not key:
            return key

        key = key.replace(' ', '-').replace('_', '-')
        if key.endswith('-CONSENSUS'):
            key = key[:-10]

        aliases = {
            'KTURN': 'K-TURN',
            'KINK-TURN': 'K-TURN',
            'CLOOP': 'C-LOOP',
            'SARCIN': 'SARCIN-RICIN',
            'SARCINRICIN': 'SARCIN-RICIN',
            'REVERSE-KTURN': 'REVERSE-K-TURN',
            'REVERSEKTURN': 'REVERSE-K-TURN',
            'ELOOP': 'E-LOOP',
        }

        no_dash = key.replace('-', '')
        return aliases.get(key, aliases.get(no_dash, key))

    def _is_placeholder_path(self, value: str) -> bool:
        """Return True when a config path is clearly a template placeholder."""
        text = str(value or '').strip()
        if not text:
            return False
        upper = text.upper()
        return ('ABSOLUTE/PATH' in upper) or ('/PATH/TO/' in upper) or upper.startswith('/ABSOLUTE/')

    def _find_first_existing_path(self, candidates: List[Path]) -> str:
        """Return first existing candidate path or empty string."""
        for candidate in candidates:
            if candidate.exists() and candidate.is_file():
                return str(candidate.resolve())
        return ''

    def _prepare_local_pdb_for_rmsx(self, pdb_id: str, output_dir: str, force_refresh: bool = False) -> str:
        """Export a local PDB from the loaded PyMOL object for RMSX preprocessing."""
        try:
            pdb_upper = str(pdb_id or '').strip().upper()
            if not pdb_upper:
                return ''
            os.makedirs(output_dir, exist_ok=True)
            pdb_out = os.path.join(output_dir, f"{pdb_upper}.pdb")
            if os.path.isfile(pdb_out) and not force_refresh:
                return pdb_out

            object_names = set(cmd.get_names('objects'))
            obj_candidates = [
                str(getattr(self, 'loaded_pdb', '') or '').strip(),
                pdb_upper.lower(),
                pdb_upper,
            ]
            for obj_name in obj_candidates:
                if obj_name and obj_name in object_names:
                    cmd.save(pdb_out, obj_name)
                    if os.path.isfile(pdb_out):
                        self.logger.debug(f"Exported local PDB for RMSX preprocessing: {pdb_out}")
                        return pdb_out
            return ''
        except Exception as exc:
            self.logger.debug(f"RMSX local PDB export step skipped: {exc}")
            return ''

    def _build_internal_rmsx_config(self) -> Dict:
        """Build Source-7 RMSX runtime config from bundled plugin assets."""
        runtime_dir = Path(self.rmsx_runtime_dir)
        bin_root = runtime_dir / 'bin'
        queries_dir = runtime_dir / 'queries'

        platform_dir = get_runtime_platform_dir()

        platform_bin_dir = bin_root / platform_dir
        candidate_bin_dirs = []
        if platform_bin_dir.exists():
            candidate_bin_dirs.append(platform_bin_dir)
        if bin_root.exists():
            for p in sorted(bin_root.iterdir()):
                if p.is_dir() and p not in candidate_bin_dirs:
                    candidate_bin_dirs.append(p)

        rmsx_exe = ''
        mc_exe = ''
        for bdir in candidate_bin_dirs:
            if not rmsx_exe:
                rmsx_exe = self._find_first_existing_path([
                    bdir / 'RNAMotifScanX',
                    bdir / 'scan',
                    bdir / 'RNAMotifScanX.exe',
                    bdir / 'scan.exe',
                ])
            if not mc_exe:
                mc_exe = self._find_first_existing_path([
                    bdir / 'MC-Annotate',
                    bdir / 'mc-annotate',
                    bdir / 'MC-Annotate.exe',
                    bdir / 'mc-annotate.exe',
                ])
            if rmsx_exe and mc_exe:
                break

        # RNAVIEW: reference preprocessing merges MC-Annotate with RNAVIEW.
        # Discover a compiled rnaview binary (if the runtime has been built) and
        # the bundled RNAVIEW base directory that holds the BASEPARS resources.
        rnaview_exe = ''
        for bdir in candidate_bin_dirs:
            rnaview_exe = self._find_first_existing_path([
                bdir / 'rnaview',
                bdir / 'rnaview.exe',
            ])
            if rnaview_exe:
                break
        rnaview_base = runtime_dir / 'src' / 'RNAMotifScanX_src' / 'ThirdParty' / 'RNAVIEW'
        rnaview_dir = str(rnaview_base.resolve()) if (rnaview_base / 'BASEPARS').exists() else ''
        if not rnaview_exe:
            bundled_rnaview = rnaview_base / 'bin' / 'rnaview'
            if bundled_rnaview.exists():
                rnaview_exe = str(bundled_rnaview.resolve())

        prebuild_archive = ''
        archive_candidates = [
            runtime_dir.parent.parent / 'database' / 'user_annotations' / 'RNAMotifScanX' / 'PDB_prebuild.tgz',
            runtime_dir.parent.parent / 'RNAMotifScanX' / 'PDB_prebuild.tgz',
            runtime_dir.parent / 'RNAMotifScanX' / 'PDB_prebuild.tgz',
        ]
        for candidate in archive_candidates:
            if candidate.exists() and candidate.is_file():
                prebuild_archive = str(candidate.resolve())
                break

        cfg = {
            'data_mode': 'preannotated',
            'jaccard_threshold': self.jaccard_threshold,
            'rmsx_executable': rmsx_exe,
            'mc_annotate_executable': mc_exe,
            'rnaview_executable': rnaview_exe,
            'rnaview_dir': rnaview_dir,
            'pdb_prebuild_archive': prebuild_archive,
            'pdb_prebuild_dir': str((Path(__file__).parent.parent / 'external' / 'rmsx_preannotated' / 'rmsx_work_default').resolve()),
            'incorporate_rnaview': True,
            'rnaview_required': False,
            'query_motifs_dir': str(queries_dir.resolve()) if queries_dir.exists() else '',
            'cif_input_dir': '',
            'output_dir': self.rmsx_output_path,
            'auto_download_cif': False,
            'motif_families': ['k-turn', 'c-loop', 'sarcin-ricin', 'reverse-kturn', 'e-loop'],
            'target_chains': ['0'],
            'max_strands': 3,
            'num_threads': 4,
        }
        self._apply_external_rmsx_config(cfg)
        return cfg

    def _apply_external_rmsx_config(self, cfg: Dict) -> None:
        """Merge a user-supplied rmsx_config.json (path/pipeline overrides + P-value cutoffs) into cfg in-place.

        Pipeline path/option keys (anything except '_'-prefixed metadata and
        'pvalue_thresholds') override the bundled runtime defaults when present
        and non-empty. 'pvalue_thresholds' seeds self.user_rmsx_custom_pvalues
        so Source 7 filtering picks up the user's cutoffs without requiring
        'rmv_db 7 MOTIF value' for every family.
        """
        path = str(getattr(self, 'rmsx_pipeline_config_path', '') or '').strip()
        if not path or not os.path.isfile(path):
            return
        try:
            with open(path, 'r', encoding='utf-8') as fh:
                external = json.load(fh)
        except Exception as e:
            self.logger.warning(f"Could not read RMSX config {path}: {type(e).__name__}: {e}")
            return

        config_dir = Path(path).resolve().parent
        path_keys = {
            'rmsx_executable', 'mc_annotate_executable', 'rnaview_executable',
            'rnaview_dir', 'pdb_prebuild_archive', 'pdb_prebuild_dir', 'query_motifs_dir',
            'cif_input_dir', 'output_dir',
        }
        for key, value in external.items():
            if key.startswith('_') or key == 'pvalue_thresholds':
                continue
            if value not in (None, '', []):
                if key in path_keys and isinstance(value, str):
                    configured_path = Path(os.path.expanduser(value))
                    if not configured_path.is_absolute():
                        value = str((config_dir / configured_path).resolve())
                cfg[key] = value

        pvalue_thresholds = external.get('pvalue_thresholds')
        if isinstance(pvalue_thresholds, dict) and pvalue_thresholds:
            overrides = {}
            for motif_name, pvalue in pvalue_thresholds.items():
                if motif_name.startswith('_'):
                    continue
                try:
                    overrides[self._normalize_filter_motif_name(motif_name)] = float(pvalue)
                except (TypeError, ValueError):
                    self.logger.warning(f"Ignoring invalid P-value for '{motif_name}' in {path}")
            if overrides:
                self.user_rmsx_custom_pvalues = overrides
                self.logger.debug(f"Loaded {len(overrides)} P-value cutoff(s) from {path}")


    def _run_rmsx_runtime_setup(self, build: bool = False) -> Dict:
        """Run integrated RMSX runtime doctor/setup script and return parsed report."""
        setup_script = Path(self.rmsx_runtime_dir) / 'setup_runtime.py'
        if not setup_script.exists():
            return {
                'ok': False,
                'missing': ['setup_runtime.py'],
                'setup_message': f'Missing setup script: {setup_script}',
            }

        cmdline = [
            sys.executable,
            str(setup_script),
            '--runtime-dir',
            self.rmsx_runtime_dir,
            '--json',
        ]
        query_file = str(getattr(self, 'rmsx_query_file', '') or '').strip()
        if query_file:
            cmdline.extend(['--query-file', query_file])
        if build:
            cmdline.append('--build')

        try:
            proc = subprocess.run(cmdline, capture_output=True, text=True, check=False)
            output = (proc.stdout or '').strip()
            if not output:
                return {
                    'ok': False,
                    'missing': ['runtime setup output'],
                    'setup_message': (proc.stderr or 'runtime setup produced no output').strip(),
                }
            report = json.loads(output)
            return report
        except Exception as e:
            return {
                'ok': False,
                'missing': ['runtime setup execution'],
                'setup_message': f'Failed to run setup script: {type(e).__name__}: {e}',
            }

    def ensure_rmsx_runtime_ready(self, auto_setup: bool = True) -> bool:
        """Ensure integrated Source-7 runtime is present. Attempts setup once per session."""
        report = self._run_rmsx_runtime_setup(build=False)
        if report.get('ok'):
            return True

        if auto_setup and not self.rmsx_setup_attempted:
            self.rmsx_setup_attempted = True
            self.logger.info('RMSX runtime missing; attempting first-run setup...')
            report = self._run_rmsx_runtime_setup(build=True)
            if report.get('ok'):
                self.logger.success('RMSX runtime setup completed successfully')
                return True

        missing = report.get('missing', []) or []
        if missing:
            self.logger.error(f"RMSX runtime is incomplete: {', '.join(str(x) for x in missing)}")
        setup_message = str(report.get('setup_message', '') or '').strip()
        if setup_message:
            self.logger.error(setup_message)
        self.logger.info('Run: rmv_rmsx_doctor')
        return False

    def rmsx_doctor(self, auto_setup: bool = False):
        """Print integrated Source-7 runtime diagnostics."""
        report = self._run_rmsx_runtime_setup(build=auto_setup)

        print('\n' + '=' * 70)
        print('RMSX Doctor')
        print('=' * 70)
        print(f"Runtime dir     : {self.rmsx_runtime_dir}")
        print(f"Platform dir    : {report.get('platform_dir', '(unknown)')}")
        print(f"Binary dir      : {report.get('bin_dir', '(unknown)')}")
        print(f"RMSX executable : {report.get('rmsx_executable', '(missing)') or '(missing)'}")
        print(f"MC-Annotate     : {report.get('mc_annotate_executable', '(missing)') or '(missing)'}")
        print(f"Status          : {'OK' if report.get('ok') else 'NOT READY'}")

        missing = report.get('missing', []) or []
        if missing:
            print('\nMissing items:')
            for item in missing:
                print(f"  - {item}")

        mq = report.get('missing_queries', []) or []
        if mq:
            print('\nMissing query files:')
            for item in mq:
                print(f"  - {item}")

        setup_message = str(report.get('setup_message', '') or '').strip()
        if setup_message:
            print(f"\nSetup message: {setup_message}")

        print('\nHints:')
        print('  - Put platform binaries into rsmviewer/tools/rmsx_runtime/bin/<platform>')
        print('  - Then run rmv_rmsx run_current or rmv_load_motif with source 7 active')
        print('=' * 70 + '\n')

    def _get_source_suffix(self):
        """Get a source suffix for PyMOL object naming.

        Named workflows use canonical textual source keys. The numeric form
        remains only as a compatibility fallback for legacy internal calls.
        """
        if self.current_source_names:
            safe_names = [
                re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
                for name in self.current_source_names
            ]
            return "_SRC_" + "_".join(safe_names)
        if self.current_source_id is not None:
            cid = str(self.current_source_id)
            if '_' in cid:  # Combine mode (e.g., "8_7")
                return f"_S_{cid}"
            return f"_S{cid}"
        return ""

    # ======================================================================
    # Source 5 — FR3D: minimal wrapper around an EXTERNAL official
    # BGSU fr3d-python checkout. Nothing here installs, vendors, patches, or
    # checks out FR3D on the user's behalf. Registration is per-session via
    #     rmv_db FR3D   (auto-registers the bundled config)
    # and the search is argument-free via `rmv_load_motif`.
    # ======================================================================

    def _fr3d_mode_aliases(self) -> Dict[str, str]:
        """Map user-facing data_mode names to internal behavior modes."""
        return {
            'cache': 'cache',
            'run_from_scratch': 'cif_local',
            # New user-requested names
            'run_fr3d_pipeline': 'cif_local',
            'fr3d_local_data': 'fr3d_native_data',
            'rna3dhub_web_interactions': 'rna3dhub_interactions',
            # Backward-compatible legacy names
            'cif_local': 'cif_local',
            'fr3d_native_data': 'fr3d_native_data',
            'rna3dhub_interactions': 'rna3dhub_interactions',
        }

    def _fr3d_reset_state(self):
        """Reset Source-5 session state to an unregistered baseline."""
        self.fr3d_config_path = ''
        self.fr3d_python_path = ''          # fr3d-python repo root
        self.fr3d_python_exe = ''           # validated interpreter
        self.fr3d_commit = ''
        self.fr3d_branch = ''
        self.fr3d_dirty = None
        self.fr3d_data_mode = 'run_fr3d_pipeline'
        self.fr3d_cache_path = str((Path(__file__).parent.parent / 'external' / 'fr3d' / 'fr3d_cache').resolve())
        self.fr3d_query_path = ''
        self.fr3d_query_selection = 'all'
        self.fr3d_default_query = ''
        self.fr3d_allow_network = False
        self.fr3d_run_output_path = str((Path(self.fr3d_output_dir) / 'runs').resolve())
        self.fr3d_query_timeout_seconds = 180
        self.fr3d_interactions_path = ''
        self.fr3d_on_missing_capability = 'error'
        self.fr3d_registered = False

    def _fr3d_sanitize_query_name(self, name: str) -> str:
        """Filesystem-safe query name (also used as the ingested motif type)."""
        safe = re.sub(r'[^A-Za-z0-9]+', '_', str(name or '').strip()).strip('_')
        return safe or 'query'

    def _fr3d_validate_repo(self, repo_root: str) -> str:
        """Return '' if repo_root is an official fr3d-python checkout, else an error string."""
        result = validate_fr3d_checkout(repo_root)
        if result['ok']:
            return ''
        missing = result.get('missing', [])
        if missing:
            return "not an fr3d-python checkout (missing: %s)" % ", ".join(missing)
        return ''

    def _fr3d_python_ok(self, python_exe: str, repo_root: str) -> Tuple[bool, str]:
        """Verify an interpreter can import FR3D's deps and the repo. Returns (ok, detail)."""
        if not python_exe or not os.path.isfile(python_exe):
            return False, "interpreter not found: %s" % python_exe
        search_dir = os.path.join(repo_root, 'fr3d', 'search')
        probe = (
            "import sys\n"
            "sys.path.insert(0, r'''%s''')\n"
            "sys.path.insert(0, r'''%s''')\n"
            "import numpy, scipy, pdbx\n"
            "import fr3d\n"
            "print('FR3D_OK')\n"
        ) % (search_dir, repo_root)
        try:
            result = subprocess.run(
                [python_exe, '-c', probe],
                capture_output=True, text=True, timeout=60, check=False,
            )
        except Exception as e:
            return False, "%s: %s" % (type(e).__name__, e)
        if result.returncode == 0 and 'FR3D_OK' in (result.stdout or ''):
            return True, ''
        detail = (result.stderr or result.stdout or '').strip().splitlines()
        return False, (detail[-1] if detail else 'unknown import failure')

    def _resolve_py_launcher(self, launcher: str) -> str:
        """Resolve the Windows 'py -3' launcher to a concrete python.exe path."""
        try:
            result = subprocess.run(
                [launcher, '-3', '-c', 'import sys; print(sys.executable)'],
                capture_output=True, text=True, timeout=30, check=False,
            )
        except Exception:
            return ''
        if result.returncode == 0:
            lines = (result.stdout or '').strip().splitlines()
            if lines and os.path.isfile(lines[-1]):
                return lines[-1]
        return ''

    def _fr3d_resolve_python(self, explicit: str) -> List[str]:
        """Return existing interpreter candidates, best first, cross-platform.

        Order: explicit config python_path -> previously validated interpreter
        -> PyMOL's own interpreter -> PATH (python3/python and the Windows
        'py -3' launcher) -> active venv/conda -> platform well-known locations.
        Non-existent paths, and POSIX-style paths on Windows (which abspath would
        mangle into a bogus 'C:\\usr\\bin\\python3'), are dropped.
        """
        candidates: List[str] = [explicit, self.fr3d_python_exe, sys.executable]

        for name in ('python3', 'python'):
            found = shutil.which(name)
            if found:
                candidates.append(found)

        for env_var in ('VIRTUAL_ENV', 'CONDA_PREFIX'):
            root = os.environ.get(env_var, '').strip()
            if not root:
                continue
            if os.name == 'nt':
                candidates.append(os.path.join(root, 'python.exe'))
                candidates.append(os.path.join(root, 'Scripts', 'python.exe'))
            else:
                candidates.append(os.path.join(root, 'bin', 'python3'))
                candidates.append(os.path.join(root, 'bin', 'python'))

        if os.name == 'nt':
            launcher = shutil.which('py')
            if launcher:
                resolved = self._resolve_py_launcher(launcher)
                if resolved:
                    candidates.append(resolved)
            for base in (os.environ.get('LOCALAPPDATA', ''),
                         os.environ.get('PROGRAMFILES', ''),
                         os.environ.get('PROGRAMFILES(X86)', ''),
                         'C:\\'):
                if not base:
                    continue
                candidates.extend(glob.glob(os.path.join(base, 'Programs', 'Python', 'Python3*', 'python.exe')))
                candidates.extend(glob.glob(os.path.join(base, 'Python3*', 'python.exe')))
        else:
            candidates.extend([
                '/opt/homebrew/bin/python3',
                '/usr/local/bin/python3',
                '/usr/bin/python3',
            ])

        ordered, seen = [], set()
        for c in candidates:
            c = str(c or '').strip()
            if not c:
                continue
            if os.name == 'nt' and (c.startswith('/') or c.startswith('\\')):
                continue
            c = os.path.abspath(os.path.expanduser(c))
            if c in seen or not os.path.isfile(c):
                continue
            seen.add(c)
            ordered.append(c)
        return ordered

    def _fr3d_record_commit(self, repo_root: str):
        """Record git branch/commit/dirty for provenance (best-effort, never fatal)."""
        self.fr3d_commit = ''
        self.fr3d_branch = ''
        self.fr3d_dirty = None
        git = shutil.which('git')
        if not git:
            return

        def _git(args):
            try:
                r = subprocess.run(
                    [git, '-C', repo_root] + args,
                    capture_output=True, text=True, timeout=15, check=False,
                )
                if r.returncode == 0:
                    return (r.stdout or '').strip()
            except Exception:
                return ''
            return ''

        self.fr3d_commit = _git(['rev-parse', 'HEAD'])
        self.fr3d_branch = _git(['rev-parse', '--abbrev-ref', 'HEAD'])
        if self.fr3d_commit:
            self.fr3d_dirty = bool(_git(['status', '--porcelain']))

    def register_fr3d_source(self, config_path: str) -> bool:
        """Register an external official fr3d-python checkout for Source 5 from a JSON config.

        Validates the repository layout, an importing interpreter (numpy/scipy/
        mmcif-pdbx + the repo), the query path, and records git provenance.
        Never installs, checks out, vendors, or edits the user's FR3D copy.
        """
        raw = str(config_path or '').strip()
        if not raw:
            self.command_error("Usage: rmv_fr3d register /absolute/path/to/fr3d_config.json")
            return False
        cfg_abs = os.path.abspath(os.path.expanduser(raw))
        if not os.path.isfile(cfg_abs):
            self.logger.error(f"FR3D config file not found: {cfg_abs}")
            return False
        try:
            with open(cfg_abs, 'r', encoding='utf-8') as f:
                cfg = json.load(f)
        except Exception as e:
            self.logger.error(f"Failed to parse FR3D config JSON: {type(e).__name__}: {e}")
            return False
        if not isinstance(cfg, dict):
            self.logger.error("FR3D config must be a JSON object")
            return False

        cfg_dir = os.path.dirname(cfg_abs)

        def _abspath(value) -> str:
            value = str(value or '').strip()
            if not value:
                return ''
            p = os.path.expanduser(value)
            if not os.path.isabs(p):
                p = os.path.join(cfg_dir, p)
            return os.path.abspath(p)

        mode_aliases = self._fr3d_mode_aliases()
        data_mode = str(cfg.get('data_mode', 'run_from_scratch') or 'run_from_scratch').strip()
        if data_mode not in mode_aliases:
            self.logger.error(f"FR3D config: unknown data_mode '{data_mode}'")
            self.logger.info("  Allowed: cache | run_from_scratch")
            return False
        canonical_mode = mode_aliases[data_mode]

        if canonical_mode == 'cache':
            cache_path = _abspath(cfg.get('cache_path'))
            if not cache_path:
                cache_path = str((Path(cfg_dir).parent / 'external' / 'fr3d' / 'fr3d_cache').resolve())
            if not os.path.isdir(cache_path):
                self.logger.error(f"FR3D cache directory not found: {cache_path}")
                return False
            self.fr3d_config_path = cfg_abs
            self.fr3d_data_mode = data_mode
            self.fr3d_cache_path = cache_path
            self.fr3d_allow_network = False
            self.fr3d_query_path = ''
            self.fr3d_query_selection = 'cache'
            self.fr3d_registered = True
            self.user_data_paths[5] = cache_path
            self.logger.success(f"FR3D cache registered: {cache_path}")
            return True

        # --- required: fr3d_python_path ------------------------------------
        repo_root = _abspath(cfg.get('fr3d_python_path'))
        if not repo_root or not os.path.isdir(repo_root):
            self.logger.error("FR3D config: 'fr3d_python_path' must point to an existing fr3d-python checkout")
            return False
        repo_err = self._fr3d_validate_repo(repo_root)
        if repo_err:
            self.logger.error(f"FR3D config: {repo_err}")
            self.logger.info(f"  Checked under: {repo_root}")
            return False

        # --- query path ----------------------------------------------------
        query_path = _abspath(cfg.get('query_path'))
        if not query_path or not os.path.exists(query_path):
            self.logger.error("FR3D config: 'query_path' must point to an existing .json query file or directory")
            return False
        query_selection = str(cfg.get('query_selection', 'all') or 'all').strip().lower()
        default_query = str(cfg.get('default_query', '') or '').strip()
        if os.path.isdir(query_path):
            if query_selection not in ('all', 'default'):
                self.logger.error("FR3D config: 'query_selection' must be 'all' or 'default' for a query directory")
                return False
            if query_selection == 'default' and not default_query:
                self.logger.error("FR3D config: 'default_query' is required when query_selection is 'default'")
                return False
            if query_selection == 'default' and not default_query.lower().endswith('.json'):
                self.logger.error("FR3D config: 'default_query' must be an exact .json filename when query_selection is 'default'")
                return False
        elif not query_path.lower().endswith('.json'):
            self.logger.error("FR3D config: a single-file 'query_path' must be a .json query")
            return False

        # --- interpreter ---------------------------------------------------
        explicit_python = _abspath(cfg.get('python_path'))
        chosen_python, py_detail = '', ''
        explicit_fail_detail = ''
        for cand in self._fr3d_resolve_python(explicit_python):
            ok, detail = self._fr3d_python_ok(cand, repo_root)
            if ok:
                chosen_python = cand
                break
            py_detail = detail
            if explicit_python and cand == explicit_python:
                # Keep probing fallbacks (notably PyMOL's own interpreter) in case
                # the configured python_path has an architecture mismatch.
                explicit_fail_detail = detail
        if not chosen_python:
            self.logger.error("FR3D config: no usable interpreter can import FR3D and its deps (numpy, scipy, mmcif-pdbx)")
            if explicit_python:
                self.logger.info(f"  Configured python_path: {explicit_python}")
                if explicit_fail_detail:
                    self.logger.info(f"  Configured python_path error: {explicit_fail_detail}")
            if py_detail:
                self.logger.info(f"  Last import error: {py_detail}")
            self.logger.info("  Fix automatically with: rmv_setup FR3D")
            return False
        if explicit_python and chosen_python != explicit_python:
            self.logger.warning(
                f"FR3D config: configured python_path failed; using fallback interpreter: {chosen_python}"
            )
            if explicit_fail_detail:
                self.logger.info(f"  Configured python_path error: {explicit_fail_detail}")

        # --- optional paths + policy --------------------------------------
        allow_network = bool(cfg.get('allow_network', False))
        run_output_path = _abspath(cfg.get('run_output_path')) or str((Path(self.fr3d_output_dir) / 'runs').resolve())
        timeout_raw = cfg.get('query_timeout_seconds', 180)
        try:
            query_timeout_seconds = int(timeout_raw)
            if query_timeout_seconds < 10:
                raise ValueError
        except Exception:
            self.logger.error("FR3D config: 'query_timeout_seconds' must be an integer >= 10")
            return False
        interactions_path = _abspath(cfg.get('interactions_path'))
        if interactions_path and not os.path.exists(interactions_path):
            self.logger.error(f"FR3D config: interactions_path does not exist: {interactions_path}")
            return False
        on_missing = str(cfg.get('on_missing_capability', 'error') or 'error').strip().lower()
        if on_missing == 'fallback_to_cif_local':
            on_missing = 'fallback_to_run_fr3d_pipeline'
        if on_missing not in ('error', 'fallback_to_run_fr3d_pipeline', 'ask'):
            on_missing = 'error'
        if canonical_mode == 'rna3dhub_interactions' and not interactions_path and not allow_network:
            self.logger.error("FR3D config: data_mode 'rna3dhub_web_interactions' needs 'interactions_path' or allow_network=true")
            return False

        # --- provenance + commit state ------------------------------------
        self._fr3d_record_commit(repo_root)
        self.fr3d_config_path = cfg_abs
        self.fr3d_python_path = repo_root
        self.fr3d_python_exe = chosen_python
        self.fr3d_data_mode = data_mode
        self.fr3d_query_path = query_path
        self.fr3d_query_selection = query_selection
        self.fr3d_default_query = default_query
        self.fr3d_allow_network = allow_network
        self.fr3d_run_output_path = run_output_path
        self.fr3d_query_timeout_seconds = query_timeout_seconds
        self.fr3d_interactions_path = interactions_path
        self.fr3d_on_missing_capability = on_missing
        self.fr3d_registered = True
        try:
            os.makedirs(self.fr3d_run_output_path, exist_ok=True)
        except Exception:
            pass

        self.logger.success("FR3D registered (external official fr3d-python)")
        commit = self.fr3d_commit[:12] if self.fr3d_commit else "n/a"
        if os.path.isdir(query_path):
            extra = f", default={default_query}" if query_selection == 'default' else ''
            qdesc = f"{query_path} (selection={query_selection}{extra})"
        else:
            qdesc = query_path
        self.logger.info(f"  Repo    : {repo_root}  (commit {commit})")
        self.logger.info(f"  Python  : {chosen_python}")
        self.logger.info(
            f"  Mode    : {data_mode}  |  network {'on' if allow_network else 'off'}  |  "
            f"timeout {self.fr3d_query_timeout_seconds}s"
        )
        self.logger.info(f"  Queries : {qdesc}")
        if interactions_path:
            self.logger.info(f"  Interact: {interactions_path}")
        self.logger.info(f"  Output  : {self.fr3d_run_output_path}")
        return True

    def _fr3d_discover_queries(self) -> Tuple[List[str], str]:
        """Resolve the ordered list of query .json files to run.

        Directory rules: top-level ``.json`` only, sorted, no recursion, no
        guessing. Returns (queries, error).
        """
        qp = self.fr3d_query_path
        if not qp or not os.path.exists(qp):
            return [], "registered query_path no longer exists: %s" % qp
        if os.path.isfile(qp):
            if not qp.lower().endswith('.json'):
                return [], "query_path is not a .json file: %s" % qp
            return [qp], ''
        entries = sorted(
            os.path.join(qp, name)
            for name in os.listdir(qp)
            if name.lower().endswith('.json') and os.path.isfile(os.path.join(qp, name))
        )
        if not entries:
            return [], "no top-level .json query files found in: %s" % qp
        if self.fr3d_query_selection == 'default':
            target = self.fr3d_default_query
            for path in entries:
                base = os.path.basename(path)
                if base == target:
                    return [path], ''
            return [], "default_query '%s' not found among top-level .json files in %s" % (target, qp)
        return entries, ''

    def _fr3d_resolve_target_cif(self, pdb_id: str, run_dir: str) -> Tuple[str, str]:
        """Return an absolute mmCIF path for the target (cif_local). Returns (path, error).

        Uses the genuine mmCIF the user already has (an original local file, or
        PyMOL's fetched .cif) and never re-exports from PyMOL, to preserve the
        exact structure official FR3D would annotate. Only downloads from RCSB
        when allow_network is enabled.
        """
        pdb_upper = str(pdb_id or '').strip().upper()
        names = [pdb_upper.lower(), pdb_upper, str(getattr(self, 'loaded_pdb', '') or '').strip()]

        src = str(getattr(self, 'loaded_structure_path', '') or '').strip()
        if src and os.path.isfile(src) and src.lower().endswith(('.cif', '.mmcif')):
            return os.path.abspath(src), ''

        dirs = []
        try:
            fp = cmd.get('fetch_path')
            if fp:
                dirs.append(str(fp).strip())
        except Exception:
            pass
        dirs.append(os.getcwd())
        for d in dirs:
            if not d:
                continue
            for nm in names:
                if not nm:
                    continue
                cand = os.path.join(d, nm + '.cif')
                if os.path.isfile(cand):
                    return os.path.abspath(cand), ''

        cached_dir = Path(__file__).parent.parent / 'cached_structures'
        for nm in names:
            if nm:
                cached = cached_dir / f"{nm.lower()}.cif"
                if cached.is_file():
                    return str(cached.resolve()), ''

        if self.fr3d_allow_network and len(pdb_upper) == 4 and pdb_upper.isalnum():
            try:
                os.makedirs(run_dir, exist_ok=True)
                dest = os.path.join(run_dir, pdb_upper + '.cif')
                url = "https://files.rcsb.org/download/%s.cif" % pdb_upper
                urllib.request.urlretrieve(url, dest)
                if os.path.isfile(dest) and os.path.getsize(dest) > 0:
                    return dest, ''
                return '', "downloaded file was empty: %s" % url
            except Exception as e:
                return '', "failed to download %s.cif from RCSB: %s" % (pdb_upper, e)

        return '', (
            "could not locate a local mmCIF for %s (searched PyMOL fetch_path and CWD). "
            "Enable allow_network in the FR3D config to download it, or load the structure "
            "from a .cif file." % pdb_upper
        )

    def run_fr3d_search(self, pdb_id: str = None) -> bool:
        """Argument-free Source-5 entry point.

        Reproduces official FR3D for the current target and the registered
        query/queries using cif_local (the reference default). Never guesses the
        target or the queries and stops with a clear message when anything is
        ambiguous or missing.
        """
        if not self.fr3d_registered:
            self.logger.error("FR3D Source 5 is not registered.")
            self.logger.info("Register the external official fr3d-python first:")
            self.logger.info("  rmv_db FR3D")
            return False

        target = str(pdb_id or self.loaded_pdb_id or '').strip()
        if not target:
            self.logger.error("No structure loaded. Use rmv_fetch <PDB_ID> first.")
            return False
        target_upper = target.upper()

        # Capability gating: only run_fr3d_pipeline/cif_local is implemented as the verified
        # reference. Other modes are opt-in and must be explicitly allowed.
        mode_aliases = self._fr3d_mode_aliases()
        effective_mode = mode_aliases.get(self.fr3d_data_mode, self.fr3d_data_mode)
        if effective_mode != 'cif_local':
            if self.fr3d_on_missing_capability in ('fallback_to_run_fr3d_pipeline', 'fallback_to_cif_local'):
                self.logger.warning(f"data_mode '{self.fr3d_data_mode}' not available; falling back to run_fr3d_pipeline.")
                effective_mode = 'cif_local'
            else:
                self.logger.error(
                    f"data_mode '{self.fr3d_data_mode}' is not yet supported by this wrapper."
                )
                self.logger.info("Set data_mode to 'run_fr3d_pipeline', or on_missing_capability to 'fallback_to_run_fr3d_pipeline'.")
                return False

        queries, qerr = self._fr3d_discover_queries()
        if qerr:
            self.logger.error(f"FR3D query resolution failed: {qerr}")
            return False

        # Validate queries; skip malformed JSON files with a concise warning
        # instead of aborting the whole run so valid queries still execute.
        valid_queries = []
        invalid = []
        for q in queries:
            try:
                with open(q, 'r', encoding='utf-8') as f:
                    json.load(f)
                valid_queries.append(q)
            except Exception as e:
                invalid.append((q, str(e)))
        if invalid:
            self.logger.warning(f"Skipping {len(invalid)} query file(s) with invalid JSON:")
            for q, e in invalid:
                self.logger.info(f"  - {os.path.basename(q)}: {e}")
        if not valid_queries:
            self.logger.error("FR3D aborted: no valid query JSON files remain.")
            return False
        queries = valid_queries

        run_id = time.strftime('%Y%m%d_%H%M%S')
        run_root = os.path.join(self.fr3d_run_output_path, target_upper, run_id)
        os.makedirs(run_root, exist_ok=True)

        target_cif, cif_err = self._fr3d_resolve_target_cif(target_upper, run_root)
        if cif_err:
            self.logger.error(f"FR3D target preparation failed: {cif_err}")
            return False

        runner = Path(__file__).parent / 'tools' / 'fr3d_search_runner.py'
        if not runner.is_file():
            self.logger.error(f"Missing FR3D runner script: {runner}")
            return False

        ingest_dir = os.path.join(run_root, 'ingest')
        # Persistent FR3D annotation cache shared across ALL runs. FR3D names its
        # cache files by structure (e.g. 4V9F-1-0_NA.pickle), so reusing one
        # directory is safe and matches FR3D's native DATAPATH design. This is
        # what stops PyMOL freezing: the expensive one-time annotation of the
        # large reference structures the queries embed (4V9F, 7K00, 8GLP, ...)
        # is built once and reused, instead of being rebuilt on every run.
        shared_data_dir = os.path.join(self.fr3d_run_output_path, '_fr3d_cache')
        os.makedirs(ingest_dir, exist_ok=True)
        os.makedirs(shared_data_dir, exist_ok=True)

        cache_primed = os.path.isdir(os.path.join(shared_data_dir, 'units')) and any(
            name.endswith('_NA.pickle')
            for name in os.listdir(os.path.join(shared_data_dir, 'units'))
        )
        self.logger.info(
            f"Running official FR3D on {target_upper} for {len(queries)} query(ies) [{self.fr3d_data_mode}]..."
        )
        if not cache_primed:
            self.logger.info(
                "  First FR3D run builds a one-time annotation cache for the reference "
                "structures embedded in the queries (this can take several minutes; "
                "subsequent runs reuse it and are fast)."
            )

        manifest = {
            'target': target_upper,
            'target_cif': target_cif,
            'run_id': run_id,
            'data_mode': self.fr3d_data_mode,
            'data_mode_internal': effective_mode,
            'fr3d_repo': self.fr3d_python_path,
            'fr3d_commit': self.fr3d_commit,
            'fr3d_branch': self.fr3d_branch,
            'fr3d_dirty': self.fr3d_dirty,
            'python': self.fr3d_python_exe,
            'queries': [],
        }

        produced = 0
        for idx, query_file in enumerate(queries, 1):
            qname = self._fr3d_sanitize_query_name(os.path.splitext(os.path.basename(query_file))[0])
            q_run_dir = os.path.join(run_root, 'queries', qname)
            os.makedirs(q_run_dir, exist_ok=True)
            self.logger.info(
                f"  [{idx}/{len(queries)}] Running FR3D query '{qname}' (timeout={self.fr3d_query_timeout_seconds}s)..."
            )
            spec = {
                'fr3d_root': self.fr3d_python_path,
                'query_file': query_file,
                'target_cif': target_cif,
                'run_dir': q_run_dir,
                'data_dir': shared_data_dir,
                'query_name': qname,
                'allow_network': self.fr3d_allow_network,
                'data_mode': self.fr3d_data_mode,
            }
            spec_path = os.path.join(q_run_dir, 'runner_spec.json')
            with open(spec_path, 'w', encoding='utf-8') as f:
                json.dump(spec, f, indent=2)

            try:
                proc = subprocess.run(
                    [self.fr3d_python_exe, str(runner), '--spec', spec_path],
                    capture_output=True, text=True, check=False,
                    timeout=self.fr3d_query_timeout_seconds,
                )
            except subprocess.TimeoutExpired:
                err = (
                    f"timed out after {self.fr3d_query_timeout_seconds}s; "
                    "try query_selection='default' or increase query_timeout_seconds"
                )
                self.logger.error(f"  [{qname}] FR3D search failed: {err}")
                manifest['queries'].append(
                    {
                        'query_name': qname,
                        'query_file': query_file,
                        'ok': False,
                        'timed_out': True,
                        'error': err,
                    }
                )
                continue
            except Exception as e:
                self.logger.error(f"  [{qname}] failed to launch FR3D: {type(e).__name__}: {e}")
                manifest['queries'].append(
                    {'query_name': qname, 'query_file': query_file, 'ok': False, 'error': str(e)}
                )
                continue

            report = {}
            for line in reversed((proc.stdout or '').strip().splitlines()):
                line = line.strip()
                if line.startswith('{'):
                    try:
                        report = json.loads(line)
                        break
                    except Exception:
                        continue

            if proc.returncode != 0 or not report.get('ok'):
                err = report.get('error')
                if not err:
                    tail = (proc.stderr or '').strip().splitlines()
                    err = tail[-1] if tail else 'unknown error'
                self.logger.error(f"  [{qname}] FR3D search failed: {err}")
                manifest['queries'].append(
                    {'query_name': qname, 'query_file': query_file, 'ok': False, 'error': str(err)}
                )
                continue

            result = report.get('result', {}) or {}
            csv_src = str(result.get('csv_path', '') or '')
            count = int(result.get('candidate_count', 0) or 0)
            entry = {
                'query_name': qname,
                'query_file': query_file,
                'ok': True,
                'candidate_count': count,
                'changed_fields': result.get('changed_fields', {}),
                'run_dir': result.get('run_dir', q_run_dir),
            }
            if csv_src and os.path.isfile(csv_src):
                dest = os.path.join(ingest_dir, f"{target_upper}_fr3d_query_{qname}.csv")
                try:
                    shutil.copyfile(csv_src, dest)
                    entry['csv'] = dest
                    produced += 1
                    self.logger.success(f"  [{qname}] {count} candidate(s)")
                except Exception as e:
                    entry['ok'] = False
                    entry['error'] = f"failed to ingest CSV: {e}"
                    self.logger.error(f"  [{qname}] {entry['error']}")
            else:
                self.logger.warning(f"  [{qname}] no candidate CSV produced")
            manifest['queries'].append(entry)

        try:
            with open(os.path.join(run_root, 'manifest.json'), 'w', encoding='utf-8') as f:
                json.dump(manifest, f, indent=2)
        except Exception:
            pass

        if produced == 0:
            self.logger.warning("FR3D produced no candidates for any query; nothing to load.")
            self.logger.info(f"Provenance: {run_root}")
            return False

        self.user_data_paths[5] = ingest_dir
        self._handle_source_by_id(5, ingest_dir)
        self.load_user_annotations_action('fr3d', target_upper, auto_pipeline=False)

        loaded = self.viz_manager.motif_loader.get_loaded_motifs()
        if loaded:
            total = sum(len(info.get('motif_details', [])) for info in loaded.values())
            self.logger.success(f"Loaded FR3D motifs: {len(loaded)} type(s), {total} instance(s)")
            self.logger.info(f"Provenance: {os.path.join(run_root, 'manifest.json')}")
            return True
        self.logger.warning("FR3D output generated, but no motifs were loaded into RSMViewer.")
        return False

    def print_fr3d_status(self):
        """Print Source-5 FR3D registration status and usage."""
        print("\n" + "=" * 70)
        print("FR3D Source 5 - external official fr3d-python wrapper")
        print("=" * 70)
        if not self.fr3d_registered:
            print("Status     : NOT registered")
            print("Register   : rmv_db FR3D   (or rmv_fr3d register <config>)")
            print("=" * 70 + "\n")
            return
        print("Status     : registered")
        print(f"Config     : {self.fr3d_config_path}")
        print(f"Repo       : {self.fr3d_python_path}")
        if self.fr3d_commit:
            state = self.fr3d_branch or 'detached'
            state += '  dirty' if self.fr3d_dirty else '  clean'
            print(f"Commit     : {self.fr3d_commit[:12]}  {state}")
        else:
            print("Commit     : (git metadata unavailable)")
        print(f"Python     : {self.fr3d_python_exe}")
        print(f"Data mode  : {self.fr3d_data_mode}")
        if os.path.isdir(self.fr3d_query_path):
            print(f"Queries    : {self.fr3d_query_path}  (selection={self.fr3d_query_selection})")
            if self.fr3d_query_selection == 'default':
                print(f"Default    : {self.fr3d_default_query}")
        else:
            print(f"Query      : {self.fr3d_query_path}")
        print(f"Network    : {'allowed' if self.fr3d_allow_network else 'disabled'}")
        print(f"Timeout    : {self.fr3d_query_timeout_seconds}s per query")
        if self.fr3d_interactions_path:
            print(f"Interact.  : {self.fr3d_interactions_path}")
        print(f"Output     : {self.fr3d_run_output_path}")
        print("\nRun        : rmv_fetch <PDB_ID> ; rmv_load_motif")
        print("=" * 70 + "\n")

    def _fr3d_default_repo_root(self) -> Tuple[str, str]:
        """Resolve and validate the fr3d-python checkout from the default config.

        Returns (repo_root, error). ``error`` is empty when the checkout exists
        and looks like an official fr3d-python repo.
        """
        cfg_abs = os.path.abspath(os.path.expanduser(str(self.default_fr3d_config or '')))
        if not os.path.isfile(cfg_abs):
            return '', f"default config not found: {cfg_abs}"
        try:
            with open(cfg_abs, 'r', encoding='utf-8') as f:
                cfg = json.load(f)
        except Exception as e:
            return '', f"could not parse {cfg_abs}: {type(e).__name__}: {e}"
        cfg_dir = os.path.dirname(cfg_abs)
        raw = str(cfg.get('fr3d_python_path', '') or '').strip()
        if not raw:
            return '', "config is missing 'fr3d_python_path'"
        repo = os.path.expanduser(raw)
        if not os.path.isabs(repo):
            repo = os.path.join(cfg_dir, repo)
        repo = os.path.abspath(repo)
        if not os.path.isdir(repo):
            return repo, "fr3d-python checkout folder not found"
        repo_err = self._fr3d_validate_repo(repo)
        if repo_err:
            return repo, repo_err
        return repo, ''

    def auto_setup_fr3d(self, python_exe: str = '') -> bool:
        """One-shot setup for FR3D (Source 5): the `rmv_setup FR3D` entry point.

        The user only pastes the official fr3d-python checkout; this discovers a
        working interpreter, installs numpy/scipy/mmcif-pdbx if needed, and
        registers Source 5. Every failure is reported with the precise reason.
        """
        repo_root, repo_err = self._fr3d_default_repo_root()
        if repo_err:
            self.logger.error(f"FR3D setup: {repo_err}")
            self.logger.info("  Paste the official BGSU fr3d-python checkout into:")
            self.logger.info(f"    {repo_root or 'external/fr3d/fr3d-python-latest'}")
            self.logger.info("  It must contain fr3d/__init__.py and fr3d/search/FR3D.py, then re-run: rmv_setup FR3D")
            return False

        self.logger.info("Setting up FR3D (Source 5)...")
        self.logger.info(f"  Checkout: {repo_root}")

        # Fast path: an interpreter already imports the repo + all deps.
        for cand in self._fr3d_resolve_python(python_exe):
            ok, _ = self._fr3d_python_ok(cand, repo_root)
            if ok:
                self.logger.success(f"Ready interpreter with FR3D dependencies: {cand}")
                if self.register_fr3d_source(str(self.default_fr3d_config)):
                    self.logger.success("FR3D is ready.")
                    return True
                return False

        # No ready interpreter: install the dependencies, then register.
        self.logger.info("No interpreter has FR3D's dependencies yet; installing numpy/scipy/mmcif-pdbx...")
        if not self.setup_fr3d_environment(python_exe):
            self.logger.error("FR3D setup could not install the required dependencies.")
            self.logger.info("  Common causes: no internet on the first install, or a read-only interpreter.")
            self.logger.info("  Fix: point setup at a Python that has pip + internet, e.g.:")
            self.logger.info("    rmv_setup FR3D /absolute/path/to/python")
            return False

        if self.register_fr3d_source(str(self.default_fr3d_config)):
            self.logger.success("FR3D is fully set up.")
            return True
        self.logger.error("Dependencies installed but FR3D registration still failed; see the reason above.")
        return False

    def setup_fr3d_environment(self, python_exe: str = '') -> bool:
        """Install FR3D's Python dependencies so Source 5 is ready to run.

        Installs numpy, scipy and mmcif-pdbx into the interpreter that will run
        FR3D (an explicitly given one, the registered interpreter, or the
        current Python), then verifies the imports. This NEVER touches the
        user's fr3d-python checkout -- it only prepares an interpreter.
        """
        requirements = ['numpy', 'scipy', 'mmcif-pdbx']   # import names: numpy, scipy, pdbx

        target = str(python_exe or '').strip()
        if not target:
            target = str(getattr(self, 'fr3d_python_exe', '') or '').strip()
        if not target:
            target = sys.executable
        target = os.path.abspath(os.path.expanduser(target))
        if not os.path.isfile(target):
            self.logger.error(f"Python interpreter not found: {target}")
            self.logger.info("Usage: rmv_setup FR3D [/abs/path/to/python]")
            return False

        self.logger.info(f"Installing FR3D requirements into: {target}")
        self.logger.info(f"  Packages: {', '.join(requirements)}")

        # Ensure pip exists (bootstrap via ensurepip if needed).
        try:
            pip_probe = subprocess.run(
                [target, '-m', 'pip', '--version'],
                capture_output=True, text=True, timeout=60, check=False,
            )
        except Exception as e:
            self.logger.error(f"Could not run pip with {target}: {type(e).__name__}: {e}")
            return False
        if pip_probe.returncode != 0:
            self.logger.info("pip not available; attempting to bootstrap with ensurepip...")
            try:
                subprocess.run(
                    [target, '-m', 'ensurepip', '--upgrade'],
                    capture_output=True, text=True, timeout=300, check=False,
                )
            except Exception:
                pass

        def _run_install(extra_args):
            cmd = [target, '-m', 'pip', 'install', '--upgrade'] + extra_args + requirements
            self.logger.info(f"  $ {' '.join(cmd)}")
            try:
                return subprocess.run(cmd, capture_output=True, text=True, check=False)
            except Exception as e:
                self.logger.error(f"pip install failed to launch: {type(e).__name__}: {e}")
                return None

        proc = _run_install([])
        if proc is None:
            return False
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or '').strip().splitlines()
            for line in tail[-8:]:
                self.logger.info(f"  {line}")
            self.logger.info("Retrying installation with --user ...")
            proc = _run_install(['--user'])
            if proc is None or proc.returncode != 0:
                if proc is not None:
                    tail = (proc.stderr or proc.stdout or '').strip().splitlines()
                    for line in tail[-8:]:
                        self.logger.info(f"  {line}")
                self.logger.error("FR3D dependency installation failed.")
                return False

        # Verify the interpreter can now import the FR3D dependencies.
        try:
            verify = subprocess.run(
                [target, '-c', "import numpy, scipy, pdbx; print('FR3D_DEPS_OK')"],
                capture_output=True, text=True, timeout=120, check=False,
            )
        except Exception as e:
            self.logger.error(f"Verification failed: {type(e).__name__}: {e}")
            return False
        if verify.returncode != 0 or 'FR3D_DEPS_OK' not in (verify.stdout or ''):
            self.logger.error("FR3D dependencies still not importable after installation.")
            detail = (verify.stderr or verify.stdout or '').strip().splitlines()
            if detail:
                self.logger.info(f"  {detail[-1]}")
            return False

        self.logger.success("FR3D requirements installed and verified.")

        # If a repo is registered, adopt this interpreter as the FR3D python
        # when it can fully import the checkout.
        if getattr(self, 'fr3d_registered', False) and self.fr3d_python_path:
            ok, detail = self._fr3d_python_ok(target, self.fr3d_python_path)
            if ok:
                self.fr3d_python_exe = target
                self.logger.info(f"Source 5 interpreter set to: {target}")
            elif detail:
                self.logger.info(f"Note: interpreter installed deps but repo import check said: {detail}")

        if getattr(self, 'fr3d_registered', False):
            self.logger.info("FR3D is ready.")
        else:
            self.logger.info("FR3D interpreter is ready.")
        return True

    def _get_current_source_motifs(self):
        """Return a *filtered copy* of loaded_motifs showing only instances
        belonging to the current PDB + source.

        Keys whose filtered instance list is empty are omitted.
        Counts are recalculated to match the filtered list.

        This allows rmv_summary / rmv_show to display only the active source
        while keeping all accumulated data available for cross-source
        superimposition.
        """
        loaded = self.viz_manager.motif_loader.get_loaded_motifs()
        if not loaded:
            return {}

        pdb_id = self.loaded_pdb_id
        suffix = self._get_source_suffix()

        filtered = {}
        for key, info in loaded.items():
            details = info.get('motif_details', [])
            # Keep instances matching current PDB + source
            kept = [
                d for d in details
                if d.get('_pdb_id', info.get('pdb_id', '')) == pdb_id
                and d.get('_source_suffix', info.get('source_suffix', '')) == suffix
            ]
            if kept:
                # Shallow copy the info dict, replace details & count
                new_info = dict(info)
                new_info['motif_details'] = kept
                new_info['count'] = len(kept)
                filtered[key] = new_info
        return filtered
    
    def _build_auth_label_chain_mapping(self, pdb_id):
        """Parse CIF file to build auth_asym_id -> label_asym_id chain mapping.
        
        When cif_use_auth=0, PyMOL uses label_asym_id as 'chain' and loses
        auth_asym_id from the model. But all motif databases reference auth chains.
        This method parses the CIF file from disk to recover the mapping.
        
        Returns:
            dict: {auth_asym_id: label_asym_id} for each unique chain pair
        """
        import os
        
        # Find CIF file in PyMOL's fetch path
        fetch_path = "."
        try:
            fp = cmd.get("fetch_path")
            if fp:
                fetch_path = str(fp).strip()
        except:
            pass
        
        cif_path = None
        for name in [f"{pdb_id.lower()}.cif", f"{pdb_id.upper()}.cif"]:
            candidate = os.path.join(fetch_path, name)
            if os.path.exists(candidate):
                cif_path = candidate
                break
        
        if not cif_path:
            self.logger.debug(f"CIF file not found in: {fetch_path}")
            return {}
        
        # Parse CIF _atom_site loop to extract auth_asym_id and label_asym_id columns
        auth_col = -1
        label_col = -1
        col_count = 0
        reading_headers = False
        reading_data = False
        mapping = {}
        
        try:
            with open(cif_path, 'r', encoding='utf-8') as f:
                for line in f:
                    stripped = line.strip()
                    
                    if not stripped or stripped.startswith('#'):
                        if reading_data:
                            break
                        continue
                    
                    if stripped.startswith('loop_'):
                        if reading_data:
                            break
                        reading_headers = False
                        reading_data = False
                        col_count = 0
                        auth_col = -1
                        label_col = -1
                        continue
                    
                    if stripped.startswith('_atom_site.'):
                        reading_headers = True
                        col_name = stripped.split()[0]
                        if col_name == '_atom_site.auth_asym_id':
                            auth_col = col_count
                        elif col_name == '_atom_site.label_asym_id':
                            label_col = col_count
                        col_count += 1
                        continue
                    
                    if reading_headers and not stripped.startswith('_'):
                        reading_headers = False
                        reading_data = True
                        if auth_col < 0 or label_col < 0:
                            self.logger.debug(f"CIF missing columns: auth_col={auth_col}, label_col={label_col}")
                            return {}
                    
                    if reading_data:
                        if stripped.startswith('_') or stripped.startswith('data_') or stripped.startswith('loop_'):
                            break
                        
                        tokens = stripped.split()
                        max_col = max(auth_col, label_col)
                        if len(tokens) > max_col:
                            auth_id = tokens[auth_col]
                            label_id = tokens[label_col]
                            # Keep first mapping per auth chain (ATOM records come before HETATM,
                            # so polymer chains are mapped first - correct for motif data)
                            if auth_id not in mapping:
                                mapping[auth_id] = label_id
            
            if mapping:
                self.logger.debug(f"CIF auth->label mapping ({len(mapping)} chains): "
                                 f"{dict(list(mapping.items())[:8])}")
            return mapping
            
        except Exception as e:
            self.logger.debug(f"Error parsing CIF for chain mapping: {e}")
            return {}

    def _load_rmsx_wrapper_config(self):
        """Load persisted RNAMotifScanX wrapper settings if available."""
        try:
            cfg_path = self.rmsx_config_file
            if not cfg_path.exists():
                return

            with open(cfg_path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            self.rmsx_executable_path = str(data.get('rmsx_executable_path', '') or '').strip()
            self.rmsx_working_dir = str(data.get('rmsx_working_dir', '') or '').strip()
            self.rmsx_args_template = str(data.get('rmsx_args_template', '') or '').strip()
            self.rmsx_auto_run_on_fetch = bool(data.get('rmsx_auto_run_on_fetch', False))
            self.rmsx_query_file = str(data.get('rmsx_query_file', '') or '').strip()
            output_path = str(data.get('rmsx_output_path', '') or '').strip()
            if output_path:
                self.rmsx_output_path = output_path
        except Exception as e:
            self.logger.warning(f"Could not read RNAMotifScanX wrapper config: {e}")

    def _save_rmsx_wrapper_config(self):
        """Persist RNAMotifScanX wrapper settings."""
        try:
            self.rmsx_config_file.parent.mkdir(parents=True, exist_ok=True)
            payload = {
                'rmsx_executable_path': self.rmsx_executable_path,
                'rmsx_working_dir': self.rmsx_working_dir,
                'rmsx_output_path': self.rmsx_output_path,
                'rmsx_args_template': self.rmsx_args_template,
                'rmsx_auto_run_on_fetch': self.rmsx_auto_run_on_fetch,
                'rmsx_query_file': getattr(self, 'rmsx_query_file', ''),
            }
            with open(self.rmsx_config_file, 'w', encoding='utf-8') as f:
                json.dump(payload, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Could not save RNAMotifScanX wrapper config: {e}")

    def _is_executable_runnable(self, executable_path: str) -> bool:
        """Return True if an executable path exists and can be invoked."""
        if not executable_path:
            return False

        expanded = os.path.abspath(os.path.expanduser(executable_path))
        if not os.path.isfile(expanded):
            return False

        if not os.access(expanded, os.X_OK):
            return False

        try:
            # We only care whether the process can start, not whether it exits 0.
            subprocess.run(
                [expanded],
                capture_output=True,
                text=True,
                check=False,
                timeout=5,
            )
            return True
        except subprocess.TimeoutExpired:
            return True
        except Exception:
            return False

    def print_rmsx_wrapper_status(self):
        """Print current RNAMotifScanX wrapper settings and quick usage hints."""
        print("\n" + "=" * 70)
        print("RNAMotifScanX Wrapper Status")
        print("=" * 70)
        print(f"Executable path : {self.rmsx_executable_path or '(not set)'}")
        print(f"Working dir     : {self.rmsx_working_dir or '(directory of executable)'}")
        print(f"Output path     : {self.rmsx_output_path}")
        print(f"Args template   : {self.rmsx_args_template or '(none)'}")
        print(f"Auto on fetch   : {'on' if self.rmsx_auto_run_on_fetch else 'off'}")
        print(f"Query file      : {getattr(self, 'rmsx_query_file', '') or '(built-in consensus set)'}")
        print(f"Runtime dir     : {self.rmsx_runtime_dir}")
        print("\nCommands:")
        print("  rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH] [QUERY_FILE]")
        print("  rmv_rmsx args <ARG_TEMPLATE>")
        print("  rmv_rmsx doctor")
        print("  rmv_rmsx setup")
        print("  rmv_rmsx test")
        print("  rmv_rmsx run <PDB_ID> [EXTRA_ARGS]   # Full fresh rerun")
        print("  rmv_rmsx run_current [EXTRA_ARGS]")
        print("  rmv_rmsx scan_prepared <PDB_ID> [CHAINS] [compare]  # Run scan on prepared inputs")
        print("  rmv_rmsx scan_cancel                 # Cancel an in-progress scan_prepared run")
        print("\nTemplate placeholders:")
        print("  {pdb_id} {pdb_lower} {output_dir} {work_dir}")
        print("Example:")
        print("  rmv_rmsx args --pdb {pdb_id} --out {output_dir}")
        print("=" * 70 + "\n")

    def configure_rmsx_wrapper(self, executable: str, output_dir: str = '', work_dir: str = '', auto_on_fetch: str = '', query_file: str = ''):
        """Set RNAMotifScanX wrapper settings and persist them."""
        if not executable:
            self.logger.error("RNAMotifScanX executable path is required")
            self.logger.info("Usage: rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH] [QUERY_FILE]")
            return False

        expanded_exe = os.path.abspath(os.path.expanduser(executable))
        if not os.path.isfile(expanded_exe):
            self.logger.error(f"RNAMotifScanX executable not found: {expanded_exe}")
            return False
        if not os.access(expanded_exe, os.X_OK):
            self.logger.error(f"RNAMotifScanX executable is not runnable: {expanded_exe}")
            return False

        self.rmsx_executable_path = expanded_exe

        if output_dir:
            self.rmsx_output_path = os.path.abspath(os.path.expanduser(output_dir))

        if work_dir:
            expanded_work = os.path.abspath(os.path.expanduser(work_dir))
            if not os.path.isdir(expanded_work):
                self.logger.error(f"Working directory not found: {expanded_work}")
                return False
            self.rmsx_working_dir = expanded_work

        if auto_on_fetch:
            auto_val = auto_on_fetch.strip().lower()
            if auto_val in ('1', 'on', 'true', 'yes'):
                self.rmsx_auto_run_on_fetch = True
            elif auto_val in ('0', 'off', 'false', 'no'):
                self.rmsx_auto_run_on_fetch = False
            else:
                self.logger.error("AUTO_ON_FETCH must be one of: on/off, true/false, 1/0")
                return False

        if query_file:
            expanded_query = os.path.abspath(os.path.expanduser(query_file))
            if not os.path.isfile(expanded_query):
                self.logger.error(f"RNAMotifScanX query file not found: {expanded_query}")
                return False
            self.rmsx_query_file = expanded_query

        try:
            os.makedirs(self.rmsx_output_path, exist_ok=True)
        except Exception as e:
            self.logger.error(f"Could not create RNAMotifScanX output directory {self.rmsx_output_path}: {type(e).__name__}: {e}")
            return False

        if not self._is_executable_runnable(self.rmsx_executable_path):
            self.logger.error(f"Configured RNAMotifScanX executable failed test run: {self.rmsx_executable_path}")
            return False

        self._save_rmsx_wrapper_config()
        self.logger.success("RNAMotifScanX wrapper configuration updated")
        self.logger.info(f"Executable path: {self.rmsx_executable_path}")
        self.logger.info(f"Working dir: {self.rmsx_working_dir or os.path.dirname(self.rmsx_executable_path)}")
        self.logger.info(f"Output path: {self.rmsx_output_path}")
        self.logger.info(f"Auto run on rmv_fetch: {'on' if self.rmsx_auto_run_on_fetch else 'off'}")
        self.logger.info(f"Query file: {getattr(self, 'rmsx_query_file', '') or '(built-in consensus set)'}")
        return True

    def set_rmsx_args_template(self, arg_template: str):
        """Set argument template used by rmv_rmsx run."""
        self.rmsx_args_template = str(arg_template or '').strip()
        self._save_rmsx_wrapper_config()
        self.logger.success("RNAMotifScanX argument template updated")
        self.logger.info(f"Template: {self.rmsx_args_template or '(none)'}")
        return True

    def test_rmsx_wrapper(self):
        """Run executable sanity checks for RNAMotifScanX wrapper."""
        if not self.rmsx_executable_path:
            self.logger.error("RNAMotifScanX executable is not configured")
            self.logger.info("Use: rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH] [QUERY_FILE]")
            return False

        if not self._is_executable_runnable(self.rmsx_executable_path):
            self.logger.error(f"RNAMotifScanX executable failed test run: {self.rmsx_executable_path}")
            return False

        self.logger.success("RNAMotifScanX executable is available and runnable")
        return True

    def run_rmsx_wrapper(self, pdb_id: str, extra_args: str = '', force_fresh: bool = False):
        """Run configured RNAMotifScanX executable and load results from source 7."""
        pdb_upper = str(pdb_id).strip().upper()
        if not pdb_upper:
            self.logger.error("PDB ID is required")
            self.logger.info("Usage: rmv_rmsx run <PDB_ID> [EXTRA_ARGS]")
            return False

        # data_mode=scan_prepared: run scan directly on prepared inputs.
        rmsx_cfg_probe = getattr(self, 'rmsx_pipeline_config', {}) or self._build_internal_rmsx_config()
        self.rmsx_pipeline_config = dict(rmsx_cfg_probe)
        if str(rmsx_cfg_probe.get('data_mode', '')).strip().lower() == 'scan_prepared':
            return self.run_rmsx_scan_prepared(pdb_upper, chains=str(extra_args or '').strip())

        query_override = ''
        if extra_args:
            try:
                import shlex
                tokens = shlex.split(extra_args)
            except Exception:
                tokens = [token for token in extra_args.split() if token]
            for token in tokens:
                lowered = token.lower()
                if lowered.endswith('.struct') or lowered.endswith('.txt') or lowered.startswith('query='):
                    query_override = token.split('=', 1)[1].strip() if '=' in token else token
                    break

        previous_query_file = getattr(self, 'rmsx_query_file', '')
        if query_override:
            self.rmsx_query_file = os.path.abspath(os.path.expanduser(query_override))

        if not self.ensure_rmsx_runtime_ready(auto_setup=True):
            if query_override:
                self.rmsx_query_file = previous_query_file
            return False

        # Preferred path: run integrated source-7 RMSX runtime config
        rmsx_cfg = getattr(self, 'rmsx_pipeline_config', {}) or self._build_internal_rmsx_config()
        self.rmsx_pipeline_config = dict(rmsx_cfg)
        if rmsx_cfg:
            try:
                data_mode = str(rmsx_cfg.get('data_mode', 'preannotated')).strip().lower()
                if data_mode in ('preannotated', 'cache', 'cached') and not force_fresh:
                    copied = self._copy_preannotated_rmsx(rmsx_cfg, pdb_upper)
                    if copied.get('copied', 0):
                        self.logger.info(
                            f"Loaded {copied['copied']} preannotated RNAMotifScanX result(s) for {pdb_upper}."
                        )
                        self.load_user_annotations_action('rnamotifscanx', pdb_upper, auto_pipeline=False)
                        return True
                    self.logger.warning(
                        f"No preannotated RNAMotifScanX results found for {pdb_upper}; "
                        "set data_mode to run_from_scratch in config/rmsx_config.json to execute RMSX."
                    )
                    return False

                tools_dir = str(Path(__file__).parent / 'tools')
                if tools_dir not in sys.path:
                    sys.path.insert(0, tools_dir)

                from rmsx_runner import run_pipeline as rmsx_run  # type: ignore

                # Ensure output_dir is always defined for predictable loading.
                if not str(rmsx_cfg.get('output_dir', '') or '').strip():
                    rmsx_cfg = dict(rmsx_cfg)
                    rmsx_cfg['output_dir'] = self.rmsx_output_path

                current_query_file = str(getattr(self, 'rmsx_query_file', '') or '').strip()
                if current_query_file:
                    rmsx_cfg = dict(rmsx_cfg)
                    rmsx_cfg['query_file'] = current_query_file

                if not query_override:
                    query_override = str(rmsx_cfg.get('query_file', '') or '').strip()

                if query_override:
                    rmsx_cfg = dict(rmsx_cfg)
                    rmsx_cfg['query_file'] = query_override
                    self.logger.info(f"Using RNAMotifScanX query file: {query_override}")

                # Source 7 strict mode: never auto-download CIF.
                rmsx_cfg = dict(rmsx_cfg)
                rmsx_cfg['auto_download_cif'] = False

                out_dir = os.path.abspath(os.path.expanduser(str(rmsx_cfg.get('output_dir', self.rmsx_output_path))))
                rmsx_cfg['cif_input_dir'] = out_dir
                local_pdb = self._prepare_local_pdb_for_rmsx(pdb_upper, out_dir, force_refresh=force_fresh)
                if local_pdb:
                    rmsx_cfg['auto_download_pdb'] = False

                mode_text = 'fresh' if force_fresh else 'incremental/prebuilt-aware'
                self.logger.info(f"Running RNAMotifScanX pipeline for {pdb_upper} ({mode_text})...")
                if extra_args and not query_override:
                    self.logger.info("Note: unrecognized extra args are ignored when using pipeline config mode.")

                results = rmsx_run(rmsx_cfg, pdb_upper, force_fresh=force_fresh)
                if not results:
                    out_dir = os.path.abspath(os.path.expanduser(str(rmsx_cfg.get('output_dir', self.rmsx_output_path))))
                    if sys.platform == 'darwin':
                        self.logger.warning(f"No RNAMotifScanX result files found for {pdb_upper} in {out_dir}.")
                        self.logger.warning(
                            "RNAMotifScanX execution is Linux x86-64 only. "
                            "Generate results on Linux, then copy result_0_100_withbs.log files to this output directory."
                        )
                    else:
                        self.logger.error(
                            "RNAMotifScanX pipeline produced no result files. "
                            "Check executable/query/database paths in the RMSX config."
                        )
                    return False

                out_dir = os.path.abspath(os.path.expanduser(str(rmsx_cfg.get('output_dir', self.rmsx_output_path))))
                self.rmsx_output_path = out_dir
                self.user_data_paths[7] = out_dir

                self._handle_source_by_id(7, out_dir)
                self.load_user_annotations_action('rnamotifscanx', pdb_upper, auto_pipeline=False)

                loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs()
                if loaded_motifs:
                    total_instances = sum(len(info.get('motif_details', [])) for info in loaded_motifs.values())
                    self.logger.success(
                        f"Loaded RNAMotifScanX motifs into RSMViewer: "
                        f"{len(loaded_motifs)} motif types, {total_instances} instances"
                    )
                    self.logger.info(f"RMSX families available: {len(results)}")
                    return True

                self.logger.warning("RNAMotifScanX pipeline finished, but no motifs were loaded into RSMViewer.")
                return False
            except Exception as e:
                self.logger.error(f"RNAMotifScanX pipeline execution failed: {type(e).__name__}: {e}")
                return False
            finally:
                if query_override:
                    self.rmsx_query_file = previous_query_file

        if not self.rmsx_executable_path:
            self.logger.error("RNAMotifScanX executable is not configured")
            self.logger.info("Use: rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH]")
            return False

        if not self._is_executable_runnable(self.rmsx_executable_path):
            self.logger.error(f"RNAMotifScanX executable is not runnable: {self.rmsx_executable_path}")
            return False

        try:
            os.makedirs(self.rmsx_output_path, exist_ok=True)
        except Exception as e:
            self.logger.error(f"Could not create RNAMotifScanX output directory {self.rmsx_output_path}: {type(e).__name__}: {e}")
            return False

        work_dir = self.rmsx_working_dir or os.path.dirname(self.rmsx_executable_path)
        if not work_dir or not os.path.isdir(work_dir):
            self.logger.error(f"RNAMotifScanX working directory not found: {work_dir}")
            return False

        template_values = {
            'pdb_id': pdb_upper,
            'pdb_lower': pdb_upper.lower(),
            'output_dir': self.rmsx_output_path,
            'work_dir': work_dir,
        }

        args = []
        if self.rmsx_args_template:
            try:
                rendered = self.rmsx_args_template.format(**template_values)
            except KeyError as e:
                self.logger.error(f"Unknown placeholder in RNAMotifScanX args template: {e}")
                return False
            args.extend(shlex.split(rendered))

        if extra_args:
            args.extend(shlex.split(extra_args))

        command = [self.rmsx_executable_path] + args
        self.logger.info(f"Running RNAMotifScanX wrapper for {pdb_upper}...")
        self.logger.debug(f"RNAMotifScanX command: {' '.join(command)}")
        self.logger.debug(f"RNAMotifScanX cwd: {work_dir}")

        try:
            proc = subprocess.run(
                command,
                cwd=work_dir,
                capture_output=True,
                text=True,
                check=False,
            )
        except Exception as e:
            self.logger.error(f"Failed to execute RNAMotifScanX: {type(e).__name__}: {e}")
            return False

        if proc.returncode != 0:
            stderr_tail = (proc.stderr or '').strip().splitlines()[-10:]
            stdout_tail = (proc.stdout or '').strip().splitlines()[-10:]
            self.logger.error(f"RNAMotifScanX failed (exit code {proc.returncode})")
            if stderr_tail:
                self.logger.error("RNAMotifScanX stderr (last lines):")
                for line in stderr_tail:
                    self.logger.error(f"  {line}")
            elif stdout_tail:
                self.logger.error("RNAMotifScanX output (last lines):")
                for line in stdout_tail:
                    self.logger.error(f"  {line}")
            return False

        self.user_data_paths[7] = self.rmsx_output_path
        self._handle_source_by_id(7, self.rmsx_output_path)
        self.load_user_annotations_action('rnamotifscanx', pdb_upper, auto_pipeline=False)

        loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs()
        if loaded_motifs:
            total_instances = sum(len(info.get('motif_details', [])) for info in loaded_motifs.values())
            self.logger.success(f"Loaded RNAMotifScanX motifs into RSMViewer: {len(loaded_motifs)} motif types, {total_instances} instances")
            return True

        self.logger.warning("RNAMotifScanX finished, but no motifs were loaded into RSMViewer.")
        self.logger.info("Check that output files exist under the configured output directory and match RNAMotifScanX format.")
        return False

    # ── scan_prepared mode ────────────────────────────────────────────────
    def run_rmsx_scan_prepared(self, pdb_id: str, chains: str = '', compare: bool = False):
        """Run RNAMotifScanX against locally prepared inputs, off the GUI thread.

        Skips MC-Annotate/RNAVIEW, never reads or substitutes preannotated data,
        writes to a dedicated per-run output directory, reports progress by
        chain, supports cancellation, and loads only its own fresh output.
        """
        import datetime

        pdb_upper = str(pdb_id).strip().upper()
        if not pdb_upper:
            self.logger.error("PDB ID is required")
            self.logger.info("Usage: rmv_rmsx scan_prepared <PDB_ID> [CHAINS]")
            return False

        existing = getattr(self, '_rmsx_scan_thread', None)
        if existing is not None and existing.is_alive():
            self.logger.error(
                "A scan_prepared run is already in progress. Use 'rmv_rmsx scan_cancel' to stop it."
            )
            return False

        rmsx_cfg = dict(getattr(self, 'rmsx_pipeline_config', {}) or self._build_internal_rmsx_config())
        self.rmsx_pipeline_config = dict(rmsx_cfg)
        chain_list = [c for c in re.split(r'[\s,]+', str(chains or '').strip()) if c]

        base_out = os.path.abspath(os.path.expanduser(
            str(rmsx_cfg.get('output_dir', self.rmsx_output_path) or self.rmsx_output_path)
        ))
        stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        run_out = os.path.join(base_out, 'scan_prepared', f"{pdb_upper}_{stamp}")

        self._rmsx_scan_cancel = threading.Event()
        self.logger.info(
            f"scan_prepared: starting for {pdb_upper}"
            + (f" chains={chain_list}" if chain_list else " (all prepared chains)")
        )
        self.logger.info("Using locally prepared RMSX inputs; MC-Annotate and RNAVIEW are skipped.")

        worker = threading.Thread(
            target=self._scan_prepared_worker,
            args=(rmsx_cfg, pdb_upper, run_out, chain_list, bool(compare)),
            name=f"rmsx-scan-{pdb_upper}",
            daemon=True,
        )
        self._rmsx_scan_thread = worker
        worker.start()
        return True

    def cancel_rmsx_scan(self):
        """Signal an in-progress scan_prepared run to stop."""
        event = getattr(self, '_rmsx_scan_cancel', None)
        thread = getattr(self, '_rmsx_scan_thread', None)
        if event is None or thread is None or not thread.is_alive():
            self.logger.info("No scan_prepared run is currently active.")
            return False
        event.set()
        self.logger.warning("scan_prepared: cancellation requested; finishing current step...")
        return True

    def _scan_prepared_worker(self, rmsx_cfg, pdb_upper, run_out, chain_list, compare):
        """Background worker: runs scan subprocesses off the PyMOL GUI thread."""
        try:
            tools_dir = str(Path(__file__).parent / 'tools')
            if tools_dir not in sys.path:
                sys.path.insert(0, tools_dir)
            from rmsx_runner import run_scan_prepared  # type: ignore

            def progress(message):
                self.logger.info(f"[scan_prepared] {message}")

            report = run_scan_prepared(
                rmsx_cfg, pdb_upper, run_out,
                chains=chain_list or None,
                families=rmsx_cfg.get('motif_families'),
                progress_cb=progress,
                cancel_event=getattr(self, '_rmsx_scan_cancel', None),
            )
            self._rmsx_last_scan_report = report
            self._finish_scan_prepared(report, pdb_upper, compare)
        except Exception as exc:
            self.logger.error(f"scan_prepared failed: {type(exc).__name__}: {exc}")

    def _finish_scan_prepared(self, report, pdb_upper, compare):
        """Load newly generated results (main outcome handling)."""
        out_dir = report.get('output_dir', '')

        if report.get('cancelled'):
            self.logger.warning(
                f"scan_prepared CANCELLED for {pdb_upper}. Partial output was not loaded: {out_dir}"
            )
            return False

        if report.get('problems'):
            for problem in report['problems']:
                self.logger.warning(f"scan_prepared: {problem}")

        if report.get('failed_runs'):
            self.logger.error(
                f"scan_prepared: {len(report['failed_runs'])} run(s) FAILED for {pdb_upper}. "
                "Not loading results and NOT falling back to preannotated data."
            )
            for run in report['failed_runs'][:10]:
                self.logger.error(
                    f"  chain {run.get('chain')} / {run.get('family')}: "
                    f"exit {run.get('exit_code')} ({run.get('error')}); stderr: {run.get('stderr')}"
                )
            return False

        if not report.get('runs'):
            self.logger.error(
                f"scan_prepared produced no runs for {pdb_upper}. "
                "See problems above; preannotated data is NOT substituted."
            )
            return False

        # All attempted runs exited 0: load strictly from this run's output.
        self.rmsx_output_path = out_dir
        self.user_data_paths[7] = out_dir
        self._rmsx_skip_pipeline_load = True
        try:
            try:
                self._handle_source_by_id(7, out_dir)
            except Exception:
                pass
            self.load_user_annotations_action('rnamotifscanx', pdb_upper, auto_pipeline=False)
        finally:
            self._rmsx_skip_pipeline_load = False

        loaded = {}
        try:
            loaded = self.viz_manager.motif_loader.get_loaded_motifs() or {}
        except Exception:
            loaded = {}
        total_instances = sum(len(info.get('motif_details', [])) for info in loaded.values())

        self.logger.success("Loaded annotations from the newly generated RMSX output.")
        self.logger.info(
            f"scan_prepared summary for {pdb_upper}: families with hits={len(report['families'])}, "
            f"total hits={report['total_hits']}, motif types loaded={len(loaded)}, "
            f"instances={total_instances}"
        )
        self.logger.info(f"Results saved to: {out_dir}")
        if report['total_hits'] == 0:
            self.logger.info(
                "Scan completed successfully with zero hits (a valid result, not a failure)."
            )
        if compare:
            self._compare_scan_prepared_with_preannotated(report, pdb_upper)
        return True

    def _compare_scan_prepared_with_preannotated(self, report, pdb_upper):
        """Development check only: compare fresh scan output with preannotated data.

        This never alters, supplements, or filters the newly generated results;
        it just reports differences for investigation.
        """
        try:
            from .tools.rmsx_runner import copy_preannotated_results

            rmsx_cfg = getattr(self, 'rmsx_pipeline_config', {}) or {}
            import tempfile
            with tempfile.TemporaryDirectory(prefix='rmsx_preann_cmp_') as tmp_dir:
                source_dir = str(rmsx_cfg.get('pdb_prebuild_dir', '') or '')
                source_archive = str(rmsx_cfg.get('pdb_prebuild_archive', '') or '')
                copied = {'copied': 0}
                for source in [s for s in (source_dir, source_archive) if s]:
                    copied = copy_preannotated_results(str(source), pdb_upper, tmp_dir)
                    if copied.get('copied', 0):
                        break
                if not copied.get('copied', 0):
                    self.logger.info(
                        f"[compare] No preannotated dataset available for {pdb_upper}; skipping comparison."
                    )
                    return

                from .tools.rmsx_runner import _count_alignment_hits  # type: ignore
                self.logger.info(f"[compare] Fresh vs preannotated hit counts for {pdb_upper}:")
                families = set(report['families'].keys())
                for family in sorted(families):
                    folder = f"{family}_consensus"
                    pre_log = os.path.join(tmp_dir, folder, 'result_0_100_withbs.log')
                    pre_hits = _count_alignment_hits(pre_log) if os.path.isfile(pre_log) else 0
                    new_hits = report['families'].get(family, {}).get('hits', 0)
                    flag = 'match' if pre_hits == new_hits else 'DIFF -> investigate'
                    self.logger.info(f"  {family:<18} fresh={new_hits:<4} preannotated={pre_hits:<4} [{flag}]")
        except Exception as exc:
            self.logger.info(f"[compare] Comparison skipped: {type(exc).__name__}: {exc}")

    def load_structure_action(self, pdb_id_or_path, background_color=None,
                              database=None):
        """
        Load structure and automatically visualize all motifs.
        
        Args:
            pdb_id_or_path (str): PDB ID or file path
            background_color (str): Color for RNA backbone (default: 'gray80')
            database (str): Database to use ('atlas', 'rfam', or None for active)
        """
        try:
            self.logger.info(f"Loading structure: {pdb_id_or_path}")
            
            # Load and visualize with specified database
            motifs = self.viz_manager.load_and_visualize(
                pdb_id_or_path, 
                background_color,
                provider_id=database
            )
            
            if not motifs:
                self.logger.warning("No motifs found or error loading structure")
                return
            
            # Update UI state
            self.motif_visibility = {}
            for motif_type, info in motifs.items():
                self.motif_visibility[motif_type] = True
            
            self.logger.success(f"Loaded {len(motifs)} motif types")
            
        except Exception as e:
            self.logger.error(f"Failed to load structure: {e}")
    
    def fetch_motif_data_action(self, pdb_id, background_color=None):
        """
        Load motif data for a structure WITHOUT creating PyMOL objects (for rmv_load_motif).
        
        Handles multi-source loading (sources kept separate; no merging at rmv_db).
        Uses self.loaded_pdb as the structure name (set by rmv_fetch).
        
        Args:
            pdb_id (str): PDB ID already loaded in PyMOL
            background_color (str): Optional background color
        """
        try:
            # Set background color if specified
            if background_color:
                cmd.bg_color(background_color)
            
            # Structure name is the raw PDB name (set by rmv_fetch, no source suffix)
            # Motif objects get source suffixes, but the PDB structure is shared
            source_suffix = self._get_source_suffix()
            structure_name = self.loaded_pdb or pdb_id.lower()
            self.loaded_pdb_id = pdb_id.upper()
            pdb_id_upper = pdb_id.upper()
            
            # Check if we're in combine mode
            if self.current_source_mode == 'combine' and self.combined_source_ids:
                # Load and merge from multiple sources
                available_motifs = self._load_combined_motifs(
                    pdb_id_upper,
                    self.combined_source_ids
                )
                source_name = f"combined ({len(self.combined_source_ids)} sources)"
            else:
                # Load from single source
                from .database import get_source_selector
                source_selector = get_source_selector()
                
                if source_selector:
                    # Check if we're in user mode with specific tool selected
                    if self.current_source_mode == 'user' and self.current_user_tool:
                        # If custom data path is set, override the tool directory
                        _udp = self.user_data_paths.get(self.current_source_id)
                        if _udp and 'user' in source_selector.providers:
                            from pathlib import Path
                            user_prov_ref = source_selector.providers['user']
                            # Map GUI tool name to provider's internal tool name
                            tool_name_map = {
                                'rnamotifscanx': 'RNAMotifScanX',
                                'rmsx': 'RNAMotifScanX',
                                'rnamotifscan': 'RNAMotifScan',
                                'rms': 'RNAMotifScan',
                                'fr3d': 'fr3d',
                                'nobias': 'NoBIAS',
                            }
                            internal_tool = tool_name_map.get(self.current_user_tool.lower(), self.current_user_tool)
                            user_prov_ref.override_tool_dirs[internal_tool] = Path(_udp)
                            self.logger.debug(f"Override tool dir for {internal_tool}: {_udp}")
                        
                        # Apply filtering settings to the provider
                        user_prov = source_selector.providers.get('user')
                        if user_prov:
                            tool_lower = self.current_user_tool.lower() if self.current_user_tool else ''
                            if tool_lower in ['rms', 'rnamotifscan']:
                                user_prov.apply_rms_filtering = self.user_rms_filtering_enabled
                                user_prov.set_rms_custom_pvalues(self.user_rms_custom_pvalues)
                            elif tool_lower in ['rmsx', 'rnamotifscanx']:
                                user_prov.apply_rmsx_filtering = self.user_rmsx_filtering_enabled
                                user_prov.set_rmsx_custom_pvalues(self.user_rmsx_custom_pvalues)
                            elif tool_lower in ['nobias']:
                                user_prov.apply_nobias_filtering = self.user_nobias_filtering_enabled
                                user_prov.set_nobias_custom_pvalues(self.user_nobias_custom_pvalues)
                        
                        # Use tool-specific method to filter data
                        self.logger.debug(f"Using tool-specific loading: {self.current_user_tool}")
                        available_motifs = source_selector.get_motifs_for_pdb_and_tool(
                            pdb_id_upper, self.current_user_tool
                        )
                        source_name = self.current_user_tool.upper()
                    else:
                        # Get motif data from source selector (default)
                        self.logger.debug(f"Using default loading (mode={self.current_source_mode}, tool={self.current_user_tool})")
                        
                        # Pass specific source if one is selected (local or web)
                        source_override = None
                        if self.current_source_mode == 'local' and self.current_local_source:
                            source_override = self.current_local_source
                            self.logger.debug(f"Using specific local source: {source_override}")
                        elif self.current_source_mode == 'web' and self.current_web_source:
                            # Map web source to provider ID
                            web_source_map = {'bgsu': 'bgsu_api', 'rfam_api': 'rfam_api'}
                            source_override = web_source_map.get(self.current_web_source)
                            if source_override:
                                self.logger.debug(f"Using specific web source: {source_override}")
                        
                        # VERIFICATION: Check config.specific_source matches
                        from .database import get_config
                        config = get_config()
                        if config.specific_source and source_override != config.specific_source:
                            self.logger.debug(f"Adjusting source_override={source_override} to config.specific_source={config.specific_source}")
                            source_override = config.specific_source  # Use config value as authoritative
                        
                        self.logger.debug(f"Final source_override={source_override}, config.specific_source={config.specific_source}")
                        
                        available_motifs, source_used = source_selector.get_motifs_for_pdb(
                            pdb_id_upper,
                            source_override=source_override
                        )
                        source_name = source_used or "unknown"
                else:
                    # Fall back to active provider
                    registry = self.viz_manager.motif_loader._registry
                    provider = registry.get_active_provider()
                    if not provider:
                        self.logger.error("No database provider available")
                        return
                    
                    available_motifs = provider.get_motifs_for_pdb(pdb_id_upper)
                    source_name = provider.info.name if hasattr(provider, 'info') else 'unknown'
            
            if not available_motifs:
                self.logger.warning(f"No motifs found for {pdb_id}")
                return
            
            # --- Normalise generic keys (HL, IL, J3…) -> semantic types ----------
            # Atlas returns generic keys; this re-categorises them based on
            # per-instance annotation (same logic BGSU API uses).  Providers
            # that already return semantic keys (e.g. BGSU) pass through unchanged.
            available_motifs = _normalize_motif_groups(available_motifs)

            # Maintain a structure-local, source-aware table independently of
            # the legacy display dictionaries.  Combined-source consolidation
            # will move here once its provider labels are textual as well.
            if len(self.current_source_names) == 1:
                table = self.annotation_tables.setdefault(
                    pdb_id_upper, ConsolidatedAnnotationTable(self.jaccard_threshold, merge_enabled=False)
                )
                table.jaccard_threshold = self.jaccard_threshold
                table.merge_enabled = False
                table.remove_source(pdb_id_upper, self.current_source_names[0])
                table.add_annotations(
                    pdb_id_upper,
                    self.current_source_names[0],
                    available_motifs,
                    provenance={"provider": source_name},
                )
            
            # Use pre-built auth->label chain mapping (from CIF parsing in rmv_fetch)
            # When cif_use_auth=0, motif data has auth_asym_id chains but PyMOL has label_asym_id
            auth_to_label = self.auth_to_label_map if self.cif_use_auth == 0 else {}
            
            # Count total motifs
            total_count = sum(len(instances) for instances in available_motifs.values())
            # Report the public source name, never the internal adapter id.
            public_source = (
                self.current_source_names[0]
                if len(self.current_source_names) == 1
                else source_name
            )
            if len(self.current_source_names) == 1:
                self.logger.debug(
                    f"Found {total_count} motifs in {pdb_id} (source: {public_source})")
            
            # Process motifs for data access (WITHOUT creating PyMOL objects)
            motif_summary = {}
            from .utils.parser import SelectionParser
            
            for motif_type, instances in available_motifs.items():
                display_type = motif_type.split(':')[-1] if ':' in motif_type else motif_type
                display_type_upper = display_type.upper()
                
                # Convert to motif details format (same as _load_motif_type does)
                motif_details = []
                motif_list = []
                
                for instance in instances:
                    if hasattr(instance, 'residues') and instance.residues:
                        # Remap chain IDs if in label mode (cif_use_auth=0)
                        if auth_to_label:
                            for r in instance.residues:
                                if r.chain in auth_to_label:
                                    r.chain = auth_to_label[r.chain]
                        
                        # Include metadata from instance (contains chainbreak info)
                        from copy import deepcopy
                        metadata_to_store = deepcopy(instance.metadata) if hasattr(instance, 'metadata') and instance.metadata else {}
                        
                        motif_details.append({
                            'motif_id': instance.motif_id,
                            'instance_id': instance.instance_id,
                            'residues': [r.to_tuple() for r in instance.residues],
                            'annotation': instance.annotation,
                            'metadata': metadata_to_store,
                            '_source_suffix': source_suffix,
                            '_pdb_id': pdb_id_upper,
                            '_structure_name': structure_name,
                        })
                        
                        # Also build motif_list for selection string creation
                        legacy_entries = instance.to_legacy_format()
                        motif_list.extend(legacy_entries)
                
                # Build main_selection string (needed for show_motif_type to work)
                main_motif_sel = None
                if motif_list:
                    all_selections = []
                    for motif in motif_list:
                        chain = motif.get('chain')
                        residues = motif.get('residues')
                        sel = SelectionParser.create_selection_string(chain, residues)
                        if sel:
                            all_selections.append(f"({sel})")
                    
                    if all_selections:
                        combined_sel = " or ".join(all_selections)
                        main_motif_sel = f"(model {structure_name}) and ({combined_sel})"
                
                if motif_details:
                    if display_type_upper in motif_summary:
                        # Accumulate into existing entry (handles key casing variants)
                        existing = motif_summary[display_type_upper]
                        existing['motif_details'].extend(motif_details)
                        existing['motifs'].extend(motif_list)
                        existing['count'] += len(motif_details)
                        existing['display_names'][display_type] += len(motif_details)
                        # Rebuild combined selection
                        if main_motif_sel:
                            if existing['main_selection']:
                                existing['main_selection'] = f"{existing['main_selection']} or {main_motif_sel}"
                            else:
                                existing['main_selection'] = main_motif_sel
                        self.logger.debug(f"Loaded {len(motif_details)} more {display_type_upper} motifs (total: {existing['count']})")
                    else:
                        from collections import Counter as _Counter
                        motif_summary[display_type_upper] = {
                            'object_name': None,  # Will be created when rmv_show is called
                            'structure_name': structure_name,
                            'pdb_id': pdb_id_upper,
                            'count': len(motif_details),
                            'visible': False,
                            'motif_details': motif_details,
                            'motifs': motif_list,  # Needed to create PyMOL objects later
                            'main_selection': main_motif_sel,
                            'source_suffix': source_suffix,
                            'display_names': _Counter({display_type: len(motif_details)}),
                        }
                        self.logger.debug(f"Loaded {len(motif_details)} {display_type_upper} motifs")
            
            # Sort motif_details within each type by minimum residue number
            # (same as _load_motif_type does in loader.py)
            def _get_min_residue(detail):
                residues = detail.get('residues', [])
                if not residues:
                    return float('inf')
                min_resi = float('inf')
                for res in residues:
                    if isinstance(res, tuple) and len(res) >= 2:
                        resi = res[1]
                        if isinstance(resi, int):
                            min_resi = min(min_resi, resi)
                return min_resi if min_resi != float('inf') else float('inf')
            
            for mtype_key, mtype_info in motif_summary.items():
                mtype_info['motif_details'].sort(key=_get_min_residue)
            
            # --- Persist to the residue-based hierarchy cache (single-source loads) ---
            # Covers the plain rmv_db <source> path so rmv_summary/rmv_show <TYPE>
            # can also be driven by residue/chain matching.
            self._persist_single_source_hierarchy(pdb_id_upper, motif_summary)
            
            # Accumulate into existing loaded_motifs (supports cross-PDB and
            # multi-source workflows).  If the same motif type was already loaded
            # from the same PDB+source, those instances are replaced; otherwise
            # the new instances are appended.
            # NOTE: We intentionally do NOT delete data from other sources for
            # the same PDB.  All accumulated data is kept so that cross-source
            # superimposition (e.g. rmv_super K-TURN 1S72_S7, 1S72_S3) works.
            # The display commands (rmv_summary, rmv_show) filter to the
            # current source instead.
            existing_loaded = self.viz_manager.motif_loader.loaded_motifs

            # --- Clean sweep: remove ALL instances matching this PDB+suffix ---
            # This prevents stale categories from lingering when the new data
            # has fewer motif types than a previous load for the same combo.
            for key in list(existing_loaded.keys()):
                ex = existing_loaded[key]
                ex['motif_details'] = [
                    d for d in ex.get('motif_details', [])
                    if not (d.get('_pdb_id', ex.get('pdb_id', '')) == pdb_id_upper
                            and d.get('_source_suffix', ex.get('source_suffix', '')) == source_suffix)
                ]
                ex['count'] = len(ex['motif_details'])
                # Remove the key entirely if no instances remain
                if not ex['motif_details']:
                    del existing_loaded[key]

            for key, new_info in motif_summary.items():
                if key in existing_loaded:
                    ex = existing_loaded[key]
                    # Remove stale instances from the same PDB+source
                    new_suffix = new_info.get('source_suffix', '')
                    new_pdb = new_info.get('pdb_id', '')
                    ex['motif_details'] = [
                        d for d in ex['motif_details']
                        if not (d.get('_pdb_id', ex.get('pdb_id', '')) == new_pdb
                                and d.get('_source_suffix', ex.get('source_suffix', '')) == new_suffix)
                    ]
                    # Append new instances
                    ex['motif_details'].extend(new_info['motif_details'])
                    ex['motifs'].extend(new_info.get('motifs', []))
                    ex['count'] = len(ex['motif_details'])
                    # Rebuild combined selection
                    new_sel = new_info.get('main_selection', '')
                    if new_sel:
                        if ex.get('main_selection'):
                            ex['main_selection'] = f"{ex['main_selection']} or {new_sel}"
                        else:
                            ex['main_selection'] = new_sel
                    # Invalidate cached PyMOL object (will be recreated on rmv_show)
                    ex['object_name'] = None
                    ex['visible'] = False
                else:
                    existing_loaded[key] = new_info
            
            # Re-sort after accumulation
            for mtype_key, mtype_info in existing_loaded.items():
                mtype_info['motif_details'].sort(key=_get_min_residue)
            
            # CRITICAL: Also set structure_loader fields so rmv_save can find them
            self.viz_manager.structure_loader.current_structure = structure_name
            self.viz_manager.structure_loader.current_pdb_id = pdb_id_upper
            
            # Register this PDB+source combo for cross-PDB tracking
            if source_suffix:
                self.loaded_sources.add((pdb_id_upper, source_suffix))
            
            if motif_summary:
                self._report_loaded_families(
                    pdb_id, pdb_id_upper, source_suffix, "RNA3DMotifAtlas", motif_summary
                )
            else:
                self.logger.warning(f"No valid motifs found for {pdb_id}")
                
        except Exception as e:
            self.logger.error(f"Failed to load motif data: {str(e)}")

    def _print_family_table(self, family_rows: List[Tuple[str, int]]) -> None:
        """Print the three-column rmv_db family table.

        Columns: [converted name, no spaces] | [annotation name] | [# loaded].
        The converted name is the copy-paste-friendly form users pass to
        rmv_select; the annotation name is the source's own spelling.
        """
        from .database.motif_aliases import converted_family_name

        rows = [
            (converted_family_name(name), name, count) for name, count in family_rows
        ]
        conv_header, name_header, count_header = (
            "SELECTABLE FAMILY NAME", "SOURCE ANNOTATION NAME", "# LOADED MOTIFS",
        )
        conv_width = max([len(conv) for conv, _, _ in rows] + [len(conv_header), len("Total")])
        name_width = max([len(name) for _, name, _ in rows] + [len(name_header)])
        print("")
        print(f"  {conv_header.ljust(conv_width)}   {name_header.ljust(name_width)}   {count_header}")
        print(f"  {'-' * conv_width}   {'-' * name_width}   {'-' * len(count_header)}")
        for conv, name, count in rows:
            print(f"  {conv.ljust(conv_width)}   {name.ljust(name_width)}   {count}")
        total = sum(count for _, _, count in rows)
        print(f"  {'Total'.ljust(conv_width)}   {' ' * name_width}   {total}")
        print("")
        print("  (copy this exact name into rmv_select)")
        print("  Input is normalized, so minor differences in spaces, hyphens,")
        print("  capitalization, or parentheses are accepted. For the clearest")
        print("  results, copy the first-column name.")

    def _family_counts_from_table(self, structure_id: str, source_filter: Optional[List[str]] = None):
        """Return [(display_name, row_count), ...] for one structure's table.

        Each consolidated row is counted under every canonical family that any
        of its source labels denotes, so the totals match what rmv_select would
        select. The representative display name is the most common original
        spelling collected for that family.

        When ``source_filter`` is given, only rows touched by one of those
        sources are considered, and only the labels contributed by those
        sources are used to determine family membership. This is what makes a
        single-source load (e.g. just Rfam) report its own family breakdown
        instead of the whole cross-source table accumulated so far.
        """
        from collections import Counter
        from .database.motif_aliases import canonical_motif

        table = self.annotation_tables.get(structure_id)
        if not table:
            return []
        filter_set = set(source_filter) if source_filter else None
        counts: Counter = Counter()
        label_pool: Dict[str, Counter] = {}
        for row in table.rows:
            if getattr(row, 'structure_id', structure_id) != structure_id:
                continue
            if filter_set is not None and not (set(row.source_annotations) & filter_set):
                continue
            families_in_row = set()
            for source_map in (row.source_hierarchy, row.source_annotations):
                for source_name, labels in source_map.items():
                    if filter_set is not None and source_name not in filter_set:
                        continue
                    for label in labels:
                        if not label:
                            continue
                        canonical = canonical_motif(label)
                        label_pool.setdefault(canonical, Counter())[label] += 1
                        families_in_row.add(canonical)
            for canonical in families_in_row:
                counts[canonical] += 1
        result = []
        # Canonical families whose raw source spelling is inconsistent
        # (a source typo or an index suffix that rmv_select already unifies).
        canonical_display = {
            "RIBOSOMAL LSU H95": "Ribosomal LSU H95",
            "RIGHT-ANGLE": "Right-angle",
            "TWIST-UP": "Twist-up",
        }
        for canonical, count in counts.items():
            representative = canonical_display.get(
                canonical, label_pool[canonical].most_common(1)[0][0])
            result.append((representative, count))
        result.sort(key=lambda item: (-item[1], item[0].upper()))
        return result

    def _family_source_breakdown(self, structure_id: str, source_names: Optional[List[str]] = None):
        """Return source-specific family rows without changing annotation data."""
        from collections import Counter, defaultdict
        from .database.motif_aliases import canonical_motif, converted_family_name

        table = self.annotation_tables.get(structure_id)
        if not table:
            return [], [], {}

        present: List[str] = []
        for row in table.rows:
            if getattr(row, 'structure_id', structure_id) != structure_id:
                continue
            for source in row.source_annotations:
                if source not in present:
                    present.append(source)
        source_order = list(source_names or present)
        for source in present:
            if source not in source_order:
                source_order.append(source)
        filter_set = set(source_order)

        fam_source_count: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        fam_source_labels: Dict[str, Dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
        totals: Dict[str, int] = defaultdict(int)
        canonical_display = {
            "RIBOSOMAL LSU H95": "Ribosomal LSU H95",
            "RIGHT-ANGLE": "Right-angle",
            "TWIST-UP": "Twist-up",
        }

        for row in table.rows:
            if getattr(row, 'structure_id', structure_id) != structure_id:
                continue
            for source in set(row.source_annotations) & filter_set:
                totals[source] += 1
                labels = list(row.source_annotations.get(source, ())) + list(
                    row.source_hierarchy.get(source, ())
                )
                families_here = set()
                for label in labels:
                    if not label:
                        continue
                    canonical = canonical_motif(label)
                    families_here.add(canonical)
                    fam_source_labels[canonical][source][label] += 1
                for canonical in families_here:
                    fam_source_count[canonical][source] += 1

        source_families = {}
        for source in source_order:
            rows = []
            for canonical, per_source_counts in fam_source_count.items():
                count = per_source_counts.get(source, 0)
                if not count:
                    continue
                labels = sorted(fam_source_labels[canonical][source])
                rows.append({
                    "selectable": converted_family_name(
                        canonical_display.get(canonical, canonical)
                    ),
                    "annotation": " / ".join(labels),
                    "count": count,
                    "canonical": canonical,
                })
            rows.sort(key=lambda row: (-row["count"], row["selectable"].upper()))
            source_families[source] = rows
        return source_order, source_families, dict(totals)

    @staticmethod
    def _print_wrapped_family_table(source, rows, total) -> None:
        """Print one compact source table with aligned wrapped text columns."""
        # Selectable names must stay on one line (they are copied verbatim), so
        # this column is never capped/wrapped; only the annotation column wraps.
        selectable_width = max([len("SELECTABLE NAME")] + [len(r["selectable"]) for r in rows])
        annotation_width = min(42, max([len("ANNOTATION NAME")] + [len(r["annotation"]) for r in rows]))
        count_width = max(len("COUNT"), len(str(total)))
        line = f"{'SELECTABLE NAME':<{selectable_width}}  {'ANNOTATION NAME':<{annotation_width}}  {'COUNT':>{count_width}}"
        rule = f"{'-' * selectable_width}  {'-' * annotation_width}  {'-' * count_width}"
        print(f"\nSource: {source}")
        if not rows:
            print("No annotations reported for this source.")
            return
        print(f"Raw annotations: {total} | Families: {len(rows)}\n")
        print(line)
        print(rule)
        for row in rows:
            left = [row["selectable"]]
            right = [row["annotation"][i:i + annotation_width] for i in range(0, len(row["annotation"]), annotation_width)]
            height = max(len(left), len(right))
            for index in range(height):
                selectable = left[index] if index < len(left) else ""
                annotation = right[index] if index < len(right) else ""
                count = str(row["count"]) if index == 0 else ""
                print(f"{selectable:<{selectable_width}}  {annotation:<{annotation_width}}  {count:>{count_width}}")
        print(f"{'Total':<{selectable_width}}  {'':<{annotation_width}}  {total:>{count_width}}")

    def _report_loaded_families(self, pdb_id, pdb_id_upper, source_suffix, source_fallback, motif_summary) -> None:
        """Log the load summary and print the per-source family table."""
        tag = f"{pdb_id_upper}{source_suffix}" if source_suffix else pdb_id_upper
        source_list = ", ".join(self.current_source_names) if self.current_source_names else source_fallback
        source_order, source_families, totals = self._family_source_breakdown(
            pdb_id_upper, self.current_source_names
        )

        row_count = sum(totals.values())
        displayed_rows = [
            row for source in source_order for row in source_families.get(source, [])
        ]

        if displayed_rows:
            self.logger.success(
                f"Loaded {row_count} raw annotation(s) from {len(source_order)} "
                f"source(s) for {pdb_id} via {source_list} (tag: {tag})"
            )
            for source in source_order:
                self._print_wrapped_family_table(
                    source,
                    source_families.get(source, []),
                    totals.get(source, 0),
                )
            print("\nCopy a SELECTABLE NAME into rmv_select.")
            print("Input normalization accepts minor differences in capitalization, spaces,")
            print("hyphens, and parentheses.")
            first_source = next(
                (source for source in source_order if source_families.get(source)),
                None,
            )
        else:
            # Fallback for providers that return raw motifs before table rows exist.
            family_rows = []
            for key, info in motif_summary.items():
                names = info.get('display_names')
                display_name = names.most_common(1)[0][0] if names else key
                family_rows.append((display_name, info.get('count', 0)))
            family_rows.sort(key=lambda item: (-item[1], item[0].upper()))
            if not row_count:
                row_count = sum(count for _, count in family_rows)
            self.logger.success(
                f"Loaded {row_count} raw annotation(s) from {len(source_order)} "
                f"source(s) for {pdb_id} via {source_list} (tag: {tag})"
            )
            fallback_source = source_order[0] if source_order else source_fallback
            fallback_rows = [
                {
                    "selectable": converted_family_name(name),
                    "annotation": name,
                    "count": count,
                }
                for name, count in family_rows
            ]
            self._print_wrapped_family_table(fallback_source, fallback_rows, row_count)
            print("\nCopy a SELECTABLE NAME into rmv_select.")
            print("Input normalization accepts minor differences in capitalization, spaces,")
            print("hyphens, and parentheses.")

    def _auto_color_motifs_on_structure(self, structure_name: str):
        """Color all loaded motif residues on the base PDB structure.

        Each motif type gets its unique color.  No PyMOL objects are created;
        the coloring is applied directly on the structure.
        """
        from .utils.parser import SelectionParser

        loaded = self.viz_manager.motif_loader.get_loaded_motifs()
        if not loaded:
            return

        # Ensure structure is visible with gray base
        self.viz_manager.cmd.enable(structure_name)
        self.viz_manager.cmd.show('cartoon', f"model {structure_name} and polymer.nucleic")
        self.viz_manager.cmd.set('cartoon_nucleic_acid_mode', 4, f"model {structure_name}")

        self.viz_manager.cmd.set('cartoon_tube_radius', 0.37, f"model {structure_name}")
        self.viz_manager.cmd.color('gray80', f"model {structure_name}")

        # Filter to current PDB + source
        fpdb = (self.loaded_pdb_id or '').upper()
        fsuf = self._get_source_suffix()

        for mtype, info in loaded.items():
            details = info.get('motif_details', [])
            if fpdb:
                details = [
                    d for d in details
                    if d.get('_pdb_id', info.get('pdb_id', '')) == fpdb
                    and d.get('_source_suffix', info.get('source_suffix', '')) == fsuf
                ]
            for detail in details:
                residues = detail.get('residues', [])
                if not residues:
                    continue
                chain_residues = {}
                for res in residues:
                    if isinstance(res, tuple) and len(res) >= 3:
                        chain = res[2]
                        chain_residues.setdefault(chain, []).append(res[1])
                for chain, resi_list in chain_residues.items():
                    sel = SelectionParser.create_selection_string(chain, sorted(resi_list))
                    if sel:
                        instance_sel = f"(model {structure_name}) and ({sel})"
                        colors.set_motif_color_in_pymol(
                            self.viz_manager.cmd, instance_sel, mtype)

        self.logger.info("Motif regions highlighted on structure")

    def _persist_single_source_hierarchy(self, pdb_id_upper: str, motif_summary: dict) -> None:
        """Record every instance from a single-source load into the hierarchy cache.

        Stores the final (possibly semantic) label at hierarchy level 1, and
        - when the instance also carries its original generic parent type
        (``_generic_type`` from ATLAS/BGSU-fallback categorisation, or
        ``loop_type`` set directly by the BGSU provider) - the generic label
        at level 1 and the specific label at level 2, so both survive in the
        table (e.g. BGSU hierarchy 1 = IL, BGSU hierarchy 2 = K-TURN).
        """
        try:
            from .database.motif_hierarchy_cache import get_hierarchy_cache, residue_key
            if len(self.current_source_names) != 1:
                return
            source_key = self.current_source_names[0]
            src_label = source_key
            cache = get_hierarchy_cache()
            for mtype_key, mtype_info in motif_summary.items():
                for detail in mtype_info['motif_details']:
                    pairs = [
                        (r[0], r[1]) for r in detail.get('residues', [])
                        if isinstance(r, tuple) and len(r) >= 2
                    ]
                    if not pairs:
                        continue
                    meta = detail.get('metadata') or {}
                    generic = meta.get('_generic_type') or meta.get('loop_type')
                    r_key = residue_key(pairs)
                    if generic and not _is_generic_equivalent(generic, mtype_key):
                        cache.record(pdb_id_upper, r_key, source_key, src_label,
                                     generic.upper(), rank=0, hierarchy_level=1,
                                     jaccard_threshold=self.jaccard_threshold)
                        cache.record(pdb_id_upper, r_key, source_key, src_label,
                                     mtype_key, rank=0, hierarchy_level=2,
                                     jaccard_threshold=self.jaccard_threshold)
                    else:
                        cache.record(pdb_id_upper, r_key, source_key, src_label,
                                     mtype_key, rank=0, hierarchy_level=1,
                                     jaccard_threshold=self.jaccard_threshold)
        except Exception as exc:
            self.logger.debug(f"Hierarchy cache update skipped: {exc}")
    
    def _copy_preannotated_rmsx(self, rmsx_cfg: Dict, pdb_id: str) -> Dict:
        """Copy preannotated RMSX results, preferring the extracted folder.

        Tries the extracted ``rmsx_work_default`` directory first (fast) and
        falls back to the compressed archive when the folder is missing or has
        no data for this PDB.
        """
        from .tools.rmsx_runner import copy_preannotated_results
        configured_dir = str(rmsx_cfg.get('pdb_prebuild_dir', '') or '')
        configured_archive = str(rmsx_cfg.get('pdb_prebuild_archive', '') or '')
        pdb_upper = pdb_id.strip().upper()
        ordered_sources = []
        if configured_dir and os.path.isdir(configured_dir):
            ordered_sources.append(configured_dir)
        if configured_archive:
            ordered_sources.append(configured_archive)
        copied = {'copied': 0, 'output_dir': self.rmsx_output_path}
        for source in ordered_sources:
            copied = copy_preannotated_results(str(source), pdb_upper, self.rmsx_output_path)
            if copied.get('copied', 0):
                break
        return copied

    def _ensure_rmsx_preannotated(self, pdb_id: str) -> bool:
        """Copy preannotated RNAMotifScanX results for ``pdb_id`` into the
        working output directory and point ``user_data_paths[7]`` at it.

        Combined-source loading (``rmv_db RNA3DMotifAtlas,RNAMotifScanX``)
        fetches each source directly via ``_fetch_from_single_source`` and
        therefore bypasses the ingestion that ``load_user_annotations_action``
        performs for single-source mode. This helper reproduces just the
        preannotated-copy step so RMSX rows are available in combine mode too.
        Returns True when data was made available.
        """
        try:
            rmsx_cfg = getattr(self, 'rmsx_pipeline_config', {}) or self._build_internal_rmsx_config()
            self.rmsx_pipeline_config = dict(rmsx_cfg)
            data_mode = str(rmsx_cfg.get('data_mode', 'preannotated')).strip().lower()
            if data_mode not in ('preannotated', 'cache', 'cached'):
                # Non-preannotated modes are handled by the single-source path.
                return bool(self.user_data_paths.get(7))
            copied = self._copy_preannotated_rmsx(rmsx_cfg, pdb_id)
            if copied.get('copied', 0):
                self.user_data_paths[7] = self.rmsx_output_path
                self.logger.info(
                    f"Using {copied['copied']} preannotated RNAMotifScanX result(s) for {pdb_id}."
                )
                return True
            self.logger.warning(
                f"No preannotated RNAMotifScanX results found for {pdb_id}."
            )
        except Exception as error:
            self.logger.warning(f"Could not prepare preannotated RMSX data: {error}")
        return False

    def _load_combined_motifs(self, pdb_id: str, source_ids: List[int]):
        """
        Load motifs from multiple sources WITHOUT merging.

        Pipeline:
        1. Fetch raw motifs from each source (labels used exactly as reported)
        2. annotation merging (right-to-left, priority = source order)
        
        Args:
            pdb_id: PDB ID to fetch
            source_ids: List of source IDs in priority order (first = highest)
        
        Returns:
            Merged motif dictionary: {motif_type: [MotifInstance, ...]}
        """
        try:
            from .database.config import SOURCE_ID_MAP

            # RNAMotifScanX (source 7) needs its preannotated results ingested
            # before _fetch_from_single_source can read them in combine mode.
            if 7 in source_ids:
                self._ensure_rmsx_preannotated(pdb_id)
            
            pdb_id = pdb_id.upper()
            
            # --- Step 1: Fetch raw motifs from each source ---
            self.logger.debug(f"Fetching motifs from {len(source_ids)} sources...")
            raw_sources = {}
            source_labels = []
            
            for sid in source_ids:
                info = SOURCE_ID_MAP.get(sid, {})
                label = info.get('name', f'Source {sid}')
                source_labels.append(label)
                
                motifs = self._fetch_from_single_source(pdb_id, sid)
                if motifs:
                    raw_sources[sid] = motifs
                    total = sum(len(v) for v in motifs.values())
                    self.logger.info(
                        f"  {label}: {total} raw annotations in {len(motifs)} categories")
                else:
                    raw_sources[sid] = {}
                    self.logger.warning(f"  {label}: no motifs found")
            
            if not any(raw_sources.values()):
                self.logger.error("No motifs found from any source")
                return {}
            
            # No homolog/representative-based renaming: labels are shown exactly as fetched.
            
            # --- Step 1.5: Stamp source origin on each instance ---
            for sid in source_ids:
                info = SOURCE_ID_MAP.get(sid, {})
                label = info.get('name', f'Source {sid}')
                for mtype, instances in raw_sources.get(sid, {}).items():
                    for inst in instances:
                        if inst.metadata is None:
                            inst.metadata = {}
                        inst.metadata['_source_id'] = sid
                        inst.metadata['_source_label'] = label

            # No within-source filtering at rmv_db: every source's annotations
            # are kept exactly as fetched (professor's "keep all records
            # separately"). The raw consolidated table indexes by residue-set
            # identity, so byte-for-byte identical annotations still map to one
            # row without dropping any distinct/contained annotation.
            self.dedup_stats = {}

            # Populate the new canonical table from each raw source before
            # the legacy annotation mergingr collapses source-specific records.
            canonical_by_legacy_id = {
                3: "RNA3DMotifAtlas",
                4: "Rfam",
                5: "FR3D",
                7: "RNAMotifScanX",
            }
            table = self.annotation_tables.setdefault(
                pdb_id, ConsolidatedAnnotationTable(self.jaccard_threshold, merge_enabled=False)
            )
            table.merge_enabled = False
            selected_names = self.current_source_names
            for index, sid in enumerate(source_ids):
                source_name = (
                    selected_names[index]
                    if index < len(selected_names)
                    else canonical_by_legacy_id.get(sid)
                )
                if source_name and raw_sources.get(sid):
                    table.jaccard_threshold = self.jaccard_threshold
                    table.remove_source(pdb_id, source_name)
                    table.add_annotations(
                        pdb_id,
                        source_name,
                        raw_sources[sid],
                        provenance={"legacy_provider": str(sid)},
                    )
            
            # --- Step 2: NO merging at rmv_db ---
            # Professor's design: rmv_db keeps every source's annotations
            # separate. Containment/Jaccard merging is deferred to rmv_select
            # (same family) and rmv_combine_groups (across families). Here we simply
            # concatenate the raw per-type instances so nothing is dropped and
            # a contained annotation (e.g. an isolated pair inside a sarcin-
            # ricin loop) survives until the user selects/combines.
            combined: Dict[str, List] = {}
            for sid in source_ids:
                for mtype, instances in raw_sources.get(sid, {}).items():
                    combined.setdefault(mtype.upper(), []).extend(instances)

            total = sum(len(v) for v in combined.values())
            self.logger.debug(
                f"Loaded {total} raw annotation(s) in {len(combined)} categories "
                f"from {len(source_ids)} sources (no merge at rmv_db)")
            return combined
            
        except Exception as e:
            self.logger.error(f"Failed to combine motifs: {e}")
            import traceback
            traceback.print_exc()
            return {}
    
    def _fetch_from_single_source(self, pdb_id: str, source_id: int):
        """
        Fetch motifs from a single source by its ID.
        
        Maps source IDs to the correct provider and fetches motifs.
        
        Args:
            pdb_id: PDB ID to fetch
            source_id: Source ID (1-8)
        
        Returns:
            Dict mapping motif_type -> [MotifInstance, ...], or empty dict
        """
        from .database.config import SOURCE_ID_MAP
        
        info = SOURCE_ID_MAP.get(source_id)
        if not info:
            return {}
        
        source_type = info['type']
        
        try:
            if source_type in ('local', 'web'):
                # Use source selector with explicit override
                from .database import get_source_selector
                source_selector = get_source_selector()
                if not source_selector:
                    return {}
                
                # Map source ID to provider ID
                if source_type == 'local':
                    provider_id = info.get('subtype')  # 'atlas' or 'rfam'
                else:
                    # Web sources
                    web_map = {'bgsu': 'bgsu_api', 'bgsu_api': 'bgsu_api', 'rfam_api': 'rfam_api'}
                    provider_id = web_map.get(info.get('subtype'))
                
                if provider_id and provider_id in source_selector.providers:
                    return source_selector.providers[provider_id].get_motifs_for_pdb(pdb_id)
                else:
                    # Try via source_override
                    motifs, _ = source_selector.get_motifs_for_pdb(
                        pdb_id, source_override=provider_id
                    )
                    return motifs
                    
            elif source_type == 'user':
                # User annotations (FR3D, RMS, RMSX)
                from .database.user_annotations import UserAnnotationProvider
                plugin_dir = Path(__file__).parent
                user_dir = plugin_dir / 'database' / 'user_annotations'
                provider = UserAnnotationProvider(str(user_dir))
                tool = info.get('tool')
                if tool:
                    provider.set_active_tool(tool)
                    # If custom data path is set for this source, override the tool directory
                    _udp = self.user_data_paths.get(source_id)
                    if _udp:
                        tool_name_map = {
                            'rnamotifscanx': 'RNAMotifScanX',
                            'rmsx': 'RNAMotifScanX',
                            'rnamotifscan': 'RNAMotifScan',
                            'rms': 'RNAMotifScan',
                            'fr3d': 'fr3d',
                            'nobias': 'NoBIAS',
                        }
                        internal_tool = tool_name_map.get(tool.lower(), tool)
                        provider.override_tool_dirs[internal_tool] = Path(_udp)
                    # Apply p-value filtering settings for RMS/RMSX/NoBIAS
                    tool_lower = tool.lower()
                    if tool_lower in ['rms', 'rnamotifscan']:
                        provider.apply_rms_filtering = self.user_rms_filtering_enabled
                        provider.set_rms_custom_pvalues(self.user_rms_custom_pvalues)
                    elif tool_lower in ['rmsx', 'rnamotifscanx']:
                        provider.apply_rmsx_filtering = self.user_rmsx_filtering_enabled
                        provider.set_rmsx_custom_pvalues(self.user_rmsx_custom_pvalues)
                    elif tool_lower in ['nobias']:
                        provider.apply_nobias_filtering = self.user_nobias_filtering_enabled
                        provider.set_nobias_custom_pvalues(self.user_nobias_custom_pvalues)
                        provider.set_rmsx_custom_pvalues(self.user_rmsx_custom_pvalues)
                return provider.get_motifs_for_pdb(pdb_id)
            
        except Exception as e:
            self.logger.warning(f"Error fetching from source {source_id}: {e}")
        
        return {}
    
    def load_user_annotations_action(
        self,
        tool,
        pdb_id,
        auto_pipeline: bool = True,
        force_pipeline_refresh: bool = False,
        rmsx_query_models_dir: str = '',
    ):
        """
        Load motifs from user-uploaded annotation files.
        
        Args:
            tool (str): Tool name ('fr3d', 'rnamotifscan')
            pdb_id (str): PDB ID to load annotations for
        """
        try:
            from .database.user_annotations import UserAnnotationProvider
            
            # Initialize user annotation provider (always use default root)
            plugin_dir = Path(__file__).parent
            user_annotations_dir = plugin_dir / 'database' / 'user_annotations'
            provider = UserAnnotationProvider(str(user_annotations_dir))
            self.logger.debug(f"User annotation loader initialized at {user_annotations_dir}")
            
            # SET ACTIVE TOOL FILTER BEFORE LOADING!
            provider.set_active_tool(tool)
            self.logger.debug(f"Active user annotation tool set to: {tool}")

            tool_lower = tool.lower() if tool else ''

            # -- FR3D (Source 5) ---------------------------------------------
            # FR3D searches are driven by run_fr3d_search(), which always calls
            # this method with auto_pipeline=False and points source 5 at the
            # freshly ingested CSV directory. There is no in-line auto-run here.

            # -- RMSX pipeline run (prebuilt/cache-aware by default) -----------
            # Source 7 reuses available/prebuilt results unless explicitly
            # forced to run fresh. scan_prepared sets _rmsx_skip_pipeline_load so
            # this block never copies preannotated data over its fresh output.
            if tool_lower in ['rmsx', 'rnamotifscanx'] and not getattr(self, '_rmsx_skip_pipeline_load', False):
                rmsx_cfg = getattr(self, 'rmsx_pipeline_config', {}) or self._build_internal_rmsx_config()
                self.rmsx_pipeline_config = dict(rmsx_cfg)
                data_mode = str(rmsx_cfg.get('data_mode', 'preannotated')).strip().lower()
                if data_mode in ('preannotated', 'cache', 'cached') and not force_pipeline_refresh:
                    try:
                        copied = self._copy_preannotated_rmsx(rmsx_cfg, pdb_id)
                        if copied.get('copied', 0):
                            self.user_data_paths[7] = self.rmsx_output_path
                            self.logger.debug(
                                f"Using {copied['copied']} preannotated RNAMotifScanX result(s); no executable run."
                            )
                        else:
                            self.logger.warning(
                                "No preannotated RMSX results found. Set data_mode to run_from_scratch "
                                "in config/rmsx_config.json to execute RNAMotifScanX."
                            )
                            return
                    except Exception as error:
                        self.logger.error(f"Could not load preannotated RMSX data: {error}")
                        return
                else:
                    if not self.ensure_rmsx_runtime_ready(auto_setup=True):
                        return
                if rmsx_cfg and not (data_mode in ('preannotated', 'cache', 'cached') and not force_pipeline_refresh):
                    try:
                        tools_dir = str(Path(__file__).parent / 'tools')
                        if tools_dir not in sys.path:
                            sys.path.insert(0, tools_dir)
                        from rmsx_runner import run_pipeline as rmsx_run  # type: ignore
                        pdb_upper = pdb_id.strip().upper()
                        force_fresh = bool(force_pipeline_refresh)

                        run_cfg = dict(rmsx_cfg)
                        if rmsx_query_models_dir:
                            run_cfg['query_motifs_dir'] = os.path.abspath(os.path.expanduser(rmsx_query_models_dir))
                        run_cfg['auto_download_cif'] = False
                        out_dir = os.path.abspath(os.path.expanduser(str(run_cfg.get('output_dir', self.rmsx_output_path))))
                        run_cfg['cif_input_dir'] = out_dir

                        local_pdb = self._prepare_local_pdb_for_rmsx(
                            pdb_upper, out_dir, force_refresh=force_fresh
                        )
                        if local_pdb:
                            run_cfg['auto_download_pdb'] = False

                        families = rmsx_cfg.get('motif_families', [])
                        mode_text = 'fresh' if force_fresh else 'incremental/prebuilt-aware'
                        if rmsx_query_models_dir:
                            self.logger.info(f"Using RMSX query model directory override: {run_cfg.get('query_motifs_dir')}")
                        self.logger.info(
                            f"Running RMSX pipeline for {pdb_upper} ({mode_text}, no external download) in {out_dir}..."
                        )
                        existing = (
                            {}
                            if data_mode in ('preannotated', 'cache', 'cached') and not force_fresh
                            else rmsx_run(run_cfg, pdb_upper, force_fresh=force_fresh)
                        )

                        if existing:
                            self.logger.info(
                                f"RMSX run complete: {len(existing)}/{len(families) or len(existing)} families"
                            )
                        else:
                            self.logger.warning(
                                "RMSX fresh run produced no result files. "
                                "Ensure rmsx_executable, mc_annotate_executable, query_motifs_dir are configured "
                                "and CIF is available locally (auto_download_cif is forced off for source 7)."
                            )
                            return
                    except Exception as _rmsx_e:
                        self.logger.debug(f"RMSX pipeline check error: {_rmsx_e}")
            # ----------------------------------------------------------------

            # If custom data path is set for this source, override the tool directory
            _udp = self.user_data_paths.get(self.current_source_id)
            if _udp:
                tool_name_map = {
                    'rnamotifscanx': 'RNAMotifScanX',
                    'rmsx': 'RNAMotifScanX',
                    'rnamotifscan': 'RNAMotifScan',
                    'rms': 'RNAMotifScan',
                    'fr3d': 'fr3d',
                    'nobias': 'NoBIAS',
                }
                internal_tool = tool_name_map.get(tool_lower, tool)
                provider.override_tool_dirs[internal_tool] = Path(_udp)
                if tool_lower == 'fr3d':
                    self.logger.info(f"Using custom FR3D data path: {_udp}")
                else:
                    self.logger.debug(f"Using custom data path for {tool.upper()}: {_udp}")
            
            # Set filtering state based on current settings (for RMS, RMSX, and NoBIAS)
            if tool_lower in ['rms', 'rnamotifscan']:
                provider.apply_rms_filtering = self.user_rms_filtering_enabled
                provider.set_rms_custom_pvalues(self.user_rms_custom_pvalues)
            elif tool_lower in ['rmsx', 'rnamotifscanx']:
                provider.apply_rmsx_filtering = self.user_rmsx_filtering_enabled
                provider.set_rmsx_custom_pvalues(self.user_rmsx_custom_pvalues)
            elif tool_lower in ['nobias']:
                provider.apply_nobias_filtering = self.user_nobias_filtering_enabled
                provider.set_nobias_custom_pvalues(self.user_nobias_custom_pvalues)
            
            # Get motifs
            pdb_id_upper = pdb_id.upper()
            available_motifs = provider.get_motifs_for_pdb(pdb_id_upper)
            
            if not available_motifs:
                self.logger.warning(f"No {tool.upper()} annotation files found for {pdb_id}")
                if tool_lower == 'fr3d':
                    fr3d_path = Path(_udp) if _udp else (user_annotations_dir / 'fr3d')
                    has_matching_csv = False
                    try:
                        pdb_prefix = str(pdb_id or '').strip().upper()
                        has_matching_csv = any(
                            p.is_file() and p.suffix.lower() == '.csv' and p.name.upper().startswith(pdb_prefix)
                            for p in fr3d_path.iterdir()
                        ) if fr3d_path.exists() else False
                    except Exception:
                        has_matching_csv = False

                    self.logger.info(f"Checked FR3D annotation path: {fr3d_path}")
                    if has_matching_csv:
                        self.logger.info(
                            "Matching FR3D CSV file was found, but it contained 0 loadable motif candidates for this query."
                        )
                    self.logger.info("Expected either FR3D CSV motif files or FR3D pairwise TXT output.")
                    self.logger.info("If you just ran rmv_fr3d, confirm the output directory contains a matching file.")
                else:
                    rmsx_cfg = getattr(self, 'rmsx_pipeline_config', {}) if tool_lower in ['rmsx', 'rnamotifscanx'] else {}
                    configured_out = str(rmsx_cfg.get('output_dir', '') or '').strip()
                    expected_root = _udp or configured_out or str(user_annotations_dir / 'RNAMotifScanX')
                    if tool_lower in ['rmsx', 'rnamotifscanx']:
                        self.logger.info(f"Checked RMSX path: {expected_root}")
                        self.logger.info(
                            "Expected files per family, e.g. "
                            "<output_dir>/k-turn_consensus/result_0_100_withbs.log"
                        )
                    else:
                        self.logger.info(f"Please place files in: {expected_root}")
                return
            
            # Structure name is the raw PDB name (set by rmv_fetch, no source suffix)
            source_suffix = self._get_source_suffix()
            structure_name = self.loaded_pdb or pdb_id_upper.lower()
            self.loaded_pdb_id = pdb_id_upper
            
            # Map numeric chain IDs to actual PyMOL chain IDs
            # FR3D uses numeric chains like "1", but PyMOL uses letters like "A"
            # RMSX/RMS use "0" to represent the chain in annotation data
            # When cif_use_auth=0 (label mode), PyMOL chains are label_asym_id (AA, BA, CA)
            # and annotations still use auth chains, so we MAP auth -> label via actual chains
            chain_mapping = {}
            try:
                actual_chains = cmd.get_chains(structure_name)
                if actual_chains:
                    if tool.lower() == 'fr3d':
                        # FR3D: Map numeric chains (1, 2, 3...) to actual chains
                        for idx, actual_chain in enumerate(sorted(actual_chains), 1):
                            chain_mapping[str(idx)] = actual_chain
                    elif tool.lower() in ['rnamotifscan', 'rnamotifscanx']:
                        # RMSX/RMS: Map "0" to first chain (works for both auth and label mode)
                        sorted_chains = sorted(actual_chains)
                        if sorted_chains:
                            chain_mapping['0'] = sorted_chains[0]
                            # If label mode (cif_use_auth=0), map other common auth IDs too
                            if self.cif_use_auth == 0 and len(sorted_chains) > 1:
                                # Map sequential auth IDs to label chains
                                for idx, label_chain in enumerate(sorted_chains):
                                    chain_mapping[str(idx)] = label_chain
                    
                    if self.cif_use_auth == 0:
                        self.logger.debug(f"Label mode chain mapping: {chain_mapping}")
            except Exception as e:
                self.logger.debug(f"Could not get chains from structure: {e}")
            
            # Apply the chain remap directly onto residue objects so the
            # consolidated table (and Jaccard matching against Atlas/Rfam rows)
            # sees real PyMOL chain IDs instead of tool-internal placeholders.
            if chain_mapping:
                for instances in available_motifs.values():
                    for instance in instances:
                        for r in getattr(instance, 'residues', None) or []:
                            if r.chain in chain_mapping:
                                r.chain = chain_mapping[r.chain]

            # Process motifs (same as fetch_motif_data_action)
            motif_summary = {}
            from .utils.parser import SelectionParser
            
            total_count = sum(len(instances) for instances in available_motifs.values())
            self.logger.debug(f"Found {total_count} motifs in {pdb_id} (source: {tool.upper()})")

            # Feed the same structure-local consolidated table used by
            # fetch_motif_data_action, so RMSX/FR3D/RMS rows participate in
            # cross-source family counts and rmv_select alongside Atlas/Rfam.
            if len(self.current_source_names) == 1:
                table = self.annotation_tables.setdefault(
                    pdb_id_upper, ConsolidatedAnnotationTable(self.jaccard_threshold, merge_enabled=False)
                )
                table.jaccard_threshold = self.jaccard_threshold
                table.merge_enabled = False
                table.remove_source(pdb_id_upper, self.current_source_names[0])
                table.add_annotations(
                    pdb_id_upper,
                    self.current_source_names[0],
                    available_motifs,
                    provenance={"provider": tool},
                )
            
            for motif_type, instances in available_motifs.items():
                display_type_upper = motif_type.upper()
                
                # Convert to motif details format
                motif_details = []
                motif_list = []
                
                for instance in instances:
                    if hasattr(instance, 'residues') and instance.residues:
                        # Convert residues to tuple format for display
                        residues_to_use = []
                        for res in instance.residues:
                            if hasattr(res, 'to_tuple'):
                                # ResidueSpec object - convert to tuple
                                residues_to_use.append(res.to_tuple())
                            else:
                                # Already a tuple
                                residues_to_use.append(res)
                        
                        # Apply chain mapping if needed (FR3D)
                        if chain_mapping:
                            remapped = []
                            for nuc, resi, chain in residues_to_use:
                                remapped.append((nuc, resi, chain_mapping.get(chain, chain)))
                            residues_to_use = remapped
                        
                        # CRITICAL: Include metadata (contains aligned_regions for RMSX)
                        from copy import deepcopy
                        instance_metadata = deepcopy(instance.metadata) if hasattr(instance, 'metadata') and instance.metadata else {}
                        
                        motif_details.append({
                            'motif_id': instance.motif_id,
                            'instance_id': instance.instance_id,
                            'residues': residues_to_use,
                            'annotation': instance.annotation,
                            'metadata': instance_metadata,
                            '_source_suffix': source_suffix,
                            '_pdb_id': pdb_id_upper,
                            '_structure_name': structure_name,
                        })
                        
                        
                        # Build motif_list for selection string with remapped chains
                        from .database.user_annotations.converters import MotifInstanceSimple
                        temp_instance = MotifInstanceSimple(
                            instance.motif_id,
                            instance.instance_id,
                            residues_to_use,  # Already converted to tuples above
                            instance.annotation
                        )
                        legacy_entries = temp_instance.to_legacy_format()
                        motif_list.extend(legacy_entries)
                
                # Build main_selection string
                main_motif_sel = None
                if motif_list:
                    all_selections = []
                    for motif in motif_list:
                        chain = motif.get('chain')
                        residues = motif.get('residues')
                        sel = SelectionParser.create_selection_string(chain, residues)
                        if sel:
                            all_selections.append(f"({sel})")
                    
                    if all_selections:
                        combined_sel = " or ".join(all_selections)
                        main_motif_sel = f"(model {structure_name}) and ({combined_sel})"
                
                if motif_details:
                    motif_summary[display_type_upper] = {
                        'object_name': None,
                        'structure_name': structure_name,
                        'pdb_id': pdb_id_upper,
                        'count': len(motif_details),
                        'visible': False,
                        'motif_details': motif_details,
                        'motifs': motif_list,
                        'main_selection': main_motif_sel,
                        'source_suffix': source_suffix,
                    }
                    self.logger.debug(f"Loaded {len(motif_details)} {display_type_upper} motifs")
            
            # Sort motif_details within each type by minimum residue number
            def _get_min_residue(detail):
                residues = detail.get('residues', [])
                if not residues:
                    return float('inf')
                min_resi = float('inf')
                for res in residues:
                    if isinstance(res, tuple) and len(res) >= 2:
                        resi = res[1]
                        if isinstance(resi, int):
                            min_resi = min(min_resi, resi)
                return min_resi if min_resi != float('inf') else float('inf')
            
            for mtype_key, mtype_info in motif_summary.items():
                mtype_info['motif_details'].sort(key=_get_min_residue)
            
            # --- Persist to the residue-based hierarchy cache (single-source loads) ---
            self._persist_single_source_hierarchy(pdb_id_upper, motif_summary)
            
            # Accumulate into existing loaded_motifs (supports cross-PDB and
            # multi-source workflows).  Loading a user-annotation source must NOT
            # wipe motif data previously loaded from other sources/PDBs, otherwise
            # cross-source superimposition (e.g. rmv_super K-TURN 1S72_S7, 1S72_S3)
            # would silently lose instances while the tag registry still
            # advertises them.  Display commands filter to the current source.
            existing_loaded = self.viz_manager.motif_loader.loaded_motifs

            # --- Clean sweep: remove ALL instances matching this PDB+suffix ---
            for key in list(existing_loaded.keys()):
                ex = existing_loaded[key]
                ex['motif_details'] = [
                    d for d in ex.get('motif_details', [])
                    if not (d.get('_pdb_id', ex.get('pdb_id', '')) == pdb_id_upper
                            and d.get('_source_suffix', ex.get('source_suffix', '')) == source_suffix)
                ]
                ex['count'] = len(ex['motif_details'])
                if not ex['motif_details']:
                    del existing_loaded[key]

            for key, new_info in motif_summary.items():
                if key in existing_loaded:
                    ex = existing_loaded[key]
                    new_suffix = new_info.get('source_suffix', '')
                    new_pdb = new_info.get('pdb_id', '')
                    ex['motif_details'] = [
                        d for d in ex['motif_details']
                        if not (d.get('_pdb_id', ex.get('pdb_id', '')) == new_pdb
                                and d.get('_source_suffix', ex.get('source_suffix', '')) == new_suffix)
                    ]
                    ex['motif_details'].extend(new_info['motif_details'])
                    ex['motifs'].extend(new_info.get('motifs', []))
                    ex['count'] = len(ex['motif_details'])
                    new_sel = new_info.get('main_selection', '')
                    if new_sel:
                        if ex.get('main_selection'):
                            ex['main_selection'] = f"{ex['main_selection']} or {new_sel}"
                        else:
                            ex['main_selection'] = new_sel
                    ex['object_name'] = None
                    ex['visible'] = False
                else:
                    existing_loaded[key] = new_info

            # Re-sort after accumulation
            for mtype_key, mtype_info in existing_loaded.items():
                mtype_info['motif_details'].sort(key=_get_min_residue)

            # Register this PDB+source combo for cross-PDB tracking
            if source_suffix:
                self.loaded_sources.add((pdb_id_upper, source_suffix))
            
            # CRITICAL: Also set structure_loader fields so rmv_save can find them
            self.viz_manager.structure_loader.current_structure = structure_name
            self.viz_manager.structure_loader.current_pdb_id = pdb_id_upper
            
            if motif_summary:
                self._report_loaded_families(
                    pdb_id, pdb_id_upper, source_suffix, tool.upper(), motif_summary
                )
            
        except Exception as e:
            self.logger.error(f"Failed to load user annotations for tool={tool}, pdb_id={pdb_id}: {type(e).__name__}: {e}")
            if tool and str(tool).lower() == 'fr3d':
                self.logger.info("FR3D loading failed after execution; check the output file format and selected Python interpreter.")
            import traceback
            traceback.print_exc()
    
    def _list_user_annotations(self):
        """List all available user annotation files."""
        try:
            from pathlib import Path
            plugin_dir = Path(__file__).parent
            user_annotations_dir = plugin_dir / 'database' / 'user_annotations'
            
            print("\n" + "="*60)
            print("Available User Annotation Files")
            print("="*60)
            
            found_any = False
            
            # Check each tool directory
            for tool_dir in user_annotations_dir.iterdir():
                if not tool_dir.is_dir():
                    continue
                
                tool_name = tool_dir.name
                files = list(tool_dir.glob('*.csv')) + list(tool_dir.glob('*.tsv'))
                
                if files:
                    found_any = True
                    print(f"\n{tool_name.upper()}:")
                    for f in files:
                        print(f"  - {f.name}")
            
            if not found_any:
                print("\nNo annotation files found.")
                print("Place files in:")
                print("  - database/user_annotations/fr3d/")
                print("  - database/user_annotations/rnamotifscan/")
            
            print("\n" + "="*60 + "\n")
            
        except Exception as e:
            print(f"Error listing user annotations: {e}")
    
    def switch_database_action(self, database_id):
        """
        Switch to a different database and reload motifs.
        
        Args:
            database_id (str): Database ID to switch to
        """
        try:
            # Check if structure is loaded
            info = self.viz_manager.get_structure_info()
            if not info.get('pdb_id'):
                # Just switch without reloading
                registry = get_registry()
                if registry.set_active_provider(database_id):
                    self.logger.success(f"Switched to database: {database_id}")
                else:
                    self.logger.error(f"Database not found: {database_id}")
                return
            
            # Reload with new database
            motifs = self.viz_manager.reload_with_database(database_id)
            
            if not motifs:
                self.logger.warning(f"No motifs found in {database_id}")
                return
            
            # Update UI state
            self.motif_visibility = {}
            for motif_type, info in motifs.items():
                self.motif_visibility[motif_type] = True
            
            self.logger.success(f"Reloaded with {len(motifs)} motif types from {database_id}")
            
        except Exception as e:
            self.logger.error(f"Failed to switch database: {e}")
    
    def toggle_motif_action(self, motif_type, visible):
        """
        Toggle visibility of a motif type.
        
        Args:
            motif_type (str): Motif type
            visible (bool): Visibility state
        """
        try:
            success = self.viz_manager.motif_loader.toggle_motif_type(motif_type, visible)
            if success:
                self.motif_visibility[motif_type] = visible
                status = "shown" if visible else "hidden"
                self.logger.info(f"Motif {motif_type} {status}")
            else:
                self.logger.warning(f"Could not toggle motif {motif_type}")
        except Exception as e:
            self.logger.error(f"Failed to toggle motif visibility: {e}")
    
    def save_all_motif_images_action(self, representation='cartoon'):
        """
        Save images of all loaded motif instances.
        
        Creates folder structure: plugin_dir/motif_images/pdb_id/MOTIF_TYPE/instance_*_info.png
        
        Args:
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        """
        try:
            success = self.viz_manager.save_all_motif_images(representation=representation)
            if success:
                self.logger.success("All motif images saved successfully")
                self._print_save_location("image", self.loaded_pdb_id)
            else:
                self.logger.error("Failed to save motif images")
        except Exception as e:
            self.logger.error(f"Failed to save motif images: {e}")
    
    def save_motif_type_images_action(self, motif_type, representation='cartoon'):
        """
        Save images for a specific motif type.
        
        Creates folder structure: plugin_dir/motif_images/pdb_id/MOTIF_TYPE/instance_*_info.png
        
        Args:
            motif_type (str): Motif type to save (e.g., 'HL', 'IL')
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        """
        try:
            motif_type = motif_type.upper().strip()
            loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs()
            
            if not loaded_motifs:
                self.logger.error("No motifs loaded")
                return
            
            if motif_type not in loaded_motifs:
                self.logger.error(f"Motif type '{motif_type}' not found")
                self.logger.info(f"Available: {', '.join(sorted(loaded_motifs.keys()))}")
                return
            
            success = self.viz_manager.save_motif_type_images(motif_type, representation=representation)
            if success:
                self.logger.success(f"Saved {motif_type} images successfully")
                self._print_save_location("image", self.loaded_pdb_id)
            else:
                self.logger.error(f"Failed to save {motif_type} images")
        except Exception as e:
            self.logger.error(f"Failed to save motif images: {e}")
    
    def save_motif_instance_by_id_action(self, motif_type, instance_id, representation='cartoon'):
        """
        Save image for a specific motif instance.
        
        Args:
            motif_type (str): Motif type (e.g., 'HL', 'IL')
            instance_id (int): Instance number (1-based, as shown in rmv_summary)
            representation: Display representation ('cartoon', 'sticks', 'spheres', etc.)
                          Default: 'cartoon'
        """
        try:
            motif_type = motif_type.upper().strip()
            loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs()
            
            if not loaded_motifs:
                self.logger.error("No motifs loaded")
                return
            
            if motif_type not in loaded_motifs:
                self.logger.error(f"Motif type '{motif_type}' not found")
                self.logger.info(f"Available: {', '.join(sorted(loaded_motifs.keys()))}")
                return
            
            # Check if instance ID is valid
            motif_instances = loaded_motifs[motif_type]['motif_details']
            if instance_id < 1 or instance_id > len(motif_instances):
                self.logger.error(f"Instance ID {instance_id} out of range (1-{len(motif_instances)})")
                return
            
            success = self.viz_manager.save_motif_instance_by_id(motif_type, instance_id, 
                                                               representation=representation)
            if success:
                self.logger.success(f"Saved {motif_type} instance #{instance_id} successfully")
                self._print_save_location("image", self.loaded_pdb_id)
            else:
                self.logger.error(f"Failed to save {motif_type} instance #{instance_id}")
        except Exception as e:
            self.logger.error(f"Failed to save motif instance: {e}")
    
    def save_current_view_action(self, filename):
        """
        Save the current PyMOL view to high-resolution PNG.
        Preserves exact rotation, angle, and zoom at 2400x1800 / 300 dpi.
        
        Args:
            filename (str): Output filename (e.g., 'my_structure.png')
        """
        try:
            from pathlib import Path
            success = self.viz_manager.save_current_view(filename)
            if success:
                # Show full path
                filepath = Path(filename).resolve()
                self.logger.success(f"Saved current view to: {filepath}")
                self.logger.info(f"  Resolution: 2400x1800 px, 300 dpi")
            else:
                self.logger.error(f"Failed to save current view")
        except Exception as e:
            self.logger.error(f"Failed to save current view: {e}")

    # ------------------------------------------------------------------
    # Structure export (mmCIF) action methods
    # ------------------------------------------------------------------

    def export_all_motif_structures_action(self):
        """Export all loaded motif instances as mmCIF files (original coordinates)."""
        try:
            success = self.viz_manager.export_all_motif_structures()
            if success:
                self.logger.success("All motif structures exported as mmCIF")
                self._print_save_location("structure", self.loaded_pdb_id)
            else:
                self.logger.error("Failed to export motif structures")
        except Exception as e:
            self.logger.error(f"Failed to export motif structures: {e}")

    def export_motif_type_structures_action(self, motif_type):
        """Export all instances of a specific motif type as mmCIF."""
        try:
            motif_type = motif_type.upper().strip()
            loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs()

            if not loaded_motifs:
                self.logger.error("No motifs loaded")
                return

            if motif_type not in loaded_motifs:
                self.logger.error(f"Motif type '{motif_type}' not found")
                self.logger.info(f"Available: {', '.join(sorted(loaded_motifs.keys()))}")
                return

            success = self.viz_manager.export_motif_type_structures(motif_type)
            if success:
                self.logger.success(f"Exported {motif_type} structures as mmCIF")
                self._print_save_location("structure", self.loaded_pdb_id)
            else:
                self.logger.error(f"Failed to export {motif_type} structures")
        except Exception as e:
            self.logger.error(f"Failed to export motif structures: {e}")

    def export_motif_instance_by_id_action(self, motif_type, instance_id):
        """Export a specific motif instance as mmCIF."""
        try:
            motif_type = motif_type.upper().strip()
            loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs()

            if not loaded_motifs:
                self.logger.error("No motifs loaded")
                return

            if motif_type not in loaded_motifs:
                self.logger.error(f"Motif type '{motif_type}' not found")
                self.logger.info(f"Available: {', '.join(sorted(loaded_motifs.keys()))}")
                return

            motif_instances = loaded_motifs[motif_type]['motif_details']
            if instance_id < 1 or instance_id > len(motif_instances):
                self.logger.error(f"Instance ID {instance_id} out of range (1-{len(motif_instances)})")
                return

            success = self.viz_manager.export_motif_instance_structure(motif_type, instance_id)
            if success:
                self.logger.success(f"Exported {motif_type} instance #{instance_id} as mmCIF")
                self._print_save_location("structure", self.loaded_pdb_id)
            else:
                self.logger.error(f"Failed to export {motif_type} instance #{instance_id}")
        except Exception as e:
            self.logger.error(f"Failed to export motif instance: {e}")

    def get_available_motifs(self):
        """
        Get list of available motif types for current PDB.
        
        Returns:
            list: Motif type names
        """
        try:
            pdb_id = self.viz_manager.structure_loader.get_current_pdb_id()
            if not pdb_id:
                return []
            
            motif_types = self.viz_manager.motif_loader.get_available_motif_types(pdb_id)
            return motif_types
        except Exception as e:
            self.logger.error(f"Failed to get motif types: {e}")
            return []
    
    def get_motif_summary(self, pdb_id):
        """
        Get human-readable summary of available motifs for a PDB.
        
        Args:
            pdb_id (str): PDB ID
            
        Returns:
            str: Summary text
        """
        try:
            return self.viz_manager.get_available_motif_summary(pdb_id)
        except Exception as e:
            self.logger.error(f"Failed to get motif summary: {e}")
            return "Error retrieving motif information"
    
    def set_background_color(self, color_name):
        """
        Change the background color of non-motif residues.
        
        Args:
            color_name (str): PyMOL color name (e.g., 'gray80', 'white', 'lightgray')
        """
        try:
            colors.set_background_color(color_name)
            # Recolor the current structure if one is loaded
            current_structure = self.viz_manager.structure_loader.get_current_structure()
            if current_structure:
                cmd.color(color_name, current_structure)
                self.logger.success(f"Background color changed to {color_name}")
            else:
                self.logger.info(f"Background color preference set to {color_name}")
        except Exception as e:
            self.logger.error(f"Failed to change background color: {e}")
    
    def get_motif_info(self, motif_type):
        """
        Get information about a motif type.
        
        Args:
            motif_type (str): Motif type
        
        Returns:
            dict: Motif information
        """
        motif_type_upper = motif_type.upper()
        
        loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs()
        
        if motif_type_upper not in loaded_motifs:
            return {
                'type': motif_type_upper,
                'loaded': False,
                'count': 0,
                'visible': False,
            }
        
        info = loaded_motifs[motif_type_upper]
        
        return {
            'type': motif_type_upper,
            'loaded': True,
            'count': info.get('count', 0),
            'visible': info.get('visible', False),
            'color': colors.get_color_name(motif_type_upper),
            'description': colors.MOTIF_LEGEND.get(motif_type_upper, {}).get('description', ''),
        }
    
    def list_databases(self):
        """
        List all available databases.
        
        Returns:
            list: Database information dictionaries
        """
        return self.viz_manager.get_available_databases()
    
    def print_status(self):
        """Print current status to PyMOL console."""
        info = self.viz_manager.get_structure_info()
        
        print("\n" + "="*60)
        print("RSMViewer - STATUS")
        print("="*60)
        
        # Database info
        databases = self.list_databases()
        print("\nAvailable Databases:")
        for db in databases:
            active_marker = " [ACTIVE]" if db.get('active') else ""
            print(f"  {db['id']:10s} - {db['name']}{active_marker}")
            print(f"              {db['motif_types']} motif types, {db['pdb_count']} PDB structures")
        
        if info['structure']:
            print(f"\nLoaded Structure: {info['structure']}")
            print(f"PDB ID: {info['pdb_id']}")
            print(f"Using database: {info.get('database', 'N/A')}")
        else:
            print("\nNo structure loaded")
            print("\nTo get started:")
            print("  rmv_load <PDB_ID>")
            print("  rmv_load <PDB_ID>, database=rfam")
            return
        
        if info['motifs']:
            print(f"\nLoaded Motifs ({len(info['motifs'])}):")
            for motif_type, data in info['motifs'].items():
                visible_str = "[visible]" if data['visible'] else "[hidden]"
                print(f"  {motif_type:20s} ({data['count']:2d} instances) {visible_str}")
        else:
            print("\nNo motifs loaded for this structure")
        
        print("="*60 + "\n")
    
    def print_sources(self):
        """Print the four canonical public annotation sources."""
        print("\n" + "="*80)
        print("  AVAILABLE DATA SOURCES")
        print("="*80)
        from .database.source_registry import get_source_registry

        for source in get_source_registry().get_all_sources().values():
            print(f"  {source.name:<20} {source.description}")

        print("\n  Use these exact source names (case-insensitive):")
        print("    RNA3DMotifAtlas, Rfam, FR3D, RNAMotifScanX")
        print("\n  Usage:")
        print("    rmv_db <source>[,<source>...]")
        print("\n  Examples:")
        print("    rmv_db RNA3DMotifAtlas")
        print("    rmv_db RNA3DMotifAtlas,Rfam")
        print("\n  rmv_db loads annotations immediately for the active structure.")
        print("="*80 + "\n")
    
    def print_help(self):
        """Print all available commands, grouped in workflow order."""
        print("\n" + "=" * 80)
        print("RSMViewer v2.0.0 - COMMAND REFERENCE")
        print("=" * 80)
        print("Tip: if a command reports a usage error, run rmv_help to see the syntax.")
        print("Annotation sources (case-insensitive): RNA3DMotifAtlas, Rfam, FR3D, RNAMotifScanX")

        print("\n1. LOAD A STRUCTURE")
        print("  rmv_fetch <PDB_ID>[, <ID2> ...]     Fetch one or more structures from the PDB")
        print("  rmv_fetch /path/to/file.cif         Load a local .pdb or .cif file")
        print("  rmv_fetch <ID>, cif_use_auth=0      Use label_asym_id chains (default: auth)")

        print("\n2. LOAD ANNOTATIONS")
        print("  rmv_db                              List the available annotation sources")
        print("  rmv_db <source>[,<source>...]       Load one or more sources for the structure")
        print("                                      e.g. rmv_db RNA3DMotifAtlas,Rfam")
        print("  rmv_refresh                         Re-fetch, bypassing the online cache")

        print("\n3. QUERY AND GROUP MOTIFS")
        print("  rmv_select <motif>, <structures>, <sources>, as <group>")
        print("      Save a Boolean query as a reusable group of motif instances.")
        print("      Source operators: not, and, or  (precedence: not > and > or).")
        print("  rmv_list                            List every loaded motif row")
        print("  rmv_list <FAMILY>                   List every row in a family (e.g. Sarcin-Ricin)")
        print("  rmv_list <MOTIF_ID>                 List one stable motif ID (e.g. 1S72_00016)")
        print("  rmv_list <group>                    List a saved group")
        print("  rmv_combine_groups <group>[, <group> ...], as <group>")
        print("      Merge saved groups into one, keeping each database's labels.")

        print("\n4. VISUALIZE")
        print("  rmv_view <MOTIF_ID|group>[, ...]    Highlight residues on the structure")
        print("  rmv_view <target>, color=<name>     Highlight in a chosen color")
        print("  rmv_view <target>, padding=<N>      Include N neighboring residues")
        print("  rmv_view all                        Highlight every loaded motif")
        print("  rmv_hide <MOTIF_ID|group|all>       Remove a highlight (reset to gray)")
        print("  rmv_create_object <MOTIF_ID|group>  Create selectable PyMOL objects")
        print("  rmv_bg_color <color>                Change the background color")

        print("\n5. COLOR")
        print("  rmv_colors                          List the supported color names")
        print("  rmv_color <motif>, <color>          Set a motif-family color")
        print("  rmv_set_color <group>, <color>      Set a group color (kept through rmv_combine_groups)")

        print("\n6. COMPARE AND EXPORT")
        print("  rmv_super <MOTIF_ID|group>          Superimpose instances (sequence-independent)")
        print("  rmv_align <MOTIF_ID|group>          Superimpose instances (sequence-dependent)")
        print("  rmv_pair <selection>                Inspect base-pair interactions")
        print("  rmv_pair_batch <selection>          Inspect base-pair interactions in batch")
        print("  rmv_save <group|MOTIF_ID|ALL> cif   Export motifs as minimal mmCIF")
        print("  rmv_save current [file.png]         Save the current view as a PNG image")
        print("  rmv_save ALL [representation]       Save an image of every motif (cartoon default)")
        print("      Every save prints its output directory.")

        print("\n7. EXTERNAL PIPELINES (FR3D, RNAMotifScanX)")
        print("  rmv_setup FR3D                      One-shot: install deps and register FR3D")
        print("  rmv_db FR3D                         Run FR3D and load its results")
        print("  rmv_fr3d status|register|run        Inspect or manage the FR3D integration")
        print("  rmv_db RNAMotifScanX                Load RNAMotifScanX results")
        print("  rmv_rmsx status|doctor|test|run     Inspect or run the RNAMotifScanX integration")
        print("  rmv_rmsx_doctor                     Diagnose the RNAMotifScanX runtime")
        print("      Settings live in config/fr3d_config.json and config/rmsx_config.json.")

        print("\n8. SESSION AND DIAGNOSTICS")
        print("  rmv_source info [<source>]          Show the active source configuration")
        print("  rmv_loaded                          List loaded structure + source tags")
        print("  rmv_chains [structure]              Show chain / auth-label diagnostics")
        print("  rmv_debug ON|OFF                    Turn diagnostic messages on or off")
        print("  rmv_reset                           Delete objects and clear all state/caches")
        print("  rmv_help                            Show this reference")

        print("\nQUICK START")
        print("  rmv_fetch 1S72")
        print("  rmv_db RNA3DMotifAtlas,Rfam")
        print("  rmv_select Sarcin-Ricin, 1S72, RNA3DMotifAtlas or Rfam, as group_SR")
        print("  rmv_list group_SR")
        print("  rmv_view group_SR, color=red")
        print("  rmv_create_object group_SR")
        print("  rmv_super group_SR")
        print("  rmv_save group_SR cif")
        print("=" * 80 + "\n")
        return

    def get_available_motifs(self):
        """Get list of available motif types for current PDB+source."""
        filtered = self._get_current_source_motifs()
        return list(filtered.keys()) if filtered else []

    def _current_active_source_ids(self) -> List[int]:
        """Source ID(s) actually selected right now via rmv_db (combine list, or the single active source)."""
        if self.combined_source_ids and len(self.combined_source_ids) >= 2:
            return list(self.combined_source_ids)
        if isinstance(self.current_source_id, int):
            return [self.current_source_id]
        return []

    def _effective_db_ops(self, ops: List[Tuple[str, int]]) -> List[Tuple[str, int]]:
        """Scope a plain 'rmv_summary/rmv_show/rmv_view TYPE' (no explicit
        'db N ...') to the currently active source(s), instead of showing
        every source ever cached for this PDB across past sessions - e.g.
        loading only db 3 this session must not surface leftover db 7 rows
        from an earlier exploratory run. Explicit db-expressions from the
        user are never touched.
        """
        if ops:
            return ops
        active = self._current_active_source_ids()
        if not active:
            return ops
        return [('OR', sid) for sid in active]

    def _register_query_alias(self, alias: str, pdb_id: str, motif_filter: str,
                               ops: List[Tuple[str, int]], rows, source_command: str) -> bool:
        """Freeze *rows* under a globally-unique alias (Applications 5/6).

        Prints an alert and returns False when the alias name is already
        taken - callers must abort instead of silently proceeding.
        """
        from .database.motif_hierarchy_cache import get_hierarchy_cache
        cache = get_hierarchy_cache()
        if cache.alias_exists(alias):
            print(f"\nAlias '{alias}' is already in use. Please choose a different alias.\n")
            return False
        residue_keys = [r_key for _, r_key, _ in rows]
        labels_snapshot = [
            {'residue_key': r_key, 'labels_by_source': labels}
            for _, r_key, labels in rows
        ]
        cache.create_alias(
            alias, pdb_id, residue_keys,
            motif_filter=motif_filter,
            db_expression=_format_db_expression(ops),
            labels_snapshot=labels_snapshot,
            source_command=source_command,
        )
        colors.get_color(alias)  # force-assign a stable, unique color now
        print(f"Alias '{alias}' saved: {len(residue_keys)} instance(s) frozen for {pdb_id} "
              f"(a unique color has been auto-assigned).")
        return True

    def _show_alias(self, alias: str, padding: int = 0) -> bool:
        """Redisplay 'rmv_show <ALIAS>' for a previously saved alias/group.

        Re-enables an already-created PyMOL group, or recreates its instance
        objects from the frozen residue snapshot if they no longer exist.
        Returns False when *alias* isn't a known alias at all.
        """
        from .database.motif_hierarchy_cache import get_hierarchy_cache
        from .utils.parser import SelectionParser

        entry = get_hierarchy_cache().get_alias(alias)
        if not entry:
            return False

        if alias in cmd.get_names('group_objects') or alias in cmd.get_object_list():
            cmd.enable(alias)
            cmd.zoom(alias)
            print(f"\n  Re-displayed existing group '{alias}'.\n")
            return True

        pdb_id = entry['pdb_id']
        structure_name = self.viz_manager.structure_loader.get_current_structure() or self.loaded_pdb
        obj_names = []
        for idx, r_key in enumerate(entry['residue_keys'], 1):
            by_chain: Dict[str, List[int]] = {}
            for token in r_key.split(';'):
                if not token:
                    continue
                chain, _, num_str = token.partition(':')
                try:
                    by_chain.setdefault(chain, []).append(int(num_str))
                except ValueError:
                    continue
            selections = []
            for chain, nums in by_chain.items():
                sel = SelectionParser.create_selection_string(chain, nums, structure_name)
                if sel:
                    selections.append(f"({sel})")
            if not selections:
                continue
            combined_sel = f"(model {structure_name}) and ({' or '.join(selections)})"
            if padding:
                combined_sel = f"byres ({combined_sel} expand {padding})"
            obj_name = f"{alias}_{idx}_{pdb_id}"
            cmd.create(obj_name, combined_sel)
            cmd.show('cartoon', obj_name)
            colors.set_motif_color_in_pymol(cmd, obj_name, alias)
            obj_names.append(obj_name)

        if not obj_names:
            print(f"\n  No residues could be resolved for alias '{alias}' "
                  f"(is {pdb_id} currently loaded?).\n")
            return False

        cmd.group(alias, " ".join(obj_names))
        cmd.zoom(alias)
        print(f"\n  Recreated {len(obj_names)} object(s) for alias '{alias}'.\n")
        return True
    
    def _print_query_context(self, ops: List[Tuple[str, int]], pdb_id: str) -> None:
        """One-line context banner shown before every rmv_summary/rmv_show/
        rmv_view result: which source(s) contributed, which PDB, and that
        the result comes from the persistent hierarchy cache (not a fresh
        API call).
        """
        seen = []
        for _, sid in ops:
            if sid not in seen:
                seen.append(sid)
        src_str = ', '.join(f"{_short_source_name(sid)}({sid})" for sid in seen) or 'all cached sources'
        print(f"Source(s): {src_str}  |  PDB: {pdb_id}  |  Data: hierarchy cache (persistent)")

    def print_source_query_table(self, motif_filter: str, ops: List[Tuple[str, int]],
                                  alias: Optional[str] = None, source_command: str = 'rmv_summary'):
        """Print the residue-based hierarchy table for a motif name and/or 'db N [and/or/not db M]' filter."""
        from .database.motif_hierarchy_cache import get_hierarchy_cache

        ops = self._effective_db_ops(ops)
        pdb_id, rows = self._matching_query_rows(motif_filter, ops)
        if not pdb_id:
            print("\nNo structure loaded. Use 'rmv_fetch <PDB_ID>' first.\n")
            return
        if not rows:
            print(f"\nNo motif instances match this query for {pdb_id}. Run 'rmv_load_motif' first if nothing has been cached yet.\n")
            return

        if alias and not self._register_query_alias(alias, pdb_id, motif_filter, ops, rows, source_command):
            return

        print()
        self._print_query_context(ops, pdb_id)

        # Explicit 'db N' filter -> show exactly those columns (even if '-').
        # No filter (plain 'rmv_summary TYPE') -> show every source that
        # actually contributed a label to at least one matched instance.
        requested_ids = []
        for _, sid in ops:
            if sid not in requested_ids:
                requested_ids.append(sid)
        seen_ids = sorted({sid for _, _, labels in rows for sid in labels})
        column_ids = requested_ids if ops else seen_ids

        # Paper spec: a hierarchical source (e.g. RNA 3D Motif Atlas: generic
        # secondary-structure context, then 3D-geometry subclass) gets one
        # dedicated column per level; a flat source gets a single column.
        # The level count is the max seen for that source across all matched
        # rows, so one instance lacking a refinement doesn't collapse the
        # column for the whole table.
        levels_by_source = {
            sid: max((len(labels_by_source.get(sid, [])) for _, _, labels_by_source in rows), default=1) or 1
            for sid in column_ids
        }

        # Internal custom_id (e.g. 1S720001) stays in the SQLite cache only;
        # the user-facing table shows a plain sequential instance number.
        # Short source names (BGSU, RMSX, ...) keep headers readable - the
        # full descriptive name is still shown by rmv_db/rmv_source info.
        # No PDB column: _matching_query_rows only ever covers the single
        # currently-loaded PDB, so repeating it on every row is redundant -
        # it's already stated once in the context banner above.
        headers = ["NO.", "RESIDUES"]
        for sid in column_ids:
            name = _short_source_name(sid)
            if levels_by_source[sid] >= 2:
                headers.append(f"{name}({sid}).L1")
                headers.append(f"{name}({sid}).L2")
            else:
                headers.append(f"{name}({sid})")
        col_widths = [len(h) for h in headers]

        # Cap the RESIDUES column so a handful of many-chain instances don't
        # stretch the whole table off-screen. Rows over the cap show a
        # truncated cell and their full residue ranges on an indented
        # continuation line right below, so nothing is lost - only wrapped.
        MAX_RESIDUE_COL_WIDTH = 40
        table_rows = []
        for idx, (custom_id, r_key, labels_by_source) in enumerate(rows, 1):
            residues_str = _format_residue_key_ranges(r_key)
            if len(residues_str) > MAX_RESIDUE_COL_WIDTH:
                display_residues = residues_str[:MAX_RESIDUE_COL_WIDTH - 3] + '...'
                full_residues = residues_str
            else:
                display_residues = residues_str
                full_residues = None
            row = [idx, display_residues]
            for sid in column_ids:
                lvls = labels_by_source.get(sid) or []
                if levels_by_source[sid] >= 2:
                    row.append(lvls[0] if len(lvls) >= 1 else '-')
                    row.append(lvls[1] if len(lvls) >= 2 else '-')
                else:
                    row.append(lvls[0] if lvls else '-')
            table_rows.append((row, full_residues))
            for i, cell in enumerate(row):
                col_widths[i] = max(col_widths[i], len(str(cell)))

        # Table name matches the alias when the user set one with 'as ALIAS'.
        print(f"Table: {alias}" if alias else f"Table: {motif_filter or 'ALL'} ({pdb_id})")
        print("  ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers)))
        print("  ".join("-" * w for w in col_widths))
        for row, full_residues in table_rows:
            print("  ".join(str(c).ljust(col_widths[i]) for i, c in enumerate(row)))
            if full_residues:
                print(f"{' ' * col_widths[0]}  -> {full_residues}")
        print(f"\n{len(table_rows)} instance(s) matched.")
        print(f"Hierarchy cache: {get_hierarchy_cache().db_path}\n")

    def _matching_query_rows(self, motif_filter: str, ops: List[Tuple[str, int]]):
        """Shared resolver: residue/Jaccard-clustered rows matching a query for the loaded PDB."""
        from .database.motif_hierarchy_cache import get_hierarchy_cache, residue_key

        pdb_id = self.loaded_pdb_id
        if not pdb_id:
            return pdb_id, []

        hierarchy = get_hierarchy_cache().get_hierarchy_for_pdb(pdb_id)
        if not hierarchy:
            return pdb_id, []

        # Cluster by residue set (exact match, subset, or Jaccard >= threshold)
        # instead of by exact hierarchy key, so near-duplicate residue spans
        # reported by different sources collapse into one instance/row - the
        # same rule residue_merger.py uses to decide "shared" during merge.
        clusters = _cluster_hierarchy_entries(hierarchy, self.jaccard_threshold)

        motif_filter_norm = _normalize_motif_alias_text(motif_filter) if motif_filter.strip() else ''
        rows = []
        for cluster in clusters:
            labels_by_source = cluster['labels_by_source']
            if motif_filter_norm:
                # Forward-only, per-source: a source counts for the predicate
                # only when THAT source labels the cluster as the queried motif
                # (e.g. 'SR' in 'SARCIN-RICIN'). Matching the other direction
                # too would let a short generic label like 'HL' match ANY filter
                # that happens to end in "hl", e.g. 'HAIRPIN LOOP (HL)' wrongly
                # pulling in every GNRA/UNCG/PSEUDOKNOT/etc. instance whose
                # level-1 label is just 'HL'.
                present_ids = {
                    sid for sid, labels in labels_by_source.items()
                    if any(
                        motif_filter_norm in _normalize_motif_alias_text(lbl)
                        for lbl in labels
                    )
                }
                if not present_ids:
                    continue
            else:
                present_ids = set(labels_by_source.keys())
            if not _evaluate_db_expression(ops, present_ids):
                continue
            r_key = residue_key(list(cluster['residues']))
            rows.append((cluster['custom_id'], r_key, labels_by_source))
        rows.sort(key=lambda r: r[0])
        return pdb_id, rows

    def show_source_query(self, motif_filter: str, ops: List[Tuple[str, int]], padding: int = 0,
                           alias: Optional[str] = None, source_command: str = 'rmv_show'):
        """PyMOL rendering for a 'rmv_show [MOTIF] db N [and/or/not db M ...]' query."""
        from .utils.parser import SelectionParser

        ops = self._effective_db_ops(ops)
        pdb_id, rows = self._matching_query_rows(motif_filter, ops)
        if not pdb_id:
            print("\nNo structure loaded. Use 'rmv_fetch <PDB_ID>' first.\n")
            return
        if not rows:
            print(f"\nNo motif instances match this query for {pdb_id}.\n")
            return

        if alias and not self._register_query_alias(alias, pdb_id, motif_filter, ops, rows, source_command):
            return

        print()
        self._print_query_context(ops, pdb_id)
        # Object/group name matches the alias when the user set one with 'as ALIAS'.
        print(f"Object group: {alias}" if alias else f"Objects: {motif_filter or 'ALL'} ({pdb_id})")
        structure_name = self.viz_manager.structure_loader.get_current_structure() or self.loaded_pdb
        created = 0
        obj_names = []
        # Object names use the sequential instance number, never the internal
        # custom_id, so the PyMOL object panel only ever shows 1, 2, 3...
        # An explicit alias becomes the object name base (App 6: "alias will
        # be set as object name"), taking priority over the motif filter text.
        name_base = alias or (re.sub(r'[^A-Za-z0-9_]', '_', motif_filter.strip()) if motif_filter.strip() else None)
        for idx, (custom_id, r_key, labels_by_source) in enumerate(rows, 1):
            by_chain: Dict[str, List[int]] = {}
            for token in r_key.split(';'):
                if not token:
                    continue
                chain, _, num_str = token.partition(':')
                try:
                    by_chain.setdefault(chain, []).append(int(num_str))
                except ValueError:
                    continue
            selections = []
            for chain, nums in by_chain.items():
                sel = SelectionParser.create_selection_string(chain, nums, structure_name)
                if sel:
                    selections.append(f"({sel})")
            if not selections:
                continue
            combined_sel = f"(model {structure_name}) and ({' or '.join(selections)})"
            if padding:
                combined_sel = f"byres ({combined_sel} expand {padding})"
            # Most specific (deepest) label per source, e.g. K-TURN over IL.
            most_specific = next((lvls[-1] for lvls in labels_by_source.values() if lvls), 'MOTIF')
            base = name_base or re.sub(r'[^A-Za-z0-9_]', '_', most_specific)
            obj_name = f"{base}_{idx}_{pdb_id}"
            cmd.create(obj_name, combined_sel)
            cmd.show('cartoon', obj_name)
            # An alias gets its own distinct color for the whole group
            # (App 6: novel vs known coloring); otherwise color per label.
            colors.set_motif_color_in_pymol(cmd, obj_name, alias or most_specific)
            obj_names.append(obj_name)
            created += 1
            # Show every hierarchy level per source (context -> subclass),
            # matching the paper's per-level column design.
            label_str = ", ".join(f"{' > '.join(v)}({k})" for k, v in labels_by_source.items())
            print(f"  Created instance #{idx}  [{label_str}]  {_format_residue_key_ranges(r_key)}")

        if alias and obj_names:
            cmd.group(alias, " ".join(obj_names))
            print(f"Grouped {len(obj_names)} object(s) under '{alias}'.")

        print(f"\n{created} object(s) created for this query.\n")

    def view_source_query(self, motif_filter: str, ops: List[Tuple[str, int]],
                          alias: Optional[str] = None, color_override: Optional[str] = None,
                          source_command: str = 'rmv_view'):
        """In-place structure coloring for a 'rmv_view [MOTIF] db N [and/or/not db M ...]' query."""
        from .utils.parser import SelectionParser

        ops = self._effective_db_ops(ops)
        pdb_id, rows = self._matching_query_rows(motif_filter, ops)
        if not pdb_id:
            print("\nNo structure loaded. Use 'rmv_fetch <PDB_ID>' first.\n")
            return
        if not rows:
            print(f"\nNo motif instances match this query for {pdb_id}.\n")
            return

        if alias and not self._register_query_alias(alias, pdb_id, motif_filter, ops, rows, source_command):
            return

        print()
        self._print_query_context(ops, pdb_id)
        # Highlight name matches the alias when the user set one with 'as ALIAS'.
        print(f"Highlight group: {alias}" if alias else f"Highlight: {motif_filter or 'ALL'} ({pdb_id})")
        structure_name = self.viz_manager.structure_loader.get_current_structure() or self.loaded_pdb
        # Gray out the base structure so the highlighted motifs stand out.
        self._gray_out_base_structures({structure_name})
        colored = 0
        for _custom_id, r_key, labels_by_source in rows:
            by_chain: Dict[str, List[int]] = {}
            for token in r_key.split(';'):
                if not token:
                    continue
                chain, _, num_str = token.partition(':')
                try:
                    by_chain.setdefault(chain, []).append(int(num_str))
                except ValueError:
                    continue
            selections = []
            for chain, nums in by_chain.items():
                sel = SelectionParser.create_selection_string(chain, nums, structure_name)
                if sel:
                    selections.append(f"({sel})")
            if not selections:
                continue
            instance_sel = f"(model {structure_name}) and ({' or '.join(selections)})"
            most_specific = next((lvls[-1] for lvls in labels_by_source.values() if lvls), 'MOTIF')
            if color_override:
                try:
                    cmd.color(color_override, instance_sel)
                except Exception:
                    colors.set_motif_color_in_pymol(cmd, instance_sel, alias or most_specific)
            else:
                # An alias gets one distinct color for the whole group;
                # otherwise fall back to the per-label motif color.
                colors.set_motif_color_in_pymol(cmd, instance_sel, alias or most_specific)
            colored += 1

        print(f"\n{colored} instance(s) highlighted on {structure_name} for this query.\n")

    def print_motif_summary(self):
        """Print detailed motif summary table to console."""
        info = self.viz_manager.get_structure_info()
        
        # If no info from viz_manager, check if we loaded via rmv_fetch
        if not info.get('pdb_id') and not self.loaded_pdb_id:
            print("\nNo structure loaded. Use 'rmv_fetch <PDB_ID>' or 'rmv_load <PDB_ID>' first.\n")
            return
        
        # Use viz_manager info if available, otherwise use our stored data
        pdb_id = info.get('pdb_id') or self.loaded_pdb_id

        # Filter to current PDB+source so we only display the active view
        motifs = self._get_current_source_motifs()
        
        if not motifs:
            print(f"\nNo motifs loaded for {pdb_id}.\n")
            return
        
        # Determine the database name to display based on current source mode
        source_names = {
            'atlas': 'RNA 3D Motif Atlas (Local)',
            'rfam': 'Rfam (Local)',
            'bgsu_api': 'BGSU RNA 3D Hub (API)',
            'rfam_api': 'Rfam (API)',
            'fr3d': 'FR3D User Annotations',
            'rnamotifscan': 'RNAMotifScan User Annotations',
            'rnamotifscanx': 'RNAMotifScanX User Annotations',
        }
        
        database_id = "Unknown"
        
        # Determine database name based on current source mode
        if self.current_source_mode == 'user':
            if self.current_user_tool:
                tool_display_names = {
                    'fr3d': 'FR3D User Annotations',
                    'rnamotifscan': 'RMS (RNAMotifScan) User Annotations',
                    'rnamotifscanx': 'RMSX (RNAMotifScanX) User Annotations',
                    'nobias': 'NoBIAS User Annotations',
                }
                database_id = tool_display_names.get(self.current_user_tool, f"{self.current_user_tool} User Annotations")
            else:
                database_id = "User Annotations"
        
        elif self.current_source_mode == 'local':
            if self.current_local_source:
                database_id = source_names.get(self.current_local_source, f"{self.current_local_source} (Local)")
            else:
                database_id = 'Local Databases (Atlas + Rfam)'
        
        elif self.current_source_mode == 'web':
            if self.current_web_source:
                database_id = source_names.get(self.current_web_source, f"{self.current_web_source} (API)")
            else:
                database_id = 'Online APIs'
        
        elif self.current_source_mode == 'auto':
            database_id = 'Auto-selected (Local First -> API)'
        
        elif self.current_source_mode == 'combine':
            database_id = 'Combined (Multiple Sources)'
        
        # Fallback: map provider_id if available
        if database_id == "Unknown":
            provider_id = info.get('database_id')
            if provider_id:
                database_id = source_names.get(provider_id, provider_id)
        
        # Use the visualization manager's summary printer
        self.viz_manager._print_motif_summary_table(pdb_id, motifs, database_id)
    
    def show_motif_summary_for_type(self, motif_type: str, alias: Optional[str] = None):
        """Print residue-based hierarchy table for a motif name (e.g., 'HL', 'K-TURN').

        Matches any cached instance where at least one source's label
        contains this name, regardless of exact label spelling per source.
        """
        motif_arg = motif_type.upper().strip()
        self.print_source_query_table(motif_arg, [], alias=alias)
        print(f"  rmv_summary {motif_arg} <NO>      Show details of specific instance")
        print(f"  rmv_show {motif_arg}              Render & create objects for {motif_arg}")
        print(f"  rmv_super {motif_arg}             Superimpose all {motif_arg} instances\n")

    @staticmethod
    def _instance_from_residue_set(residue_set, family, source_label, instance_id, source_labels=()):
        """Build a lightweight MotifInstance from a consolidated row residue set."""
        from .database.base_provider import MotifInstance, ResidueSpec
        residues = [
            ResidueSpec(
                chain=str(chain),
                residue_number=int(number),
                insertion_code=str(insertion or ''),
                model=int(model),
            )
            for chain, number, insertion, model in residue_set
        ]
        # _source_family_labels keeps each database's original wording so the
        # displayed value is not collapsed to the normalized query family.
        family_labels = {source_label: list(source_labels)} if source_label and source_labels else {}
        return MotifInstance(
            instance_id=str(instance_id),
            motif_id=str(family),
            pdb_id='',
            residues=residues,
            annotation='',
            metadata={
                '_source_label': source_label,
                '_source_family_labels': family_labels,
            },
        )

    def _merge_source_instances(self, structure_id, ordered_sources, ordered_labels, group_name):
        """Residue-merge per-source instance dicts into consolidated AnnotationRows.

        Uses :class:`ResidueMerger` so the merge is residue-set based and keeps
        the larger set: strict subset/superset removal (100% containment) plus
        Jaccard >= threshold. This is the professor's rule applied at
        rmv_select (same family) and rmv_combine_groups (across families) time.
        """
        from .database.residue_merger import ResidueMerger
        from .database.consolidated_table import normalize_residue_set, AnnotationRow

        merger = ResidueMerger(jaccard_threshold=self.jaccard_threshold)
        # pdb_id=None so the transient select/combine merge does not write to
        # the shared hierarchy cache; per-source labels still propagate via
        # each instance's _also_found_in metadata.
        merged = merger.merge_sources(ordered_sources, ordered_labels, pdb_id=None)

        rows = []
        index = 0
        for motif_type, instances in merged.items():
            for instance in instances:
                residue_set = normalize_residue_set(instance.residues)
                if not residue_set:
                    continue
                index += 1
                meta = instance.metadata or {}
                labels = [meta.get('_source_label', '')] + list(meta.get('_also_found_in', []))
                family_labels = meta.get('_source_family_labels', {}) or {}
                source_annotations = {}
                for label in labels:
                    if not label:
                        continue
                    original = tuple(family_labels.get(label, ()))
                    source_annotations.setdefault(label, original or (str(motif_type),))
                rows.append(
                    AnnotationRow(
                        motif_id=f"{group_name}_{index:03d}",
                        structure_id=structure_id,
                        residue_set=residue_set,
                        source_annotations=source_annotations,
                    )
                )
        return rows

    def _merge_selection_rows(self, structure_id, entries, family, group_name):
        """Merge matched raw rows of ONE family into consolidated group rows.

        ``entries`` is a list of ``(AnnotationRow, matching_sources)`` for a
        single structure; each contributing source becomes a merge input so the
        surviving instances carry every source that labelled them.
        """
        from collections import defaultdict

        per_source: Dict[str, Dict[str, list]] = defaultdict(dict)
        discovered: List[str] = []
        for row, sources in entries:
            for source in sources:
                if source not in discovered:
                    discovered.append(source)
                original_labels = list(row.source_annotations.get(source, ()))
                for label in row.source_hierarchy.get(source, ()):
                    if label and label not in original_labels:
                        original_labels.append(label)
                instance = self._instance_from_residue_set(
                    row.residue_set, family, source, f"{source}:{row.motif_id}",
                    source_labels=original_labels,
                )
                per_source[source].setdefault(family, []).append(instance)

        ordered_labels = [name for name in self.current_source_names if name in per_source]
        for source in discovered:
            if source not in ordered_labels:
                ordered_labels.append(source)
        ordered_sources = [per_source[name] for name in ordered_labels]
        return self._merge_source_instances(structure_id, ordered_sources, ordered_labels, group_name)

    def _combine_group_rows(self, structure_id, ordered_group_rows, new_alias):
        """Merge rows from several selected groups while keeping each database's labels.

        ``ordered_group_rows`` is a list of ``(group_name, [AnnotationRow])`` in
        precedence order. Overlapping instances are merged (strict containment
        keeps the larger residue set; otherwise Jaccard >= threshold keeps the
        earlier group's residues), and every contributing database's labels are
        unioned onto the surviving row - so a combined row can show, e.g.,
        RNAMotifScanX=Kink-turn and FR3D=Sarcin-Ricin side by side, including
        when the sources disagree about the family.

        Returns a list of ``(AnnotationRow, contributing_group_names)``.
        """
        from .database.consolidated_table import residue_jaccard, AnnotationRow

        clusters: List[dict] = []
        for group_name, rows in ordered_group_rows:
            for row in rows:
                rset = row.residue_set
                if not rset:
                    continue
                incoming = {
                    source: set(labels)
                    for source, labels in row.source_annotations.items()
                    if source
                }
                incoming_set = set(rset)
                target = None
                for cluster in clusters:
                    cluster_set = set(cluster["rset"])
                    contained = incoming_set <= cluster_set or cluster_set <= incoming_set
                    if contained or residue_jaccard(cluster["rset"], rset) >= self.jaccard_threshold:
                        target = cluster
                        break
                if target is None:
                    clusters.append({
                        "rset": rset,
                        "sources": incoming,
                        "groups": [group_name],
                    })
                    continue
                for source, labels in incoming.items():
                    target["sources"].setdefault(source, set()).update(labels)
                if group_name not in target["groups"]:
                    target["groups"].append(group_name)
                # Keep the strictly larger (superset) residue set as representative.
                if incoming_set > set(target["rset"]):
                    target["rset"] = rset

        merged: List[Tuple] = []
        for cluster in clusters:
            source_annotations = {
                source: tuple(sorted(labels))
                for source, labels in cluster["sources"].items()
                if source
            }
            merged.append((
                AnnotationRow(
                    motif_id=f"{new_alias}_pending",
                    structure_id=structure_id,
                    residue_set=cluster["rset"],
                    source_annotations=source_annotations,
                ),
                list(cluster["groups"]),
            ))
        return merged

    def _rows_for_target(self, target: str):
        """Return AnnotationRows for a stable motif ID or a saved group.

        Groups produced by rmv_select/rmv_combine_groups store their merged rows
        directly; legacy groups and bare motif IDs resolve against the raw
        load-time tables. A bare id may also be a single merged row inside a
        saved group (e.g. ``group_SR_001``).
        """
        group = self.query_groups.get(target)
        if group is not None:
            if group.get("rows"):
                return list(group["rows"])
            motif_ids = set(group.get("motif_ids", []))
            return [
                row
                for table in self.annotation_tables.values()
                for row in table.rows
                if row.motif_id in motif_ids
            ]
        # Bare id: raw load-time tables first...
        matches = [
            row
            for table in self.annotation_tables.values()
            for row in table.rows
            if row.motif_id == target
        ]
        if matches:
            return matches
        # ...then an individual merged row inside a saved group.
        for saved in self.query_groups.values():
            for row in saved.get("rows", []):
                if row.motif_id == target:
                    return [row]
        return []

    def select_annotation_query(self, query_text: str) -> None:
        """Save a consolidated-table query as a motif-ID result snapshot."""
        try:
            query = parse_selection_query(query_text)
        except (QuerySyntaxError, ValueError) as exc:
            self.command_error(str(exc))
            return

        if query.group in self.query_groups:
            self.logger.error(f"Group '{query.group}' already exists. Choose another name.")
            return

        from .database.motif_aliases import labels_match_motif

        def _row_labels(row) -> List[str]:
            return [
                label for values in row.source_hierarchy.values() for label in values
            ] + [
                label for values in row.source_annotations.values() for label in values
            ]

        def _fr3d_labels(row) -> List[str]:
            # FR3D motif names embed the family loosely (e.g. query file names),
            # so they are matched by keyword rather than exact canonicalization.
            return list(row.source_hierarchy.get("FR3D", ())) + list(
                row.source_annotations.get("FR3D", ())
            )

        def _sources_matching_motif(row) -> Tuple[str, ...]:
            # A source is "true" for the predicate only when THAT source labels
            # this row as the queried motif, so benchmarking TP/FP/FN stay exact.
            matching = []
            for source in row.source_annotations:
                source_labels = list(row.source_annotations.get(source, ())) + list(
                    row.source_hierarchy.get(source, ())
                )
                if source == "FR3D":
                    if labels_match_motif(query.motif, (), source_labels):
                        matching.append(source)
                elif labels_match_motif(query.motif, source_labels, ()):
                    matching.append(source)
            return tuple(matching)

        family = query.motif or "MOTIF"

        # Step 1: gather every raw row of this family, per structure, with the
        # sources that label it. NO predicate filtering yet - the cross-source
        # support needed by 'A and B' (TP), 'not A and B' (FP) and
        # 'A and not B' (FN) only exists after the residue merge unions the
        # per-source labels onto one instance.
        family_by_structure: Dict[str, List] = {}
        for structure_id in sorted(self.annotation_tables):
            if not query.matches_structure(structure_id):
                continue
            for row in self.annotation_tables[structure_id].rows:
                if query.motif:
                    row_sources = _sources_matching_motif(row)
                else:
                    row_sources = tuple(row.source_annotations)
                if not row_sources:
                    continue
                family_by_structure.setdefault(structure_id, []).append((row, row_sources))

        if not family_by_structure:
            self.logger.error(f"No motifs matched '{query.motif}'.")
            return

        # Step 2: containment + Jaccard merge over this family only (across all
        # loaded sources), so overlapping instances collapse into one row that
        # carries every source which annotated it.
        merged_all = []
        for structure_id, entries in family_by_structure.items():
            merged_all.extend(
                self._merge_selection_rows(structure_id, entries, family, query.group)
            )

        # Step 3: apply the source predicate to each merged instance using the
        # set of sources that actually contributed to it. This is where the
        # benchmark groups split: TP (A and B), FP (not A and B), FN (A and not B).
        final_rows = [
            row
            for row in merged_all
            if query.sources.matches(tuple(row.source_annotations.keys()))
        ]

        # Step 4: assign unique sequential IDs across all structures.
        for index, row in enumerate(final_rows, 1):
            row.motif_id = f"{query.group}_{index:03d}"

        if not final_rows:
            self.logger.error(
                f"'{query.motif}' motifs exist, but none satisfy the requested source "
                f"condition. Run 'rmv_list {query.motif}' to inspect per-source support."
            )
            return

        self.query_groups[query.group] = {
            "rows": final_rows,
            "motif_ids": [row.motif_id for row in final_rows],
            "query": query.text,
            "motif": query.motif,
            "structures": query.structures,
        }
        self.logger.success(
            f"Saved group '{query.group}' with {len(final_rows)} merged motif instance(s)."
        )

    def list_annotation_results(self, target: str = "") -> None:
        """List loaded motif rows, a stable motif ID, or a saved group."""
        target = target.strip()
        motif_ids = None
        group = None
        family_label = None
        if target in self.query_groups:
            group = self.query_groups[target]
            motif_ids = set(group["motif_ids"])
        elif target:
            known_id = any(
                self.annotation_tables[structure_id].get(target) is not None
                for structure_id in self.annotation_tables
            )
            if known_id:
                motif_ids = {target}
            else:
                # Accept a motif-family name (e.g. SARCIN-RICIN, SR): list every
                # row any source labels as that family - the same union the
                # family count after rmv_db reports.
                from .database.motif_aliases import labels_match_motif, canonical_motif
                family_label = canonical_motif(target) or None
                family_ids = set()
                if family_label:
                    for structure_id in self.annotation_tables:
                        for row in self.annotation_tables[structure_id].rows:
                            row_labels = [
                                lbl for values in row.source_hierarchy.values() for lbl in values
                            ] + [
                                lbl for values in row.source_annotations.values() for lbl in values
                            ]
                            fr3d_labels = list(row.source_hierarchy.get("FR3D", ())) + list(
                                row.source_annotations.get("FR3D", ())
                            )
                            if labels_match_motif(family_label, row_labels, fr3d_labels):
                                family_ids.add(row.motif_id)
                if family_ids:
                    motif_ids = family_ids
                else:
                    motif_ids = {target}
                    family_label = None

        if group is not None and group.get("rows"):
            rows = list(group["rows"])
        else:
            rows = [
                row
                for structure_id in sorted(self.annotation_tables)
                for row in self.annotation_tables[structure_id].rows
                if motif_ids is None or row.motif_id in motif_ids
            ]
            if not rows and motif_ids:
                # Individual merged row inside a saved group (e.g. group_SR_001).
                rows = [
                    row
                    for saved in self.query_groups.values()
                    for row in saved.get("rows", [])
                    if row.motif_id in motif_ids
                ]
        if not rows:
            if target:
                print(
                    f"No motif rows found for '{target}'. Use a saved group name, a "
                    f"stable motif ID (e.g. 1S72_00016), or a motif family name "
                    f"(e.g. SARCIN-RICIN). Run rmv_list with no argument to list all rows."
                )
            else:
                print("No motif rows found.")
            return

        source_names = [
            source for source in self.current_source_names
            if any(source in row.source_annotations for row in rows)
        ]
        for source in sorted({source for row in rows for source in row.source_annotations}):
            if source not in source_names:
                source_names.append(source)
        structures = sorted({row.structure_id for row in rows})
        cache_path = Path(__file__).parent / "database" / "motif_hierarchy_cache.sqlite3"

        if group:
            print("\nRSMViewer selected motif group")
            print(f"Group:      {target}")
            print(f"Query:      {group['query']}")
            motif_line = group.get("motif") or group['query'].split(',', 1)[0].strip()
            print(f"Motif:      {motif_line}")
        elif family_label:
            print("\nRSMViewer motif family listing")
            print(f"Motif:      {family_label}")
        print(f"Structures: {', '.join(structures)}")
        print(f"Sources:    {', '.join(source_names) if source_names else '(none)'}")
        if group or family_label:
            print(f"Members:    {len(rows)}")
        print(f"SQLite:     {cache_path}")
        print()

        # Per-source hierarchy depth (context -> subclass); usually 1.
        source_levels = {}
        for source in source_names:
            source_levels[source] = max(
                (len(row.source_hierarchy.get(source, row.source_annotations.get(source, ())))
                 for row in rows),
                default=0,
            )

        headers = ["MOTIF_ID", "RESIDUES"]
        for source in source_names:
            levels = source_levels[source]
            for level in range(levels or 1):
                headers.append(f"{source}.L{level + 1}" if levels > 1 else source)

        table_rows = []
        for row in rows:
            raw_key = ";".join(
                f"{chain}:{number}"
                for chain, number, _insertion, _model in row.residue_set
            )
            residue_text = _format_residue_key_ranges(raw_key)
            cells = [row.motif_id, residue_text]
            for source in source_names:
                levels = source_levels[source] or 1
                values = list(
                    row.source_hierarchy.get(source, row.source_annotations.get(source, ()))
                )
                cells.extend(values + ["-"] * max(0, levels - len(values)))
            table_rows.append(cells)

        # Fixed-width columns so each source's annotation lines up under its header.
        widths = [len(header) for header in headers]
        for cells in table_rows:
            for index, cell in enumerate(cells):
                widths[index] = max(widths[index], len(cell))

        def _fmt_row(cells):
            return "  " + "   ".join(
                cell.ljust(widths[index]) for index, cell in enumerate(cells)
            )

        print(_fmt_row(headers))
        print("  " + "   ".join("-" * widths[index] for index in range(len(headers))))
        for cells in table_rows:
            print(_fmt_row(cells))

    def _group_member_color_key(self, group_name: str, motif_id: str) -> str:
        """Resolve the color key for one member of a (possibly combined) group.

        An explicit color on the group itself paints every member uniformly;
        otherwise a combined group colors each member by the source group it
        came from, so distinct per-source colors survive rmv_combine_groups.
        """
        if not group_name:
            return motif_id
        if colors.has_custom_color(group_name):
            return group_name
        origin = self.query_groups.get(group_name, {}).get("member_colors", {}).get(motif_id)
        return origin or group_name

    def _print_save_location(self, kind: str, pdb_id: str = "") -> None:
        """Print the absolute output directory for a save/export command."""
        plugin_dir = Path(__file__).parent.parent
        base = plugin_dir / ("motif_images" if kind == "image" else "motif_structures")
        if pdb_id:
            base = base / str(pdb_id).lower()
        self.logger.info(f"  Saved to: {base.resolve()}")

    def create_annotation_objects(self, target: str) -> None:
        """Create selectable PyMOL objects for a motif ID or saved group."""
        target = target.strip()
        if target in self.query_groups:
            group_name = target
        else:
            group_name = ""

        rows = self._rows_for_target(target)
        if not rows:
            self.logger.error(f"No consolidated motif found for '{target}'.")
            return

        object_names = []
        skipped = 0
        from .utils.parser import SelectionParser

        for row in sorted(rows, key=lambda item: item.motif_id):
            structure_name = self.loaded_structures.get(row.structure_id, row.structure_id)
            by_chain: Dict[str, List[int]] = {}
            for chain, number, _insertion, _model in row.residue_set:
                by_chain.setdefault(chain, []).append(number)

            auth_parts = []
            segi_parts = []
            for chain, numbers in sorted(by_chain.items()):
                ordered = sorted(set(numbers))
                auth_sel = SelectionParser.create_selection_string(chain, ordered, structure_name)
                if auth_sel:
                    auth_parts.append(f"({auth_sel})")
                segi_sel = SelectionParser.create_selection_string(
                    chain, ordered, structure_name, use_segi=True)
                if segi_sel:
                    segi_parts.append(f"({segi_sel})")
            if not auth_parts:
                skipped += 1
                continue

            object_name = f"motif_{row.motif_id}"
            if object_name in cmd.get_object_list():
                cmd.delete(object_name)

            selection = (
                f"(model {structure_name}) and polymer.nucleic and "
                f"({' or '.join(auth_parts)})"
            )
            cmd.create(object_name, selection)

            # Fall back to segi (label_asym_id) when the auth-chain selection
            # matched no atoms, e.g. structures loaded under label chains.
            if cmd.count_atoms(object_name) == 0 and segi_parts:
                segi_selection = (
                    f"(model {structure_name}) and polymer.nucleic and "
                    f"({' or '.join(segi_parts)})"
                )
                if cmd.count_atoms(segi_selection) > 0:
                    cmd.delete(object_name)
                    cmd.create(object_name, segi_selection)

            if cmd.count_atoms(object_name) == 0:
                cmd.delete(object_name)
                skipped += 1
                self.logger.warning(
                    f"Skipped {object_name}: no atoms matched in {row.structure_id}.")
                continue

            cmd.hide("everything", object_name)
            cmd.show("cartoon", object_name)
            cmd.set("cartoon_nucleic_acid_mode", 4, object_name, quiet=1)
            cmd.set("cartoon_tube_radius", 0.4, object_name, quiet=1)
            colors.set_motif_color_in_pymol(
                cmd, object_name, self._group_member_color_key(group_name, row.motif_id))
            object_names.append(object_name)
            print(f"Created {object_name} from {row.structure_id} ({cmd.count_atoms(object_name)} atoms)")

        if not object_names:
            self.logger.error(
                f"No objects were created for '{target}' — the annotated residues "
                f"did not match any atoms in the loaded structure(s).")
            return

        if group_name:
            cmd.group(group_name, " ".join(object_names))
            print(f"Grouped {len(object_names)} object(s) under '{group_name}'.")
        if skipped:
            self.logger.warning(f"{skipped} motif(s) skipped (no matching atoms).")

    def _gray_out_base_structures(self, structure_names) -> None:
        """Color the given base structures gray80 so highlighted motifs stand out."""
        live = set(cmd.get_object_list())
        for name in {n for n in structure_names if n}:
            if name in live:
                cmd.color("gray80", f"model {name}")

    def view_annotation_results(
        self, target: str, color_override: Optional[str] = None,
        hide: bool = False, padding: int = 0, gray_base: bool = True
    ) -> None:
        """Highlight a stable motif ID or saved group on its parent structure."""
        if target in self.query_groups:
            color_name = color_override or colors.get_color(target)
        else:
            color_name = color_override or colors.get_color("MOTIF")

        rows = self._rows_for_target(target)
        if not rows:
            self.logger.error(f"No consolidated motif found for '{target}'.")
            return

        from .utils.parser import SelectionParser

        # Gray out the base structure(s) so the highlighted motifs stand out.
        if gray_base and not hide:
            self._gray_out_base_structures(
                {self.loaded_structures.get(row.structure_id, row.structure_id)
                 for row in rows})

        for row in rows:
            structure_name = self.loaded_structures.get(row.structure_id, row.structure_id)
            by_chain: Dict[str, List[int]] = {}
            for chain, number, _insertion, _model in row.residue_set:
                by_chain.setdefault(chain, []).append(number)
            selections = [
                SelectionParser.create_selection_string(chain, sorted(numbers), structure_name)
                for chain, numbers in sorted(by_chain.items())
            ]
            selections = [selection for selection in selections if selection]
            if not selections:
                continue
            selection = f"(model {structure_name}) and ({' or '.join(selections)})"
            if padding:
                selection = f"byres ({selection} expand {padding})"
            if hide:
                cmd.color("gray80", selection)
            elif color_override:
                cmd.color(color_override, selection)
            elif target in self.query_groups:
                colors.set_motif_color_in_pymol(
                    cmd, selection, self._group_member_color_key(target, row.motif_id))
            else:
                colors.set_motif_color_in_pymol(cmd, selection, target)
        self.logger.info(f"{'Reset' if hide else 'Highlighted'} motif target '{target}'.")

    def export_annotation_rows(self, target: str) -> int:
        """Export a stable motif ID or saved group as minimal mmCIF files."""
        target = target.strip()
        rows = self._rows_for_target(target)
        if not rows:
            self.logger.error(f"No consolidated motif found for '{target}'.")
            return 0

        from .structure_exporter import MotifStructureExporter

        exporter = MotifStructureExporter(cmd)
        saved = 0
        output_dirs = set()
        for row in sorted(rows, key=lambda item: item.motif_id):
            cif_path = None
            if row.structure_id == self.loaded_pdb_id and getattr(self, "loaded_structure_path", ""):
                candidate = Path(self.loaded_structure_path)
                if candidate.suffix.lower() in (".cif", ".mmcif") and candidate.exists():
                    cif_path = str(candidate)
            cif_path = cif_path or exporter._find_cif_file(row.structure_id)
            if not cif_path:
                self.logger.warning(f"Original CIF not found for {row.structure_id}; skipped {row.motif_id}.")
                continue
            folder = exporter.create_folder_hierarchy(row.structure_id)
            output_dirs.add(str(Path(folder).resolve()))
            folder = exporter.create_motif_type_folder(folder, row.motif_id)
            details = {
                "residues": [("", number, chain) for chain, number, _insertion, _model in row.residue_set]
            }
            if exporter.export_instance(folder, 1, row.motif_id, details, cif_path, row.structure_id):
                saved += 1
        self.logger.info(f"Exported {saved} stable motif mmCIF file(s) for '{target}'.")
        for directory in sorted(output_dirs):
            self.logger.info(f"  Saved to: {directory}")
        return saved
    
    def _print_source_attribution_report(self, motif_type: str, motif_details: list):
        """Print unique merged-instance IDs grouped by source.

        Only prints when instances carry _source_label metadata (combine mode).
        Three categories: unique-to-source-A, unique-to-source-B, and shared
        (instances where the annotation mergingr detected overlap from both sources).
        Also shows within-source deduplication counts when available.
        """
        # Collect per-instance source info
        source_only = {}   # source_label -> [instance_numbers]  (unique to that source)
        shared_ids = []    # instance numbers that were found in multiple sources
        shared_labels = {} # idx -> combined label string

        for idx, detail in enumerate(motif_details, 1):
            meta = detail.get('metadata', {})
            label = meta.get('_source_label', '')
            also = meta.get('_also_found_in', [])
            if not label:
                continue
            if also:
                shared_ids.append(idx)
                shared_labels[idx] = ' + '.join([label] + also)
            else:
                source_only.setdefault(label, []).append(idx)

        # Need at least two distinct source labels to show the report
        all_labels = set(source_only.keys())
        for combo in shared_labels.values():
            for part in combo.split(' + '):
                all_labels.add(part.strip())
        if len(all_labels) < 2:
            return

        # Keep report order aligned with selected source priority
        ordered_labels = []
        if self.current_source_mode == 'combine' and self.combined_source_ids:
            from .database.config import SOURCE_ID_MAP
            for sid in self.combined_source_ids:
                lbl = SOURCE_ID_MAP.get(sid, {}).get('name', f'Source {sid}')
                if lbl in source_only:
                    ordered_labels.append(lbl)
        for lbl in sorted(source_only.keys()):
            if lbl not in ordered_labels:
                ordered_labels.append(lbl)

        # Build source_id lookup for dedup stats
        sid_for_label = {}
        if self.current_source_mode == 'combine' and self.combined_source_ids:
            from .database.config import SOURCE_ID_MAP
            for sid in self.combined_source_ids:
                lbl = SOURCE_ID_MAP.get(sid, {}).get('name', f'Source {sid}')
                sid_for_label[lbl] = sid

        print("\n" + "-" * 70)
        print(f"  INSTANCES BY SOURCE - {motif_type}")
        print("-" * 70)

        # Show within-source deduplication counts if available
        if self.dedup_stats:
            print("\n  Within-source deduplication:")
            for label in ordered_labels:
                sid = sid_for_label.get(label)
                if sid and sid in self.dedup_stats:
                    before, after = self.dedup_stats[sid]
                    removed = before - after
                    if removed > 0:
                        print(f"    {label}: {before} -> {after} (removed {removed} duplicates)")
                    else:
                        print(f"    {label}: {after} (no duplicates)")
            # Also show labels not in source_only (they may only appear in shared)
            for label in all_labels:
                if label not in ordered_labels:
                    sid = sid_for_label.get(label)
                    if sid and sid in self.dedup_stats:
                        before, after = self.dedup_stats[sid]
                        removed = before - after
                        if removed > 0:
                            print(f"    {label}: {before} -> {after} (removed {removed} duplicates)")
                        else:
                            print(f"    {label}: {after} (no duplicates)")

        total = len(motif_details)
        for label in ordered_labels:
            ids = source_only.get(label, [])
            count = len(ids)
            id_str = ', '.join(str(i) for i in ids) if ids else '-'
            print(f"\n  Unique in {label}: {count} instance(s)")
            print(f"    IDs: {id_str}")

        if shared_ids:
            id_str = ', '.join(str(i) for i in shared_ids)
            # Use the combo label from the first shared instance as header
            combo = shared_labels.get(shared_ids[0], 'Shared')
            print(f"\n  Shared ({combo}): {len(shared_ids)} instance(s)")
            print(f"    IDs: {id_str}")

        print(f"\n  Total merged instances: {total}")
        print("-" * 70)

    def _resolve_loaded_motif_type(self, motif_type: str, loaded_motifs: Optional[Dict] = None) -> str:
        """Resolve user motif aliases (HL/IL/Jn) to an available loaded motif key.

        This keeps CLI examples stable across sources where keys may be short
        aliases (e.g., HL) or semantic names (e.g., HAIRPIN LOOP).
        """
        raw = str(motif_type or '').strip()
        if not raw:
            return raw

        if loaded_motifs is None:
            loaded_motifs = self._get_current_source_motifs() or {}

        if not loaded_motifs:
            return raw.upper()

        def _norm(s: str) -> str:
            return ''.join(ch for ch in (s or '').upper() if ch.isalnum())

        target = raw.upper()
        target_norm = _norm(target)

        # 1) Exact canonical match (case/space/punctuation-insensitive)
        for key in loaded_motifs.keys():
            if _norm(key) == target_norm:
                return key

        # 1b) Cross-check alignment.py's canonical alias table (KINK-TURN ->
        # K-TURN, SARCIN -> SARCIN-RICIN, etc.) so rmv_show/rmv_view/
        # rmv_summary accept the same synonyms rmv_super/rmv_align already do.
        try:
            from .alignment import MOTIF_ALIASES
        except Exception:
            MOTIF_ALIASES = {}
        for alias_key, canonical in MOTIF_ALIASES.items():
            if _norm(alias_key) == target_norm:
                canonical_norm = _norm(canonical)
                for key in loaded_motifs.keys():
                    if _norm(key) == canonical_norm:
                        return key

        # 2) Common loop aliases used in help/docs
        loop_aliases = {
            'HL': ['HAIRPIN LOOP', 'HAIRPINLOOP', 'HAIRPIN'],
            'IL': ['INTERNAL LOOP', 'INTERNALLOOP', 'INTERNAL'],
        }
        junction_match = re.fullmatch(r'J(\d+)', target)
        if junction_match:
            n = junction_match.group(1)
            loop_aliases[target] = [f'{n}-WAY JUNCTION', f'{n}WAYJUNCTION', f'{n} WAY JUNCTION']

        for alias, candidates in loop_aliases.items():
            if target != alias:
                continue
            candidate_norms = {_norm(c) for c in candidates}
            for key in loaded_motifs.keys():
                key_upper = key.upper()
                key_norm = _norm(key_upper)
                key_no_parens_norm = _norm(re.sub(r'\s*\([^)]*\)\s*', ' ', key_upper))
                if key_norm in candidate_norms or key_no_parens_norm in candidate_norms:
                    return key

        # 3) Heuristic fallback for semantic names carrying alias in parentheses
        if target in ('HL', 'IL'):
            needle = 'HAIRPIN LOOP' if target == 'HL' else 'INTERNAL LOOP'
            for key in loaded_motifs.keys():
                if needle in key.upper():
                    return key

        return target

    def _resolve_source_filter(self, motif_type: str, source_filter: str):
        """Resolve a source-filter keyword to instance numbers for a motif type.

        Args:
            motif_type: Uppercased motif type key (e.g., 'K-TURN')
            source_filter: Case-insensitive keyword - a source name/alias
                           (e.g., 'nobias', 'rmsx') or 'shared'.

        Returns:
            list[int] | None: 1-based instance numbers matching the filter,
                              or None if the filter didn't match anything
                              (caller should treat the word as part of the
                              motif name instead).
        """
        # Treat any state where two or more sources are combined as combine
        # mode for filter resolution.  current_source_mode can drift to other
        # values after `rmv_load_motif` reload sequences, but if
        # combined_source_ids still has 2+ entries the alias filter is still
        # meaningful.
        if self.current_source_mode != 'combine' and len(self.combined_source_ids or []) < 2:
            return None

        loaded_motifs = self._get_current_source_motifs()
        # Case + whitespace insensitive lookup so 'KINK-TURN'/'kink-turn'/
        # ' KINK-TURN ' all resolve to the same key as the loader stored.
        def _key_norm(s: str) -> str:
            return ''.join((s or '').upper().split())
        target_key = _key_norm(motif_type)
        resolved_motif_key = None
        for k in loaded_motifs.keys():
            if _key_norm(k) == target_key:
                resolved_motif_key = k
                break
        if resolved_motif_key is None:
            self.logger.debug(
                f"_resolve_source_filter: motif '{motif_type}' not found in "
                f"loaded_motifs (keys={list(loaded_motifs.keys())})")
            return None
        motif_type = resolved_motif_key

        motif_details = loaded_motifs[motif_type].get('motif_details', [])
        if not motif_details:
            self.logger.debug(
                f"_resolve_source_filter: motif '{motif_type}' has no details")
            return None

        # Categorise instances exactly like _print_source_attribution_report
        source_only = {}   # source_label -> [1-based idx]
        shared_ids = []    # 1-based idx list

        for idx, detail in enumerate(motif_details, 1):
            meta = detail.get('metadata', {})
            label = meta.get('_source_label', '')
            also = meta.get('_also_found_in', [])
            if not label:
                continue
            if also:
                shared_ids.append(idx)
            else:
                source_only.setdefault(label, []).append(idx)

        sf = source_filter.upper()

        # --- "shared" keyword ---
        if sf == 'SHARED':
            if shared_ids:
                return shared_ids
            return []

        # --- Match against source labels ---
        # Build a lookup: lowercase fragments -> full label
        # Accept: full name, full name w/o parentheses, tool shorthand,
        # subtype shorthand, or any single word in the name.
        import re as _re
        from .database.config import SOURCE_ID_MAP
        alias_to_label = {}  # alias (upper, normalised) -> full label

        def _norm(s: str) -> str:
            # Collapse whitespace and uppercase for case/spacing-insensitive
            # matching of multi-token source names.
            return " ".join(s.split()).upper()

        for sid in (self.combined_source_ids or []):
            info = SOURCE_ID_MAP.get(sid, {})
            full_name = info.get('name', '')
            if not full_name:
                continue
            # Full name (with parentheses), normalised
            alias_to_label[_norm(full_name)] = full_name
            # Full name with parentheses stripped - e.g. "RNAMotifScanX"
            # for "RNAMotifScanX (RMSX)"
            no_parens = _re.sub(r'\s*\([^)]*\)\s*', ' ', full_name)
            alias_to_label[_norm(no_parens)] = full_name
            # Tool shorthand (e.g., 'rmsx', 'nobias', 'rms', 'fr3d')
            tool = info.get('tool', '')
            if tool:
                alias_to_label[_norm(tool)] = full_name
            # Subtype shorthand for local/web sources (e.g., 'bgsu',
            # 'atlas', 'rfam', 'rfam_api') - these sources use 'subtype'
            # rather than 'tool' in SOURCE_ID_MAP.
            subtype = info.get('subtype', '')
            if subtype:
                alias_to_label[_norm(subtype)] = full_name
                alias_to_label[_norm(subtype.replace('_', ''))] = full_name
                alias_to_label[_norm(subtype.replace('_', ' '))] = full_name
            # Each word in the name (e.g., 'RNAMOTIFSCANX', 'RMSX')
            for word in full_name.replace('(', ' ').replace(')', ' ').split():
                alias_to_label[_norm(word)] = full_name

        matched_label = alias_to_label.get(_norm(source_filter))
        if matched_label is None:
            self.logger.debug(
                f"_resolve_source_filter: alias '{source_filter}' not in "
                f"{sorted(alias_to_label.keys())}")
            return None  # not a recognized source filter

        ids = source_only.get(matched_label, [])
        return sorted(ids)

    def show_motif_instance_summary(self, motif_type: str, instance_no: int):
        """Print residue-based details of a specific matched instance (for rmv_summary MOTIF NO).
        
        Args:
            motif_type (str): Motif type (e.g., 'HL', 'K-TURN')
            instance_no (int): Instance number (1-indexed, matches the table row order)
        """
        from .database.motif_hierarchy_cache import get_hierarchy_cache

        motif_arg = motif_type.upper().strip()
        ops = self._effective_db_ops([])
        pdb_id, rows = self._matching_query_rows(motif_arg, ops)

        if not pdb_id:
            print("\nNo structure loaded. Use 'rmv_fetch <PDB_ID>' first.\n")
            return
        if not rows:
            print(f"\nNo motif instances matching '{motif_arg}' found for {pdb_id}.\n")
            return
        if instance_no < 1 or instance_no > len(rows):
            print(f"\nInstance {instance_no} not found. Valid range: 1-{len(rows)}\n")
            return

        _custom_id, r_key, labels_by_source = rows[instance_no - 1]
        print()
        self._print_query_context(ops, pdb_id)
        print(f"\n{'='*70}")
        print(f"  {motif_arg} INSTANCE #{instance_no}")
        print('='*70)
        print(f"  PDB      : {pdb_id}")
        print(f"  Residues : {_format_residue_key_ranges(r_key)}")
        for sid, labels in labels_by_source.items():
            from .database.config import SOURCE_ID_MAP
            src_name = SOURCE_ID_MAP.get(sid, {}).get('name', f'Source {sid}')
            # Show every stored hierarchy level (context -> subclass).
            print(f"  {src_name}({sid}): {' > '.join(labels)}")
        print(f"  Hierarchy cache: {get_hierarchy_cache().db_path}")
        print('='*70)
        
    def set_source_mode(self, mode: str):
        """
        Set the motif data source mode.
        
        Args:
            mode (str): Source mode: auto, local, web, bgsu, rfam, all, user
        """
        try:
            mode_lower = mode.lower()
            
            # Handle user annotations specially
            if mode_lower == 'user':
                self._set_user_annotations_source()
                return
            
            from .database import get_config, SourceMode
            
            mode_map = {
                'auto': SourceMode.AUTO,
                'local': SourceMode.LOCAL,
                'web': SourceMode.AUTO,        # web mode uses AUTO (smart selection)
                'bgsu': SourceMode.BGSU,
                'rfam': SourceMode.RFAM,
                'all': SourceMode.ALL
            }
            
            if mode_lower not in mode_map:
                valid_modes = ['auto', 'local', 'web', 'web bgsu', 'web rfam', 'local atlas', 'local rfam', 'all', 'user fr3d', 'user rnamotifscan', 'user rnamotifscanx']
                self.logger.error(f"Invalid source mode '{mode}'.")
                self.logger.info("Valid source modes:")
                for m in valid_modes:
                    self.logger.info(f"  rmv_db {m}")
                return
            
            config = get_config()
            config.source_mode = mode_map[mode_lower]
            
            # BUG FIX: Clear specific_source when using generic modes (auto, all)
            # Note: specific source handlers (_handle_local_source_by_id, etc.) set this explicitly
            if mode_lower in ['auto', 'all']:
                config.specific_source = None
            
            mode_display = mode_lower if mode_lower != 'web' else 'web (auto-select online APIs)'
            # Internal adapter mode; the public source name is reported by the
            # caller (e.g. 'Source: RNA3DMotifAtlas'), so keep this off-console.
            self.logger.debug(f"Motif source mode set to: {mode_display}")

        except Exception as e:
            self.logger.error(f"Failed to set source mode: {e}")
    
    def _set_user_annotations_source(self):
        """Set source to user annotations with tool selection."""
        print("\n" + "="*60)
        print("USER ANNOTATIONS")
        print("="*60)
        print("\nAvailable tools:")
        print("  1. fr3d           - FR3D output format (BGSU base pairs)")
        print("  2. rnamotifscan   - RNAMotifScan output format (RMS)")
        print("  3. rnamotifscanx  - RNAMotifScanX output format (RMSX)")
        print("\nAfter selecting a tool with rmv_db user <TOOL>,")
        print("use rmv_fetch to load structures:")
        print("\nUsage:")
        print("  rmv_fetch <PDB_ID>")
        print("\nExample:")
        print("  rmv_fetch 1S72")
        print("="*60 + "\n")
        
        # Store that user annotations are selected
        self.current_source_mode = 'user'
        self.logger.success("User Annotations mode selected")
        self.logger.info("Use: rmv_fetch <PDB_ID>")
        self.logger.info("Tools: fr3d, rnamotifscan, rnamotifscanx")
    
    def _print_source_mode_info(self):
        """Print information about current source mode."""
        try:
            from .database import get_config, SourceMode
            
            config = get_config()
            mode = config.source_mode
            
            mode_descriptions = {
                SourceMode.AUTO: "Use the configured named-source adapter",
                SourceMode.LOCAL: "Legacy compatibility mode; offline databases are not shipped",
                SourceMode.BGSU: "Use the RNA 3D Motif Atlas online adapter",
                SourceMode.RFAM: "Use the Rfam online adapter",
                SourceMode.ALL: "Combine the explicitly selected named sources"
            }
            
            print(f"\nCurrent mode: {mode.value}")
            print(f"Description: {mode_descriptions.get(mode, 'Unknown')}")
            
        except ImportError:
            print("Source selector not available")
    
    def _handle_source_by_id(self, source_id: int, extra_args: str = None):
        """Compatibility handler for internal provider dispatch.
        
        Also detects multi-source mode when extra_args contains additional
        numeric source IDs (e.g., 'rmv_db 1 3' or 'rmv_db 2 5 3').
        
        Args:
            source_id (int): Source ID (1-8)
            extra_args (str): Optional arguments:
                - Additional source IDs for multi-source combine (e.g., "3" or "5 3")
                - External paths used by compatibility adapters.
        """
        if source_id not in SOURCE_ID_MAP:
            self.logger.error(f"Invalid source ID: {source_id}")
            self.logger.error("Valid source IDs:")
            for sid, info in SOURCE_ID_MAP.items():
                self.logger.error(f"  {sid} = {info['name']}")
            return
        
        # --- Multi-source detection ---
        # If extra_args contains ONLY numeric source IDs, enter combine mode.
        # e.g., rmv_db 1 3 -> source_id=1, extra_args='3'
        # e.g., rmv_db 2 5 3 -> source_id=2, extra_args='5 3'
        if extra_args:
            extra_parts = str(extra_args).strip().split()
            all_numeric = all(p.isdigit() for p in extra_parts)
            if all_numeric and extra_parts:
                # All extra args are numbers -> multi-source combine mode
                all_ids = [source_id] + [int(p) for p in extra_parts]
                # Validate all IDs
                invalid = [sid for sid in all_ids if sid not in SOURCE_ID_MAP]
                if invalid:
                    self.logger.error(f"Invalid source ID(s): {invalid}")
                    self.logger.error("Valid source IDs: " + 
                                     ", ".join(f"{k}={v['name']}" for k, v in SOURCE_ID_MAP.items()))
                    return
                if len(all_ids) != len(set(all_ids)):
                    self.logger.error("Duplicate source IDs not allowed")
                    return
                # Enter combine mode
                self._handle_multi_source(all_ids)
                return
        
        # Store the numeric source ID for object tagging
        self.current_source_id = source_id
        
        source_info = SOURCE_ID_MAP[source_id]
        source_type = source_info['type']
        
        # Handle different source types (single source mode)
        if source_type == 'local':
            self._handle_local_source_by_id(source_id, source_info, extra_args)
        elif source_type == 'web':
            self._handle_web_source_by_id(source_id, source_info, extra_args)
        elif source_type == 'user':
            self._handle_user_source_by_id(source_id, source_info, extra_args)
        elif source_type == 'analysis':
            self.logger.error(f"Source {source_id} is not currently available.")
            return
        else:
            self.logger.error(f"Unknown source type: {source_type}")
    
    def _handle_multi_source(self, source_ids: list):
        """Handle multi-source combine mode.
        
        Called when user provides multiple source IDs:
            rmv_db 1 3     -> source_ids=[1, 3]
            rmv_db 2 5 3   -> source_ids=[2, 5, 3]
        
        Priority order = left to right (first = highest priority).
        
        Args:
            source_ids: List of source IDs in priority order
        """
        self.combined_source_ids = source_ids
        self.current_source_mode = 'combine'
        self.current_local_source = None
        self.current_web_source = None
        self.current_user_tool = None
        # Use combined IDs joined for suffix (e.g., S1_3 for sources 1+3)
        self.current_source_id = '_'.join(str(s) for s in source_ids)
        
        # Clear specific_source in config
        from .database import get_config
        config = get_config()
        config.specific_source = None
        
        for i, sid in enumerate(source_ids, 1):
            info = SOURCE_ID_MAP[sid]
            public_name = (
                self.current_source_names[i - 1]
                if i - 1 < len(self.current_source_names)
                else info['name']
            )
            # Show p-value status for RMS/RMSX/NoBIAS
            pval_note = ""
            tool = info.get('tool', '')
            if tool in ['rms', 'rnamotifscan']:
                if self.user_rms_custom_pvalues:
                    pv = ", ".join(f"{m}={p}" for m, p in self.user_rms_custom_pvalues.items())
                    pval_note = f" | P-values: {pv}"
                else:
                    fs = "ON" if self.user_rms_filtering_enabled else "OFF"
                    pval_note = f" | Filtering: {fs}"
            elif tool in ['rmsx', 'rnamotifscanx']:
                if self.user_rmsx_custom_pvalues:
                    pv = ", ".join(f"{m}={p}" for m, p in self.user_rmsx_custom_pvalues.items())
                    pval_note = f" | P-values: {pv}"
                else:
                    fs = "ON" if self.user_rmsx_filtering_enabled else "OFF"
                    pval_note = f" | Filtering: {fs}"
            elif tool in ['nobias']:
                if self.user_nobias_custom_pvalues:
                    pv = ", ".join(f"{m}={p}" for m, p in self.user_nobias_custom_pvalues.items())
                    pval_note = f" | P-values: {pv}"
                else:
                    fs = "ON" if self.user_nobias_filtering_enabled else "OFF"
                    pval_note = f" | Filtering: {fs}"
            # Show custom path indicator if this source has one
            path_note = ""
            if self.user_data_paths.get(sid):
                path_note = " - custom path"
            self.logger.info(f"  {public_name}{pval_note}{path_note}")

        # The RMSX P-value tip is only relevant when RNAMotifScanX is selected.
        if 7 in source_ids:
            self.logger.info("Tip: RMSX P-value thresholds are configured in config/rmsx_config.json.")
    
    def _handle_local_source_by_id(self, source_id: int, source_info: Dict, extra_args: str = None):
        """Handle local source selection by ID."""
        subtype = source_info.get('subtype')
        
        self.current_source_mode = 'local'
        self.current_local_source = subtype
        self.current_web_source = None
        self.current_user_tool = None
        self.combined_source_ids = []
        
        # BUG FIX: Set specific_source in config to ensure ONLY this source is used
        from .database import get_config
        config = get_config()
        config.specific_source = subtype  # e.g., 'atlas' for source 1
        
        self.set_source_mode('local')
        self.logger.debug(f"Set config.specific_source = {subtype}")
        self.logger.success(f"Source: {source_info['name']}")
        
        # VERIFICATION: Print source configuration
        self.logger.debug(f"SOURCE CONFIG VERIFICATION:")
        self.logger.debug(f"  - self.current_source_mode = {self.current_source_mode}")
        self.logger.debug(f"  - self.current_local_source = {self.current_local_source}")
        self.logger.debug(f"  - config.specific_source = {config.specific_source}")
        self.logger.debug(f"  - Expected to load from: {subtype} ONLY")
    
    def _handle_web_source_by_id(self, source_id: int, source_info: Dict, extra_args: str = None):
        """Handle online source selection by ID."""
        subtype = source_info.get('subtype')
        
        self.current_source_mode = 'web'
        self.current_web_source = subtype
        self.current_local_source = None
        self.current_user_tool = None
        self.combined_source_ids = []
        
        # BUG FIX: Set specific_source in config to ensure ONLY this source is used
        from .database import get_config
        config = get_config()
        # Map subtype to provider ID
        subtype_to_provider = {'bgsu': 'bgsu_api', 'bgsu_api': 'bgsu_api', 'rfam_api': 'rfam_api'}
        provider_id = subtype_to_provider.get(subtype, subtype)
        config.specific_source = provider_id
        
        # Map subtype to SourceMode
        mode_map = {'bgsu': 'bgsu', 'bgsu_api': 'bgsu', 'rfam_api': 'rfam'}
        mode = mode_map.get(subtype, 'auto')
        self.set_source_mode(mode)
        self.logger.debug(f"Set config.specific_source = {provider_id} (from subtype={subtype})")
        
        self.logger.success(f"Source: {source_info['name']}")
        
        # VERIFICATION: Print source configuration
        self.logger.debug(f"SOURCE CONFIG VERIFICATION:")
        self.logger.debug(f"  - self.current_source_mode = {self.current_source_mode}")
        self.logger.debug(f"  - self.current_web_source = {self.current_web_source}")
        self.logger.debug(f"  - config.specific_source = {config.specific_source}")
        self.logger.debug(f"  - Expected to load from: {provider_id} ONLY")
    
    def _handle_user_source_by_id(self, source_id: int, source_info: Dict, extra_args: str = None):
        """Handle user annotation source selection by ID.
        
        Supports:
        - rmv_db 6              (RMS with default filtering ON)
        - rmv_db 6 off          (RMS with filtering OFF)
        - rmv_db 6 on           (RMS with filtering ON - explicit)
        - rmv_db 6 C-LOOP 0.05 KINK-TURN 0.02  (RMS with custom P-values)
        - rmv_db 6 /path/to/data   (RMS with custom data directory)
        - rmv_db 7 /path/to/data   (RMSX with custom data directory)
        """
        tool = source_info.get('tool')
        
        if not tool:
            self.logger.error(f"Source {source_id} is not a user annotation source")
            return
        
        self.current_source_mode = 'user'
        self.current_user_tool = tool
        self.current_local_source = None
        self.current_web_source = None
        self.combined_source_ids = []
        # Per-source custom paths: no reset needed - each source has its own entry
        
        # BUG FIX: Clear specific_source for user annotation sources (they use tool-based loading)
        from .database import get_config
        config = get_config()
        config.specific_source = None
        
        # Parse extra arguments
        custom_pvalues = {}
        filtering_enabled = True  # Default: ON
        
        if extra_args:
            extra_str = str(extra_args).strip()
            
            # Strip surrounding quotes (PyMOL may preserve them from user input)
            if len(extra_str) >= 2 and (
                (extra_str[0] == "'" and extra_str[-1] == "'") or
                (extra_str[0] == '"' and extra_str[-1] == '"')
            ):
                extra_str = extra_str[1:-1].strip()
            
            # Check if the argument looks like a file path
            import os
            if extra_str.startswith('/') or extra_str.startswith('~') or \
               extra_str.startswith('./') or extra_str.startswith('..'):
                # Expand user home directory
                expanded_path = os.path.expanduser(extra_str)

                if source_id in (5, 7):
                    self.logger.warning(
                        "External source paths and configs are fixed by this project; "
                        "the supplied path was ignored."
                    )
                elif os.path.isdir(expanded_path):
                    self.user_data_paths[source_id] = expanded_path
                    if source_id == 7:
                        self.rmsx_output_path = os.path.abspath(os.path.expanduser(expanded_path))
                        os.makedirs(self.rmsx_output_path, exist_ok=True)
                        self.rmsx_pipeline_config = self._build_internal_rmsx_config()
                    self.logger.success(f"Custom data path set: {expanded_path}")
                elif os.path.isfile(expanded_path):
                    self.user_data_paths[source_id] = expanded_path
                    self.logger.success(f"Custom data file set: {expanded_path}")
                else:
                    self.logger.warning(f"Path not found: {expanded_path}")
                    self.logger.warning("Will use default data directory instead")
            else:
                self.logger.warning(
                    "RMSX filtering and P-value thresholds are configured only in "
                    "config/rmsx_config.json; command-line overrides are ignored."
                )
        elif source_id == 5:
            if not self.register_fr3d_source(str(self.default_fr3d_config)):
                self.logger.error("FR3D (Source 5) is not ready yet - see the precise reason above.")
                self.logger.info("  Run this once to install everything required and register FR3D:")
                self.logger.info("    rmv_setup FR3D")
                return
        
        # Store filtering state and custom P-values
        if tool in ['rms', 'rnamotifscan']:
            self.user_rms_filtering_enabled = filtering_enabled
            self.user_rms_custom_pvalues = custom_pvalues
        elif tool in ['rmsx', 'rnamotifscanx']:
            self.user_rmsx_filtering_enabled = filtering_enabled
            if custom_pvalues:
                self.user_rmsx_custom_pvalues = custom_pvalues
            if source_id == 7:
                self.rmsx_pipeline_config = self._build_internal_rmsx_config()
        elif tool in ['nobias']:
            self.user_nobias_filtering_enabled = filtering_enabled
            self.user_nobias_custom_pvalues = custom_pvalues
        
        # Build status message
        public_name = self.current_source_names[0] if len(self.current_source_names) == 1 else source_info['name']
        status_msg = f"Source: {public_name}"
        if tool in ['rms', 'rnamotifscan', 'rmsx', 'rnamotifscanx', 'nobias']:
            if custom_pvalues:
                pvalue_str = ", ".join([f"{m}={p}" for m, p in custom_pvalues.items()])
                status_msg += f" | Custom P-values: {pvalue_str}"
            else:
                filter_status = "Filtering: ON" if filtering_enabled else "Filtering: OFF"
                status_msg += f" | {filter_status}"
        
        self.logger.success(status_msg)
        _udp = self.user_data_paths.get(source_id)
        if _udp:
            self.logger.info(f"  Data path: {_udp}")
        
    def _handle_source_info_command(self, source_id_str: str = None):
        """Display information about the active source, or detailed info about a specific source.
        
        Usage:
            rmv_source info      - Show currently active source
            rmv_source info <N>  - Show detailed info about source N
        """
        if not source_id_str:
            # Show currently active source only
            self._print_active_source_info()
            return
        
        try:
            source_id = int(source_id_str.strip())
            self._print_single_source_info(source_id)
        except ValueError:
            self.logger.error(f"Invalid source ID: {source_id_str}")
            self.command_error("Usage: rmv_source info [<ID>] or rmv_db")
    
    def _print_active_source_info(self):
        """Print info about the currently active source only."""
        if self.current_source_id is None or self.current_source_mode is None:
            self.logger.info("No source is currently active.")
            self.logger.info("")
            self.logger.info("Select a source first:")
            self.logger.info("  rmv_db      List all available sources")
            self.logger.info("  rmv_db <N>       Select source (1-8)")
            return
        
        # Handle combine mode (multiple sources)
        if self.current_source_mode == 'combine' and self.combined_source_ids:
            print("\n" + "=" * 60)
            print("  ACTIVE SOURCE: COMBINED MODE")
            print("=" * 60)
            print(f"\n  Sources combined: {', '.join(str(s) for s in self.combined_source_ids)}")
            for sid in self.combined_source_ids:
                if sid in SOURCE_ID_MAP:
                    info = SOURCE_ID_MAP[sid]
                    print(f"    [{sid}] {info['name']:30} | {info['coverage']}")
            print(f"\n  Mode:     combine")
            if self.loaded_pdb_id:
                print(f"  PDB:      {self.loaded_pdb_id}")
            print("\n" + "=" * 60 + "\n")
            return
        
        # Single source mode - determine source ID
        try:
            source_id = int(self.current_source_id)
        except (ValueError, TypeError):
            source_id = None
        
        if source_id and source_id in SOURCE_ID_MAP:
            info = SOURCE_ID_MAP[source_id]
            
            print("\n" + "=" * 60)
            print(f"  ACTIVE SOURCE: [{source_id}] {info['name'].upper()}")
            print("=" * 60)
            print(f"\n  Description:  {info['description']}")
            print(f"  Category:     {info.get('category', 'N/A')}")
            print(f"  Coverage:     {info['coverage']}")
            print(f"  Type:         {info['type'].upper()}")
            
            if self.loaded_pdb_id:
                print(f"  PDB loaded:   {self.loaded_pdb_id}")
            
            # Show filtering status for RMS/RMSX/NoBIAS
            if source_id == 6:
                status = "ON" if self.user_rms_filtering_enabled else "OFF"
                print(f"  Filtering:    {status}")
                if self.user_rms_custom_pvalues:
                    print(f"  Custom P-values: {self.user_rms_custom_pvalues}")
            elif source_id == 7:
                status = "ON" if self.user_rmsx_filtering_enabled else "OFF"
                print(f"  Filtering:    {status}")
                if self.user_rmsx_custom_pvalues:
                    print(f"  Custom P-values: {self.user_rmsx_custom_pvalues}")
            elif source_id == 8:
                status = "ON" if self.user_nobias_filtering_enabled else "OFF"
                print(f"  Filtering:    {status}")
                if self.user_nobias_custom_pvalues:
                    print(f"  Custom P-values: {self.user_nobias_custom_pvalues}")
            
            # Show custom data path if set
            _udp = self.user_data_paths.get(source_id)
            if _udp and source_id in (5, 6, 7, 8):
                print(f"  Data path:    {_udp}")
            
            print("\n" + "=" * 60 + "\n")
        else:
            self.logger.info(f"Active source mode: {self.current_source_mode}")
            self.logger.info(f"Source ID: {self.current_source_id}")
    
    def _print_single_source_info(self, source_id: int):
        """Print detailed information about a single source."""
        if source_id not in SOURCE_ID_MAP:
            self.logger.error(f"Source ID {source_id} not found")
            return
        
        info = SOURCE_ID_MAP[source_id]
        
        print("\n" + "="*70)
        print(f"  SOURCE {source_id}: {info['name'].upper()}")
        print("="*70)
        
        print(f"\nDescription:  {info['description']}")
        print(f"Category:     {info.get('category', 'N/A')}")
        print(f"Coverage:     {info['coverage']}")
        print(f"Type:         {info['type'].upper()}")
        
        # Source-specific information
        if info['type'] == 'user':
            if source_id == 5:
                # FR3D - external official BGSU fr3d-python (registered via config)
                print("\n--- Pipeline ---")
                print("  Source 5 wraps an external, user-installed official BGSU")
                print("  fr3d-python. RSMViewer never vendors or modifies that repo;")
                print("  it only runs it in a subprocess and ingests the resulting CSV.")
                print("\n--- Execution order ---")
                print("  1. Resolve the registered FR3D repo + Python environment")
                print("  2. Run the official FR3D search on the loaded structure")
                print("  3. Ingest the resulting CSV and load motifs into RSMViewer")
                print("\n--- Current configuration ---")
                registered = bool(getattr(self, 'fr3d_registered', False))
                print(f"  Registered   : {'yes' if registered else 'no'}")
                print(f"  Config path  : {getattr(self, 'fr3d_config_path', '') or '(not registered)'}")
                print(f"  FR3D repo    : {getattr(self, 'fr3d_python_path', '') or '(not registered)'}")
                print(f"  Python exe   : {getattr(self, 'fr3d_python_exe', '') or '(auto-detect)'}")
                print(f"  Data mode    : {getattr(self, 'fr3d_data_mode', 'cif_local')}")
                print(f"  Output dir   : {getattr(self, 'fr3d_output_dir', '')}")
                print("\n--- Setup ---")
                print("  1. Place the official checkout under external/fr3d/fr3d-python-latest")
                print("  2. Review the fixed config/fr3d_config.json")
                print("  3. Run: rmv_db FR3D")
            elif source_id == 7:
                pipe_cfg = getattr(self, 'rmsx_pipeline_config', {})
                print("\n--- Pipeline ---")
                print("  Runs RNAMotifScanX per motif family, parses result_*.log files,")
                print("  and loads source 7 motifs directly into RSMViewer.")
                print("\n--- Current configuration ---")
                print(f"  Mode         : {pipe_cfg.get('data_mode', 'preannotated')}")
                print(f"  Runtime dir  : {self.rmsx_runtime_dir}")
                print(f"  Executable   : {pipe_cfg.get('rmsx_executable', '(not found)')}")
                print(f"  Output dir   : {self.rmsx_output_path}")
                print(f"  Families     : {', '.join(pipe_cfg.get('motif_families', [])) or 'all defaults'}")
                print("\n--- Suggested workflow ---")
                print("  1. rmv_fetch <PDB_ID>")
                print("  2. rmv_db RNAMotifScanX")
                print("  3. rmv_select <MOTIF>, <PDB_ID>, RNAMotifScanX, as <GROUP>")
                print("  4. rmv_list <GROUP>   then   rmv_view <GROUP>")
            else:
                print(f"\nAvailable motif types will be shown after loading a structure")
                print(f"with rmv_fetch <PDB_ID>")
        
        # Display sample commands
        print(f"\n--- Sample commands ---")
        if source_id == 5:
            print(f"  rmv_db FR3D                            Load FR3D annotations")
            print(f"  rmv_fetch 1S72                          Load PDB structure")
            print(f"  rmv_db FR3D                             Run FR3D search, then load")
            print(f"  rmv_list                                List loaded motif rows")
            print(f"  rmv_select Sarcin-Ricin, 1S72, FR3D, as group_SR  Select a family")
            print(f"  rmv_view group_SR                       Highlight the group")
            print(f"  rmv_fr3d status                         Show FR3D registration status")
            print(f"\n--- Combine FR3D with other sources ---")
            print(f"  rmv_db FR3D,RNAMotifScanX               FR3D + RMSX combined")
            print(f"  rmv_db RNA3DMotifAtlas,FR3D             Atlas + FR3D combined")
        else:
            print(f"  rmv_fetch 1S72                  Load structure")
            print(f"  rmv_db <SOURCE>                 Load annotations")
            print(f"  rmv_list                        List loaded motif rows")
            print(f"  rmv_view <MOTIF_ID|group>       Highlight a motif or group")
        
        # RMS/RMSX/NoBIAS specific features
        if info.get('supports_filtering'):
            print(f"\nWith filtering control:")
            print(f"  rmv_db {source_id} off              Disable filtering (show all motifs)")
            print(f"  rmv_db {source_id} on               Enable filtering (default)")
            print(f"\nWith custom P-values:")
            print(f"  rmv_db {source_id} C-LOOP 0.05 KINK-TURN 0.02")
            print(f"    -> Apply custom thresholds for specific motif types")
            print(f"    -> Other motif types use default thresholds")
        
        print("\n" + "="*70 + "\n")
    
    def _print_all_source_info(self):
        """Print summary info for all sources."""
        print("\n" + "="*70)
        print("  AVAILABLE DATA SOURCES (Quick Reference)")
        print("="*70)
        
        current_category = None
        for source_id in sorted(SOURCE_ID_MAP.keys()):
            info = SOURCE_ID_MAP[source_id]
            
            if info.get('category') != current_category:
                current_category = info.get('category')
                print(f"\n{current_category}:")
                print("-" * 70)
            
            print(f"  [{source_id}] {info['name']:30} | {info['coverage']:20} | {info['description']}")
        
        print("\n" + "="*70)
        print("Usage:")
        print("  rmv_db <ID>                    Select source by ID")
        print("  rmv_source info <ID>           Show detailed info")
        print("  rmv_db                    List all sources (this display)")
        print("="*70 + "\n")
    
    def _handle_user_source(self, tool_name):
        """Handle user annotations source selection."""
        if not tool_name:
            self.command_error("Usage: rmv_db user <tool_name> [on|off]")
            self.logger.error("Available tools:")
            self.logger.error("  rmv_db user fr3d")
            self.logger.error("  rmv_db user rms [on|off]          (default: on)")
            self.logger.error("  rmv_db user rmsx [on|off]         (default: on)")
            return
        
        # Parse tool name and optional on/off parameter
        parts = str(tool_name).strip().split()
        tool = parts[0].lower()
        filtering_enabled = True  # Default: filters ON
        
        # Check for optional on/off parameter (only for rms and rmsx)
        if len(parts) > 1:
            filter_arg = parts[1].lower()
            if filter_arg in ['on', 'off']:
                filtering_enabled = (filter_arg == 'on')
            else:
                self.logger.warning(f"Unknown parameter '{filter_arg}'. Expected 'on' or 'off'. Using default: on")
        
        valid_tools = ['fr3d', 'rnamotifscan', 'rms', 'rnamotifscanx', 'rmsx']
        if tool not in valid_tools:
            self.logger.error(f"Invalid tool '{tool}'. Valid options: {', '.join(valid_tools)}")
            return
        
        # Store filtering state for RMS and RMSX
        if tool in ['rms', 'rnamotifscan']:
            self.user_rms_filtering_enabled = filtering_enabled
        elif tool in ['rmsx', 'rnamotifscanx']:
            self.user_rmsx_filtering_enabled = filtering_enabled
        
        self.current_source_mode = 'user'
        self.current_user_tool = tool
        self.current_local_source = None
        self.current_web_source = None
        
        tool_descriptions = {
            'fr3d': 'FR3D (BGSU base pair annotations)',
            'rnamotifscan': 'RNAMotifScan (RMS - structural motif search)',
            'rms': 'RNAMotifScan (RMS - structural motif search)',
            'rnamotifscanx': 'RNAMotifScanX (RMSX - extended motif search)',
            'rmsx': 'RNAMotifScanX (RMSX - extended motif search)'
        }
        
        # Build status message
        status_msg = f"Source set to user annotations: {tool_descriptions.get(tool, tool)}"
        if tool in ['rms', 'rnamotifscan', 'rmsx', 'rnamotifscanx']:
            filter_status = "Filtering: ON (default cutoffs applied)" if filtering_enabled else "Filtering: OFF (all motifs shown)"
            status_msg += f" | {filter_status}"
        
        self.logger.success(status_msg)
        
    def _handle_local_source(self, source_name):
        """Handle local source selection."""
        if not source_name:
            # Just 'rmv_source local' - use local (both atlas and rfam)
            self.current_source_mode = 'local'
            self.current_local_source = None
            self.current_web_source = None
            self.current_user_tool = None
            self.set_source_mode('local')
            self.logger.info("Using local sources (RNA 3D Motif Atlas + Rfam database)")
            return
        
        # For specific local sources
        if source_name == 'atlas':
            self.current_source_mode = 'local'
            self.current_local_source = 'atlas'
            self.current_web_source = None
            self.current_user_tool = None
            self.set_source_mode('local')
            self.logger.success("Source set to local RNA 3D Motif Atlas")
        elif source_name == 'rfam':
            self.current_source_mode = 'local'
            self.current_local_source = 'rfam'
            self.current_web_source = None
            self.current_user_tool = None
            self.set_source_mode('local')
            self.logger.success("Source set to local Rfam database")
        else:
            self.logger.error(f"Invalid local source '{source_name}'")
            self.logger.error("Valid local sources: atlas, rfam")
    
    def _handle_web_source(self, source_name):
        """Handle web/online source selection."""
        if not source_name:
            # Just 'rmv_source web' - use smart web source selection
            self.current_source_mode = 'web'
            self.current_web_source = None
            self.current_local_source = None
            self.current_user_tool = None
            self.set_source_mode('web')
            self.logger.info("Using online sources (auto-select between BGSU and Rfam APIs)")
            return
        
        # For specific online sources
        if source_name == 'bgsu':
            self.current_source_mode = 'web'
            self.current_web_source = 'bgsu_api'
            self.current_local_source = None
            self.current_user_tool = None
            self.set_source_mode('bgsu')
            self.logger.success("Source set to BGSU RNA 3D Hub API (~3000+ PDBs)")
        elif source_name == 'rfam':
            self.current_source_mode = 'web'
            self.current_web_source = 'rfam_api'
            self.current_local_source = None
            self.current_user_tool = None
            self.set_source_mode('rfam')
            self.logger.success("Source set to Rfam API (named motifs)")
        else:
            self.logger.error(f"Invalid online source '{source_name}'")
            self.logger.error("Valid online sources: bgsu, rfam")
    
    def _handle_combine_sources(self, source_ids_str: str):
        """Handle combining multiple sources.
        
        Args:
            source_ids_str: Space-separated source IDs (e.g., "1 3" or "2 5")
        """
        if not source_ids_str:
            self.command_error("Usage: rmv_db combine <ID1> <ID2> [<ID3> ...]")
            self.logger.error("Example: rmv_db combine 1 3    (combine Atlas + BGSU)")
            self.logger.error("Valid source IDs:")
            self.logger.error("  1 = RNA 3D Motif Atlas (Local)")
            self.logger.error("  2 = Rfam (Local)")
            self.logger.error("  3 = BGSU RNA 3D Hub (Online)")
            self.logger.error("  4 = Rfam API (Online)")
            self.logger.error("  5 = FR3D Annotations (User)")
            self.logger.error("  6 = RNAMotifScan (User)")
            return
        
        # Parse source IDs
        try:
            source_ids = [int(sid.strip()) for sid in source_ids_str.split()]
        except ValueError:
            self.logger.error(f"Invalid source IDs: '{source_ids_str}'")
            self.logger.error("IDs must be integers (1-6)")
            return
        
        # Validate source IDs
        try:
            from .database.source_registry import get_source_registry
            registry = get_source_registry()
            is_valid, msg = registry.validate_source_ids(source_ids)
            
            if not is_valid:
                self.logger.error(msg)
                return
            
            # Store combined source IDs
            self.combined_source_ids = source_ids
            self.current_source_mode = 'combine'
            self.current_local_source = None
            self.current_web_source = None
            self.current_user_tool = None
            
            # BUG FIX: Clear specific_source when combining multiple sources
            from .database import get_config
            config = get_config()
            config.specific_source = None
            
            # Display what we're combining
            source_names = registry.get_source_descriptions(source_ids)
            self.logger.success(f"Combining {len(source_ids)} sources:")
            for i, name in enumerate(source_names, 1):
                self.logger.info(f"  {i}. {name}")
            
            self.logger.info("Use 'rmv_fetch <PDB_ID>' to load and combine data from these sources")
            
        except Exception as e:
            self.logger.error(f"Failed to validate sources: {e}")
    
    def refresh_motifs_action(self, pdb_id: str = None):
        """
        Force refresh cache and collect motif data again.
        
        Uses all active loaded PDBs when no PDB is specified, or only the
        requested PDB when one is specified. Clears cached data and re-fetches
        motif data from the same source(s) that were last used.
        
        Args:
            pdb_id (str): PDB ID to refresh, or all active PDBs when omitted
        """
        try:
            if pdb_id:
                pdb_ids = [str(pdb_id).upper()]
            else:
                active_objects = set(cmd.get_object_list())
                pdb_ids = sorted(
                    structure_id
                    for structure_id, object_name in self.loaded_structures.items()
                    if object_name in active_objects
                )

            if not pdb_ids:
                self.logger.error("No structure loaded. Use rmv_fetch <PDB_ID> first.")
                return
            
            # Determine which source(s) to refresh from
            if not self.current_source_mode:
                self.logger.error("No source selected. Use rmv_db <N> first.")
                return
            
            # Describe what we're refreshing
            if self.current_source_mode == 'combine' and self.combined_source_ids:
                source_desc = f"combined sources {self.combined_source_ids}"
            elif self.current_source_mode == 'user':
                source_desc = f"user annotations ({self.current_user_tool or 'unknown'})"
            elif self.current_source_mode == 'local':
                source_desc = f"local ({self.current_local_source or 'auto'})"
            elif self.current_source_mode == 'web':
                source_desc = f"API ({self.current_web_source or 'auto'})"
            else:
                source_desc = self.current_source_mode
            
            from .database import get_source_selector
            source_selector = get_source_selector()

            for active_pdb_id in pdb_ids:
                self.logger.info(
                    f"Clearing cache and re-collecting motifs for {active_pdb_id} from {source_desc}..."
                )
                if source_selector and hasattr(source_selector, '_cache_manager'):
                    try:
                        source_selector._cache_manager.clear_cache_for_pdb(active_pdb_id)
                        self.logger.debug(f"Cleared cache entries for {active_pdb_id}")
                    except Exception:
                        pass  # Cache clearing is best-effort

                self.fetch_motif_data_action(active_pdb_id)
                self.logger.success(
                    f"Refresh complete for {active_pdb_id} from {source_desc}"
                )

                
        except Exception as e:
            self.logger.error(f"Failed to refresh motifs: {e}")
    
    def print_source_info(self):
        """Print the currently selected data source, loaded PDB, and motif count."""
        print("\n" + "="*70)
        print("  CURRENT SOURCE")
        print("="*70)
        
        # Show loaded PDB info
        pdb_id = self.loaded_pdb_id
        if pdb_id:
            print(f"\n  Loaded PDB: {pdb_id.upper()}")
            # Show motif counts if available
            loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs() if self.viz_manager and self.viz_manager.motif_loader else {}
            if loaded_motifs:
                total_instances = sum(len(info.get('motif_details', [])) for info in loaded_motifs.values())
                print(f"  Motifs: {len(loaded_motifs)} types, {total_instances} total instances")
            else:
                print(f"  Motifs: None loaded (run rmv_load_motif)")
        else:
            print(f"\n  Loaded PDB: None (run rmv_fetch <PDB_ID>)")
        
        # Show chain mode
        cif_mode = getattr(self, 'cif_use_auth', 1)
        chain_label = 'auth_asym_id' if cif_mode == 1 else 'label_asym_id'
        print(f"  Chain ID mode: {chain_label} (cif_use_auth={cif_mode})")
        
        # Determine and display the active source with ID
        print()
        if self.current_source_mode == 'user':
            tool_descriptions = {
                'fr3d': '[5] FR3D (BGSU base pair annotations)',
                'rnamotifscan': '[6] RNAMotifScan (RMS - structural motif search)',
                'rms': '[6] RNAMotifScan (RMS - structural motif search)',
                'rnamotifscanx': '[7] RNAMotifScanX (RMSX - extended motif search)',
                'rmsx': '[7] RNAMotifScanX (RMSX - extended motif search)',
            }
            tool_name = self.current_user_tool or 'unknown'
            description = tool_descriptions.get(tool_name, tool_name)
            print(f"  Source: {description}")
            print(f"  Type: User annotations")
            
        elif self.current_source_mode == 'local':
            if self.current_local_source == 'atlas':
                print(f"  Source: [1] RNA 3D Motif Atlas")
                print(f"  Type: Local (offline) - 759 PDB structures")
            elif self.current_local_source == 'rfam':
                print(f"  Source: [2] Rfam")
                print(f"  Type: Local (offline) - 173 PDB structures")
            else:
                print(f"  Source: [1] RNA 3D Motif Atlas + [2] Rfam")
                print(f"  Type: Local (offline)")
            
        elif self.current_source_mode == 'web':
            if self.current_web_source == 'bgsu_api':
                print(f"  Source: [3] BGSU RNA 3D Hub")
                print(f"  Type: Online API - ~3000+ PDB structures")
            elif self.current_web_source == 'rfam_api':
                print(f"  Source: [4] Rfam API")
                print(f"  Type: Online API - All Rfam motifs")
            else:
                print(f"  Source: Online API (auto-select)")
                print(f"  Type: Online API")
            
        elif self.current_source_mode == 'combine':
            ids_str = ', '.join(str(s) for s in self.combined_source_ids)
            names = []
            for sid in self.combined_source_ids:
                info = SOURCE_ID_MAP.get(sid, {})
                names.append(f"[{sid}] {info.get('name', 'Unknown')}")
            print(f"  Source: Combined - {' + '.join(names)}")
            print(f"  Type: Multi-source merge (IDs: {ids_str})")
            
        else:
            print(f"  Source: None selected")
            print(f"  Run: rmv_db   (check available sources)")
            print(f"  Run: rmv_db <N>    (1-8)")
        
        # Always show workflow steps
        print("\n" + "-"*70)
        print("   WORKFLOW:")
        print("     Step 1: rmv_fetch <PDB_ID>       # Load PDB structure")
        print("     Step 2: rmv_db               # Check available sources")
        print("     Step 3: rmv_db <N>                # Select data source (1-8)")
        print("     Step 4: rmv_load_motif            # Fetch motif data")
        print("-"*70)
        print("   AVAILABLE SOURCES:")
        print("     [1] RNA 3D Motif Atlas   [2] Rfam          (offline)")
        print("     [3] BGSU API       [4] Rfam API      (online)")
        print("     [5] FR3D           [6] RMS   [7] RMSX  [8] NoBIAS (user annotations)")
        print("")
        print("\n" + "="*70 + "\n")


# Global GUI instance
_gui_instance = None
gui = None  # Module-level gui reference (set by initialize_gui())


def get_gui():
    """Get or create global GUI instance."""
    global _gui_instance
    if _gui_instance is None:
        _gui_instance = MotifVisualizerGUI()
    return _gui_instance


def initialize_gui():
    """Initialize GUI and register commands."""
    global gui
    gui = get_gui()
    
    # Register PyMOL commands
    def fetch_raw_pdb(pdb_id='', background_color='', cif_use_auth='', *extra_args, **kwargs):
        """PyMOL command: Load raw PDB structure(s) (motif data auto-loads for
        already-known PDBs; brand-new PDBs still need rmv_db + rmv_load_motif).
        
        Downloads and loads the PDB/mmCIF structure into PyMOL. Accepts a
        single PDB ID or a comma-separated list to load several PDBs in one
        call; each one is loaded from its own last-selected source (or the
        currently active source, if it has no history yet).
        Previously loaded PDB+source motif datasets are preserved for cross-PDB
        superimposition (rmv_super / rmv_align) until rmv_reset is called.
        
        Usage:
            rmv_fetch 1S72                           # Load PDB structure
            rmv_fetch 1S72, 4V9F                     # Load multiple PDBs
            rmv_fetch 1S72, bg_color=lightgray       # With background color
            rmv_fetch 1S72, cif_use_auth=0           # Use label_asym_id chains
        
        Chain ID modes:
            cif_use_auth=1 (default)  - Use auth_asym_id (0, 9, A, B...)
            cif_use_auth=0            - Use label_asym_id (AA, BA, CA...)
        """
        if not pdb_id:
            gui.command_error("Usage: rmv_fetch <PDB_ID>[, <PDB_ID2>, ...] [, bg_color=gray80] [, cif_use_auth=0]")
            gui.logger.error("Examples:")
            gui.logger.error("  rmv_fetch 1S72")
            gui.logger.error("  rmv_fetch 1S72, 4V9F")
            gui.logger.error("  rmv_fetch 1S72, bg_color=lightgray")
            gui.logger.error("  rmv_fetch 1S72, cif_use_auth=0    (use label_asym_id)")
            return
        
        positional = [pdb_id, background_color, cif_use_auth] + list(extra_args)
        positional = [str(value).strip() for value in positional if str(value).strip()]
        pdb_arg = ",".join(positional)
        # Background color is only honored via the explicit bg_color= keyword
        # (or 'bg_color=' embedded below). A bare second positional is always a
        # second PDB token, so a typo like '$v9F' is reported as an invalid PDB
        # rather than being sent to cmd.bg_color.
        bg_arg = str(kwargs.get('bg_color', '') or '').strip() or None
        
        # Handle cif_use_auth parameter - may be embedded in pdb_id or bg_color
        # because PyMOL's cmd.extend doesn't always separate keyword args correctly.
        # User might type: rmv_fetch 1S72 cif_use_auth=0  (space, no comma)
        #              or: rmv_fetch 1S72, cif_use_auth=0  (comma-separated)
        import re
        cif_auth_val = 1  # Default: auth_asym_id
        
        # Extract cif_use_auth= from pdb_id (space-separated case)
        cif_match = re.search(r'\bcif_use_auth\s*=\s*(\S+)', pdb_arg, re.IGNORECASE)
        if cif_match:
            cif_str = cif_match.group(1).strip()
            if cif_str in ('0', 'off', 'false', 'label'):
                cif_auth_val = 0
            pdb_arg = re.sub(r'\s*\bcif_use_auth\s*=\s*\S+', '', pdb_arg, flags=re.IGNORECASE).strip()
        
        # Extract bg_color= from pdb_id (space-separated case)
        bg_match = re.search(r'\bbg_color\s*=\s*(\S+)', pdb_arg, re.IGNORECASE)
        if bg_match:
            bg_arg = bg_match.group(1).strip()
            pdb_arg = re.sub(r'\s*\bbg_color\s*=\s*\S+', '', pdb_arg, flags=re.IGNORECASE).strip()
        
        # Extract cif_use_auth= from background_color (comma-separated positional fallback)
        if bg_arg:
            cif_match_bg = re.search(r'\bcif_use_auth\s*=\s*(\S+)', bg_arg, re.IGNORECASE)
            if cif_match_bg:
                cif_str = cif_match_bg.group(1).strip()
                if cif_str in ('0', 'off', 'false', 'label'):
                    cif_auth_val = 0
                bg_arg = re.sub(r'\s*\bcif_use_auth\s*=\s*\S+', '', bg_arg, flags=re.IGNORECASE).strip()
                if not bg_arg:
                    bg_arg = None
        
        # Also check the explicit keyword argument
        if kwargs.get('cif_use_auth'):
            cif_use_auth = kwargs.get('cif_use_auth')
        if cif_use_auth and str(cif_use_auth).strip().lower() in ('0', '1', 'off', 'on', 'false', 'true', 'label', 'auth'):
            cif_str = str(cif_use_auth).strip()
            if cif_str in ('0', 'off', 'false', 'label'):
                cif_auth_val = 0
            elif cif_str in ('1', 'on', 'true', 'auth'):
                cif_auth_val = 1

        gui.cif_use_auth = cif_auth_val

        # Split into one or more PDB tokens. A single local file path never
        # contains a comma, so this is safe for both PDB IDs and file paths.
        tokens = [t.strip() for t in pdb_arg.split(',') if t.strip()]
        if not tokens:
            gui.command_error("Usage: rmv_fetch <PDB_ID>[, <PDB_ID2>, ...]")
            return
        multi = len(tokens) > 1

        # In a multi-structure fetch, skip PDB IDs already loaded this session so
        # an existing structure is not re-downloaded and its family table
        # re-printed. Use rmv_refresh to force a reload. Single-token fetches and
        # local file paths are never skipped.
        if multi:
            kept = []
            for token in tokens:
                display_id = token.strip().upper()
                is_pdb_id = len(display_id) == 4 and display_id.isalnum()
                already = (
                    is_pdb_id
                    and display_id in gui.loaded_structures
                    and gui.loaded_structures.get(display_id) in cmd.get_object_list()
                )
                if already:
                    gui.logger.info(
                        f"{display_id} is already loaded - skipping.")
                else:
                    kept.append(token)
            if not kept:
                gui.logger.info("All requested structures are already loaded.")
                return
            tokens = kept
            multi = len(tokens) > 1

        def _fetch_one(token: str):
            import os
            # Detect whether the argument is a local file path or a PDB ID
            is_file_path = (
                os.sep in token or
                token.startswith('~') or
                token.startswith('.') or
                token.endswith(('.pdb', '.cif', '.mmcif', '.pdb.gz', '.cif.gz', '.ent', '.ent.gz'))
            )

            if is_file_path:
                expanded = os.path.expanduser(token)
                if not os.path.isfile(expanded):
                    gui.logger.error(f"File not found: {expanded}")
                    return None

                # Derive a structure name from the filename (strip extension)
                base = os.path.basename(expanded)
                structure_name = base.split('.')[0].lower()
                display_id = structure_name.upper()
            else:
                # Validate PDB ID
                if not token or len(token) != 4 or not token.isalnum():
                    gui.command_error(f"Invalid PDB ID: '{token}'")
                    gui.logger.error("PDB ID must be exactly 4 alphanumeric characters (e.g., 1S72)")
                    gui.logger.error("For local files: rmv_fetch /path/to/file.cif")
                    return None
                expanded = None
                structure_name = token.lower()
                display_id = token.upper()

            # Load the structure
            try:
                # Set cif_use_auth before loading
                try:
                    cmd.set("cif_use_auth", cif_auth_val)
                except Exception:
                    pass
                
                # Set background color if specified
                if bg_arg:
                    cmd.bg_color(bg_arg)
                
                # Delete any existing object with the same name
                try:
                    cmd.delete(structure_name)
                except:
                    pass
                
                if expanded:
                    # Local file - use cmd.load
                    cmd.load(expanded, structure_name)
                else:
                    # PDB ID - use cmd.fetch
                    cmd.fetch(token, structure_name)

                # Apply uniform cartoon tube radius to the loaded structure so
                # motif objects (rendered slightly thicker at 0.4) sit clearly on
                # top.  Scoping to the structure name leaves any other loaded
                # objects untouched.
                try:
                    cmd.set('cartoon_nucleic_acid_mode', 4, structure_name, quiet=1)
                    cmd.set('cartoon_tube_radius', 0.37, structure_name, quiet=1)
                except Exception:
                    pass

                # Store loaded PDB info
                gui.loaded_pdb = structure_name
                gui.loaded_pdb_id = display_id
                gui.loaded_structures[display_id] = structure_name
                gui.annotation_tables.setdefault(
                    display_id, ConsolidatedAnnotationTable(gui.jaccard_threshold, merge_enabled=False)
                )
                gui.loaded_structure_path = os.path.abspath(expanded) if expanded else ''
                
                # Set structure_loader fields
                gui.viz_manager.structure_loader.current_structure = structure_name
                gui.viz_manager.structure_loader.current_pdb_id = display_id
                
                # Preserve previously loaded motif datasets across PDB switches.
                # Cross-PDB rmv_super / rmv_align rely on accumulated per-instance
                # metadata (_pdb_id + _source_suffix). Use rmv_reset to clear.
                
                # Build auth->label chain mapping if in label mode
                gui.auth_to_label_map = {}
                if cif_auth_val == 0:
                    gui.auth_to_label_map = gui._build_auth_label_chain_mapping(display_id)
                
                # Report chain ID mode
                chain_mode = "auth_asym_id (default)" if cif_auth_val == 1 else "label_asym_id"
                if not multi:
                    gui.logger.success(f"Loaded structure {display_id} as '{structure_name}'")
                    gui.logger.info(f"Chain ID mode: {chain_mode}")
                else:
                    gui.logger.debug(
                        f"Loaded structure {display_id} as '{structure_name}' "
                        f"with chain mode {chain_mode}"
                    )

                # Multi-PDB session: if this PDB previously had a source
                # selected (via rmv_db), restore it automatically instead of
                # silently keeping whatever source another PDB last selected.
                # Falls back to whatever source is currently active when this
                # PDB has no history yet.
                restore = gui.pdb_source_state.get(display_id)
                if restore and restore.get('mode'):
                    gui.logger.info(
                        f"Restoring previously selected source for {display_id}: "
                        f"rmv_db {restore['mode']}" + (f" {restore['tool']}" if restore.get('tool') else "")
                    )
                    gui.jaccard_threshold = restore.get('jaccard', gui.jaccard_threshold)
                    select_database(restore['mode'], restore.get('tool') or '', '')

                # When loading several PDBs in one rmv_fetch call, also fetch
                # motif data for each right away using whichever source is
                # now active for it (restored above, or already-selected).
                if multi and gui.current_source_mode:
                    gui.logger.info(f"Auto-loading motifs for {display_id} from the active source...")
                    load_motif_data()
                if gui.rmsx_auto_run_on_fetch:
                    gui.logger.info("RNAMotifScanX auto-run is ON: running RNAMotifScanX wrapper...")
                    ok = gui.run_rmsx_wrapper(display_id)
                    if ok:
                        gui.logger.success("RNAMotifScanX motifs loaded automatically after rmv_fetch")
                    else:
                        gui.logger.warning("RNAMotifScanX auto-run after rmv_fetch failed; structure is still loaded")

                return display_id

            except Exception as e:
                gui.logger.error(f"Failed to load {token}: {str(e)}")
                return None

        loaded_ids = []
        for token in tokens:
            result = _fetch_one(token)
            if result:
                loaded_ids.append(result)

        if multi:
            gui.logger.info("")
            if loaded_ids:
                gui.logger.success(f"Loaded {len(loaded_ids)}/{len(tokens)} PDBs: {', '.join(loaded_ids)}")
            else:
                gui.logger.error("No PDBs were loaded successfully.")
    
    def load_motif_data(argument=''):
        """PyMOL command: Fetch motif data from the selected source for the loaded PDB.
        
        Requires:
            1. A PDB structure must be loaded first (rmv_fetch)
            2. A source must be selected (rmv_db)
        
        Usage:
            rmv_load_motif               Fetch motifs from current source
            rmv_load_motif /path/to/dir Source-specific query model directory

        Optional path semantics:
            - Source 5 (FR3D): argument-free. Configure with
              `rmv_db FR3D`, then run `rmv_load_motif`.
            - Source 7 (RNAMotifScanX): path is treated as RMSX motif query model directory.
            - Sources 1-6 except 5: optional path is ignored.
        
        Workflow:
            rmv_fetch 1S72               # Step 1: Load PDB structure
            rmv_db                  # Step 2: Check available sources
            rmv_db 3                     # Step 3: Select BGSU API
            rmv_load_motif               # Step 4: Fetch motif data
            rmv_summary                  # Step 5: Show motif types & counts
            rmv_summary HL               # Step 6: Show instances of a type
            rmv_show HL                  # Step 7: Render hairpin loops
        """
        # Check: PDB must be loaded first
        if not gui.loaded_pdb_id:
            gui.logger.error("No PDB structure loaded.")
            gui.logger.info("Load a structure first:")
            gui.logger.info("  rmv_fetch <PDB_ID>")
            gui.logger.info("  Example: rmv_fetch 1S72")
            return
        
        # Check: Source must be selected
        if not gui.current_source_mode:
            gui.logger.error("No source selected.")
            gui.logger.info("Select a source first:")
            gui.logger.info("  rmv_db <N>                 (1-8)")
            gui.logger.info("  Example: rmv_db 3          (BGSU API)")
            gui.logger.info("  rmv_db                (list all)")
            return
        
        pdb_id = gui.loaded_pdb_id

        # Optional path argument for source-specific query model directories.
        user_path_arg = ''
        arg_text = str(argument or '').strip()
        tokens = []
        if arg_text:
            try:
                tokens = shlex.split(arg_text)
            except Exception:
                tokens = [arg_text]
            if tokens:
                user_path_arg = tokens[0]

        def _resolve_path_with_spaces(raw_tokens):
            if not raw_tokens:
                return '', []
            # Prefer the full argument string first (supports spaces in path).
            candidate_full = os.path.abspath(os.path.expanduser(str(arg_text).strip()))
            if os.path.exists(candidate_full):
                return candidate_full, []
            # Otherwise find the longest token prefix that resolves to an existing path.
            for i in range(len(raw_tokens), 0, -1):
                candidate = ' '.join(raw_tokens[:i])
                expanded = os.path.abspath(os.path.expanduser(candidate))
                if os.path.exists(expanded):
                    return expanded, raw_tokens[i:]
            return '', raw_tokens[1:]

        # Dispatch to appropriate loader based on source mode
        if gui.current_source_mode == 'user' and gui.current_user_tool:
            active_source_id = None
            try:
                active_source_id = int(gui.current_source_id)
            except Exception:
                active_source_id = None

            # Source 5 (FR3D) is fully config-driven and argument-free.
            # It is registered via `rmv_db FR3D` and run with
            # a bare `rmv_load_motif`.
            if active_source_id == 5:
                if user_path_arg:
                    gui.logger.warning(
                        "Source 5 (FR3D) does not accept rmv_load_motif path arguments; "
                        "register it with 'rmv_db FR3D' instead."
                    )
                if getattr(gui, 'fr3d_data_mode', '') == 'cache':
                    gui.load_user_annotations_action(gui.current_user_tool, pdb_id, auto_pipeline=False)
                    return
                gui.run_fr3d_search(pdb_id)
                return

            if user_path_arg:
                path_override, remaining_tokens = _resolve_path_with_spaces(tokens)
                if not path_override:
                    path_override = os.path.abspath(os.path.expanduser(user_path_arg))
                if active_source_id == 7:
                    if not os.path.isdir(path_override):
                        gui.logger.error(f"RNAMotifScanX query model directory not found: {path_override}")
                        return
                    gui.load_user_annotations_action(
                        gui.current_user_tool,
                        pdb_id,
                        rmsx_query_models_dir=path_override,
                    )
                    return

                # Sources 1-6 except 5 ignore optional rmv_load_motif path per policy.
                if active_source_id in [1, 2, 3, 4, 6]:
                    gui.logger.info(
                        f"Ignoring path argument for source {active_source_id}; this source does not use rmv_load_motif path overrides."
                    )
                else:
                    gui.logger.info(
                        "Ignoring path argument for current source; only source 5 (FR3D) and source 7 (RNAMotifScanX) use it."
                    )

            gui.load_user_annotations_action(gui.current_user_tool, pdb_id)
        else:
            if user_path_arg:
                gui.logger.info(
                    "Ignoring path argument for current source mode; only source 5 (FR3D) and source 7 (RNAMotifScanX) use it."
                )
            gui.fetch_motif_data_action(pdb_id, None)
    
    def load_structure(pdb_id_or_path='', background_color='', database=''):
        """PyMOL command: Load structure and automatically show all motifs.
        
        NOTE: This command is deprecated in the recommended workflow.
        Users should use rmv_fetch first, then rmv_db; rmv_load_motif.
        
        Usage:
            rmv_load <pdb_id_or_path>
        """
        if not pdb_id_or_path:
            gui.command_error("Usage: rmv_load <PDB_ID>")
            return
        
        pdb_arg = str(pdb_id_or_path).strip().upper()
        
        # Instead of trying to fetch + load motifs all at once (which hangs PyMOL),
        # guide the user through the proper step-by-step workflow.
        print("\n" + "=" * 60)
        print("  RECOMMENDED WORKFLOW")
        print("=" * 60)
        print(f"\n  To visualize motifs for {pdb_arg}, follow these steps:\n")
        print(f"  Step 1:  rmv_fetch {pdb_arg}        # Fetch the PDB structure")
        print(f"  Step 2:  rmv_db               # Check available data sources")
        print(f"  Step 3:  rmv_db <N>                # Select data source (1-8)")
        print(f"  Step 4:  rmv_load_motif             # Fetch motif data")
        print(f"  Step 5:  rmv_summary               # Show motif types & counts")
        print(f"  Step 6:  rmv_summary <TYPE>         # Show instances of a type")
        print(f"  Step 7:  rmv_show <TYPE>            # Render a motif type")
        print(f"  Step 8:  rmv_show <TYPE> <NO>       # Zoom to specific instance")
        print(f"\n  Type 'rmv_db' to see all available data sources.")
        print(f"  Type 'rmv_help' for full command reference.")
        print("=" * 60 + "\n")
    
    def toggle_motif(motif_type='', visible=''):
        """PyMOL command: Toggle motif visibility."""
        # PyMOL can pass arguments different ways, so handle both
        
        # Case 1: Both arguments passed separately
        if motif_type and visible:
            motif_arg = motif_type
            visible_arg = visible
        else:
            # Case 2: Everything in motif_type as a single string
            full_arg = str(motif_type).strip()
            parts = full_arg.split()
            
            if len(parts) < 2:
                gui.command_error(f"Usage: rmv_toggle MOTIF_TYPE on/off")
                gui.logger.error(f"Example: rmv_toggle HL on")
                return
            
            motif_arg = parts[0]
            visible_arg = parts[1]
        
        # Parse visibility
        visible_bool = str(visible_arg).lower() in ['on', 'true', '1', 'yes', 'show']
        motif_arg = str(motif_arg).upper().strip()
        
        gui.toggle_motif_action(motif_arg, visible_bool)
    
    def list_sources():
        """PyMOL command: Show available data sources."""
        gui.print_sources()

    def set_debug(mode=''):
        """PyMOL command: Enable or disable diagnostic debug messages."""
        value = str(mode or '').strip().lower()
        if value not in ('on', 'off'):
            gui.command_error("Usage: rmv_debug ON|OFF")
            return
        gui.logger.set_debug(value == 'on')
        gui.logger.info(f"Debug messages: {'ON' if value == 'on' else 'OFF'}")

    def hide_annotation(target='', *extra_args, **_kwargs):
        """PyMOL command: dehighlight a stable motif ID/group or all motifs."""
        parts = [str(target).strip()] + [str(part).strip() for part in extra_args]
        value = " ".join(part for part in parts if part)
        if not value or value.lower() == 'all':
            structure = gui.loaded_pdb or gui.loaded_pdb_id or ''
            if structure:
                cmd.color('gray80', f"model {structure}")
                gui.logger.info(f"Dehighlighted all motif residues on {structure}.")
            else:
                gui.command_error("No structure loaded. Use rmv_fetch first.")
            return
        if value in gui.query_groups or any(
            value in table._rows for table in gui.annotation_tables.values()
        ):
            gui.view_annotation_results(value, hide=True)
            return
        gui.command_error(
            f"Unknown hide target '{value}'. Use rmv_hide <Motif_ID|group> or rmv_hide all."
        )
    
    def show_help():
        """PyMOL command: Show all available commands."""
        gui.print_help()
    
    def set_bg_color(color_name='gray80'):
        """PyMOL command: Change background color of non-motif residues."""
        color_arg = str(color_name).strip()
        if not color_arg:
            color_arg = 'gray80'
        gui.set_background_color(color_arg)
    
    def motif_summary(motif_type='', *extra_args, **_kwargs):
        """PyMOL command: Show motif summary table (console only, no rendering).
        
        Usage:
            rmv_summary              Show all motifs summary for loaded PDB
            rmv_summary HL           Show detailed instances of HL motif
            rmv_summary HL 1         Show specific HL instance #1
            rmv_summary HL 1-5       Show instances 1 through 5
            rmv_summary HL 1,2,3     Show instances 1, 2, and 3
            rmv_summary K-TURN db 3 and db 2      Hierarchy table, source filter
            rmv_summary SR db 5 or db 3           Hierarchy table, source filter
            rmv_summary SR db 4 and not db 3       Hierarchy table, source filter
            rmv_summary SR db 3 and db 7, as alias      Freeze result under 'alias'
        """
        all_parts = (str(motif_type).split() if motif_type else []) + [str(a) for a in extra_args]
        all_parts = [p.strip() for p in all_parts if p.strip()]

        if not all_parts:
            gui.print_motif_summary()
            return

        all_parts, alias, alias_err = _extract_trailing_alias(all_parts)
        if alias_err:
            gui.logger.error(alias_err)
            return

        motif_filter, ops = _split_motif_and_db_expression(all_parts)
        ops = ops or []

        # Strip trailing instance-selector tokens ('1', '1-5', '1,2,3') from
        # whatever comes before the 'db' expression (or from all tokens if
        # there is no 'db' expression at all).
        filter_parts = motif_filter.split()
        instance_nums: List[int] = []
        while filter_parts:
            expanded = _expand_instance_token(filter_parts[-1])
            if expanded is None:
                break
            instance_nums = expanded + instance_nums
            filter_parts.pop()
        motif_arg = ' '.join(filter_parts).upper()

        if not motif_arg and not ops:
            gui.print_motif_summary()
            return

        if instance_nums:
            seen_nums = list(dict.fromkeys(instance_nums))
            for inst_no in seen_nums:
                gui.show_motif_instance_summary(motif_arg, inst_no)
        elif ops:
            gui.print_source_query_table(motif_arg, ops, alias=alias)
        else:
            gui.show_motif_summary_for_type(motif_arg, alias=alias)

    def select_annotation(query='', *extra_args, **_kwargs):
        """PyMOL command: save a motif/source/structure query as a group."""
        parts = [str(query).strip()] + [str(part).strip() for part in extra_args]
        text = ", ".join(part for part in parts if part)
        if not text:
            gui.command_error(
                "Usage: rmv_select <motif>, <structures>, <sources>, as <group>"
            )
            return
        gui.select_annotation_query(text)

    def list_annotations(target='', *extra_args, **_kwargs):
        """PyMOL command: list consolidated rows by ID or saved group."""
        parts = [str(target).strip()] + [str(part).strip() for part in extra_args]
        gui.list_annotation_results(" ".join(part for part in parts if part))

    def create_annotation_object(target='', *extra_args, **_kwargs):
        """PyMOL command: create selectable objects for an ID or group."""
        parts = [str(target).strip()] + [str(part).strip() for part in extra_args]
        value = " ".join(part for part in parts if part)
        if not value:
            gui.command_error("Usage: rmv_create_object <Motif_ID|group_name>")
            return
        gui.create_annotation_objects(value)
    
    def select_database(mode='', tool='', jaccard_threshold=''):
        """Select named annotation sources and load them for the active structure.

        Usage:
            rmv_db RNA3DMotifAtlas
            rmv_db Rfam
            rmv_db RNA3DMotifAtlas,Rfam
            rmv_db FR3D,RNAMotifScanX
        """
        if not mode:
            gui.print_sources()
            return

        if jaccard_threshold:
            gui.logger.warning("Jaccard threshold is code-defined in rsmviewer/database/consolidated_table.py; command-line overrides are ignored.")
        
        source_parts = [str(mode).strip()]
        if tool:
            source_parts.append(str(tool).strip())
        source_expression = ",".join(source_parts)

        from .database.source_registry import get_source_registry
        try:
            source_names = get_source_registry().parse_names(source_expression)
        except ValueError as exc:
            gui.logger.error(str(exc))
            return

        if not gui.loaded_pdb_id:
            gui.logger.error("No structure loaded. Run 'rmv_fetch <PDB_ID>' before rmv_db.")
            return

        # Existing loaders are migrated incrementally; numeric IDs remain an
        # internal bridge and are never accepted or displayed by rmv_db.
        internal_source_ids = {
            "RNA3DMotifAtlas": 3,
            "Rfam": 4,
            "FR3D": 5,
            "RNAMotifScanX": 7,
        }
        selected_ids = [internal_source_ids[name] for name in source_names]

        # Remember this exact selection for the currently loaded PDB so a
        # later 'rmv_fetch' back to it can auto-restore the same source
        # (multi-PDB session workflow).
        gui.pdb_source_state[gui.loaded_pdb_id] = {
            'mode': ",".join(source_names),
            'tool': None,
            'jaccard': gui.jaccard_threshold,
        }
        gui.current_source_names = list(source_names)

        if len(selected_ids) == 1:
            gui._handle_source_by_id(selected_ids[0])
        else:
            gui._handle_multi_source(selected_ids)

        gui.logger.info(f"Loading annotations from: {', '.join(source_names)}")
        structure_ids = list(gui.loaded_structures) or [gui.loaded_pdb_id]
        active_structure_id = gui.loaded_pdb_id
        active_structure_name = gui.loaded_pdb
        for structure_id in structure_ids:
            if not structure_id:
                continue
            gui.loaded_pdb_id = structure_id
            gui.loaded_pdb = gui.loaded_structures.get(structure_id, structure_id)
            gui.viz_manager.structure_loader.current_pdb_id = structure_id
            gui.viz_manager.structure_loader.current_structure = gui.loaded_pdb
            gui.logger.info(f"Loading annotations for {structure_id}...")
            load_motif_data()

        gui.loaded_pdb_id = active_structure_id
        gui.loaded_pdb = active_structure_name
        if active_structure_id:
            gui.viz_manager.structure_loader.current_pdb_id = active_structure_id
            gui.viz_manager.structure_loader.current_structure = active_structure_name
    
    def set_source(mode='', source_id=''):
        """PyMOL command: Show current source info or detailed info about a specific source.
        
        Usage:
            rmv_source info          - Show currently selected source info
            rmv_source info <N>      - Show detailed info about source N (1-8)
        """
        if not mode:
            gui.command_error("Usage: rmv_source info [<ID>]")
            gui.logger.error("  rmv_source info        Show current source info")
            gui.logger.error("  rmv_source info <N>    Show detailed info about source N")
            gui.logger.error("  rmv_db <N>             Select a source")
            return
        
        mode_arg = str(mode).strip()
        tool_arg = str(source_id).strip() if source_id else None
        
        # Handle PyMOL passing arguments as combined string
        parts = mode_arg.split(None, 1)
        first_part = parts[0].lower()
        remaining_arg = parts[1] if len(parts) > 1 else tool_arg
        
        if first_part == 'info':
            gui._handle_source_info_command(remaining_arg)
            return
        
        gui.command_error(f"Unknown subcommand: {first_part}")
        gui.logger.info("Use: rmv_source info [<ID>]")
        gui.logger.error("       rmv_db <ID>        - Select a source")
        gui.logger.error("  rmv_source info <N>  - N can be 1-8")
    
    def refresh_motifs(pdb_id=''):
        """PyMOL command: Force refresh cache and collect motif data again.
        
        Clears cached data for all active loaded PDBs and re-fetches motif
        information from the last selected source (or combined sources if
        combine mode was used). A PDB argument limits the refresh to that PDB.
        
        Usage:
            rmv_refresh        - Refresh all active PDBs from the last selected source
            rmv_refresh <PDB_ID> - Refresh one active PDB
        """
        pdb_arg = str(pdb_id).strip() if pdb_id else None
        gui.refresh_motifs_action(pdb_arg)
    
    def _resolve_motif_type_and_instance(full_arg, instance_arg=''):
        """Resolve multi-word motif type name and optional instance number.
        
        Handles cases like:
            '4-WAY JUNCTION (J4)', ''    -> ('4-WAY JUNCTION (J4)', None)
            '4-WAY JUNCTION (J4)', '1'   -> ('4-WAY JUNCTION (J4)', 1)
            '4-WAY JUNCTION (J4) 1', ''  -> ('4-WAY JUNCTION (J4)', 1)
            'HL', '1'                    -> ('HL', 1)
            'HL 1', ''                   -> ('HL', 1)
        """
        full_arg = str(full_arg).strip().upper()
        instance_arg = str(instance_arg).strip() if instance_arg else ''
        
        # If instance_arg is provided and is a number, use it directly
        if instance_arg:
            try:
                return full_arg, int(instance_arg)
            except ValueError:
                # instance_arg is actually part of the motif name
                full_arg = f"{full_arg} {instance_arg}"
        
        # Try to match against loaded motif types
        loaded_motifs = gui.viz_manager.motif_loader.get_loaded_motifs() if gui.viz_manager.motif_loader else {}
        
        # Check if full_arg exactly matches a loaded motif type
        if full_arg in loaded_motifs:
            return full_arg, None

        # Resolve user aliases (HL/IL/Jn) to loaded semantic keys when possible
        resolved_full = gui._resolve_loaded_motif_type(full_arg, loaded_motifs)
        if resolved_full in loaded_motifs:
            return resolved_full, None
        
        # Check if the last token is a number (instance ID)
        # Try removing the last word and see if the rest matches a motif type
        parts = full_arg.rsplit(None, 1)  # split from right, max 1 split
        if len(parts) == 2 and parts[1].isdigit():
            candidate_type = parts[0]
            instance_no = int(parts[1])
            resolved_candidate = gui._resolve_loaded_motif_type(candidate_type, loaded_motifs)
            if resolved_candidate in loaded_motifs:
                return resolved_candidate, instance_no
            if candidate_type in loaded_motifs:
                return candidate_type, instance_no
            # Also try without matching - maybe it's a simple type like 'HL 1'
            return candidate_type, instance_no
        
        # No instance number found
        return full_arg, None
    
    def show_motif(motif_type='', *extra_args, **_kwargs):
        """PyMOL command: Show specific motif type, all types, or specific instance.
        
        Usage:
            rmv_show ALL           - Show all loaded motif types (creates objects)
            rmv_show GNRA          - Show only GNRA motifs (all instances)
            rmv_show HL            - Show only hairpin loops (all instances)
            rmv_show HL 1          - Show specific HL instance #1 (zoom + details)
            rmv_show HL 1,3,5     - Show specific HL instances 1, 3, 5
            rmv_show GNRA 2        - Show specific GNRA instance #2
            rmv_show K-TURN nobias - Show K-TURN instances unique to NoBIAS (combine mode)
            rmv_show K-TURN rmsx   - Show K-TURN instances unique to RMSX (combine mode)
            rmv_show K-TURN shared - Show K-TURN instances found in both sources
            rmv_show 4-WAY JUNCTION (J4)      - Multi-word motif type
            rmv_show 4-WAY JUNCTION (J4) 1    - Multi-word with instance
            rmv_show K-TURN, SARCIN-RICIN, HL - Multiple families in one call
            rmv_show SR db 3 and db 7, as alias - Freeze + group under 'alias'
        """
        # PyMOL splits by both spaces and commas into separate positional
        # args.  e.g. "rmv_show HL 1,3,5" arrives as ("HL", "1", "3", "5")
        # However, multi-token args like "rmv_show K-TURN rmsx" can arrive
        # as a single positional ("K-TURN rmsx") depending on quoting and
        # the PyMOL command parser path.  Split the first arg by spaces
        # ourselves so that source filter words and multi-word motif types
        # are always recognised.
        first_parts = str(motif_type).split() if motif_type else []
        rest_parts = [str(a) for a in extra_args]
        all_parts = first_parts + rest_parts
        # Strip empty / whitespace-only parts
        all_parts = [p.strip() for p in all_parts if p.strip()]
        
        # Extract padding=N if present (positional: "padding=10" or kwarg: ", padding=10")
        padding = 0
        filtered_parts = []
        for p in all_parts:
            if p.lower().startswith('padding='):
                try:
                    padding = int(p.split('=', 1)[1])
                except ValueError:
                    gui.logger.error(f"Invalid padding value: {p}")
                    return
            else:
                filtered_parts.append(p)
        all_parts = filtered_parts
        # Also check keyword args (handles PyMOL comma syntax: "rmv_show K-TURN 1, padding=10")
        if padding == 0 and 'padding' in _kwargs:
            try:
                padding = int(_kwargs['padding'])
            except (ValueError, TypeError):
                gui.logger.error(f"Invalid padding value: {_kwargs['padding']}")
                return

        all_parts, alias, alias_err = _extract_trailing_alias(all_parts)
        if alias_err:
            gui.logger.error(alias_err)
            return
        
        # "rmv_show [MOTIF] db N [and/or/not db M ...]" - residue-based query syntax
        motif_filter, db_ops = _split_motif_and_db_expression(all_parts)
        if db_ops is not None or alias:
            gui.show_source_query(motif_filter, db_ops or [], padding=padding, alias=alias)
            return
        
        if not all_parts:
            print("\n  [rmv_show] Render motif objects")
            print("\n  Usage:")
            print("    rmv_show ALL")
            print("    rmv_show <MOTIF_TYPE>")
            print("    rmv_show <MOTIF_TYPE> <INSTANCE_NO>")
            print("    rmv_show <MOTIF_TYPE> 1,3,5")
            print("\n  Examples:")
            print("    rmv_show HL")
            print("    rmv_show HL 1")
            print("    rmv_show HL 1,3,5")
            print("    rmv_show K-TURN shared   (combine mode source attribution)")
            print("\n  Stage checks:")
            if not gui.loaded_pdb_id:
                print("    ERROR: No structure loaded. Run: rmv_fetch <PDB_ID>")
            loaded = gui.viz_manager.motif_loader.get_loaded_motifs() if gui.viz_manager and gui.viz_manager.motif_loader else {}
            if not loaded:
                print("    ERROR: No motifs loaded. Run: rmv_load_motif")
            else:
                print(f"    OK: {len(loaded)} motif types are loaded for display")
            return
        
        # Handle 'ALL' keyword
        if all_parts[0].strip().upper() == 'ALL':
            gui.viz_manager.show_all_motifs(
                filter_pdb=gui.loaded_pdb_id or '',
                filter_suffix=gui._get_source_suffix())
            return

        # Application 6: bare alias name(s), e.g. "rmv_show group_1" or
        # "rmv_show group_1, group_2, group_3" - redisplay (or recreate)
        # previously saved 'as ALIAS' groups. Only triggers when EVERY
        # token is a known alias, so it never shadows real motif types.
        from .database.motif_hierarchy_cache import get_hierarchy_cache as _ghc
        _alias_cache = _ghc()
        if all_parts and all(_alias_cache.alias_exists(p) for p in all_parts):
            for a in all_parts:
                gui._show_alias(a, padding=padding)
            return

        # Application 3: multiple motif families in one call, e.g.
        # "rmv_show K-TURN, SARCIN-RICIN, HL" (commas already collapsed to
        # plain tokens by PyMOL). Only triggers when the WHOLE token list
        # cleanly splits into 2+ distinct known types with nothing left
        # over, so single multi-word types and instance-number lists are
        # never affected.
        loaded_motifs_for_split = gui.viz_manager.motif_loader.get_loaded_motifs() if gui.viz_manager and gui.viz_manager.motif_loader else {}
        multi_types = _partition_into_known_types(all_parts, loaded_motifs_for_split, gui._resolve_loaded_motif_type)
        if multi_types:
            for one_type in multi_types:
                show_motif(one_type, padding=padding)
            return
        
        # Separate trailing numeric parts (instance numbers) from motif name.
        # Walk from the end: pure-digit parts are instance numbers.
        instance_nums = []
        while all_parts and all_parts[-1].isdigit():
            instance_nums.insert(0, int(all_parts.pop()))
        
        if not all_parts:
            gui.command_error("Usage: rmv_show <MOTIF_TYPE> [<INSTANCE_NO>]")
            return

        # --- Source filter detection (combine mode) ---
        # e.g. "rmv_show K-TURN nobias", "rmv_show K-TURN shared", or a
        # full source name "rmv_show K-TURN BGSU RNA 3D Hub".  We try
        # matching progressively LONGER suffixes of all_parts as the
        # source filter (longest first), joining them with spaces.  The
        # longest match that resolves wins, so multi-word source names
        # always take priority over a coincidental short-alias match
        # within the same token sequence.
        source_filter_ids = None
        source_filter_word = None
        if not instance_nums and len(all_parts) >= 2:
            for k in range(len(all_parts) - 1, 0, -1):
                candidate = " ".join(all_parts[-k:])
                candidate_motif = " ".join(all_parts[:-k]).upper()
                if not candidate_motif:
                    continue
                resolved = gui._resolve_source_filter(candidate_motif, candidate)
                if resolved is not None:
                    source_filter_ids = resolved
                    source_filter_word = candidate
                    all_parts = all_parts[:-k]  # Strip filter tokens
                    break
        
        raw_motif = " ".join(all_parts)
        motif_arg, inst_no = _resolve_motif_type_and_instance(raw_motif, '')
        
        # If _resolve found an instance number embedded in the motif arg
        # (e.g. "K-TURN 1" from PyMOL comma-split of "K-TURN 1,2")
        if inst_no is not None and inst_no not in instance_nums:
            instance_nums.insert(0, inst_no)

        # If a source filter was detected, use those IDs as instance_nums
        if source_filter_ids is not None:
            if not source_filter_ids:
                gui.logger.info(
                    f"No instances of {motif_arg} unique to "
                    f"'{source_filter_word}' in the combined result.")
                return
            instance_nums = source_filter_ids
        
        # Source filter params - restrict display to current PDB + source
        fpdb = gui.loaded_pdb_id or ''
        fsuf = gui._get_source_suffix()
        
        if instance_nums:
            for inum in instance_nums:
                gui.viz_manager.show_motif_instance(
                    motif_arg, inum,
                    filter_pdb=fpdb, filter_suffix=fsuf,
                    padding=padding)
        else:
            gui.viz_manager.show_motif_type(
                motif_arg,
                filter_pdb=fpdb, filter_suffix=fsuf,
                padding=padding)
    
    def view_motif(motif_type='', *extra_args, **_kwargs):
        """PyMOL command: Zoom to motif regions on the base structure (no objects).

        Usage:
            rmv_view all                      Highlight ALL motif regions on structure
            rmv_view K-TURN                   Zoom to all K-TURN instances
            rmv_view K-TURN 1                 Zoom to instance #1 and create selection
            rmv_view K-TURN, color=red       Highlight K-TURN in red for this view only
            rmv_view hide                     Reset all view coloring to gray
            rmv_view K-TURN hide              Reset only K-TURN view coloring
            rmv_view K-TURN, SARCIN-RICIN, HL Highlight multiple families at once
            rmv_view SR db 3 and db 7            Hierarchy table query, in-place
            rmv_view SR db 3 and db 7, as alias   Freeze + color as 'alias'
        """
        # Build parts list - split every arg by spaces for robustness
        # (PyMOL may pass 'group_SR color=red' or 'K-TURN hide' as one string).
        first_parts = str(motif_type).split() if motif_type else []
        rest_parts = []
        for a in extra_args:
            rest_parts.extend(str(a).split())
        all_parts = first_parts + rest_parts
        all_parts = [p.strip(" ,") for p in all_parts if p.strip(" ,")]

        color_override = _kwargs.get('color')
        if color_override is None:
            color_tokens = [p for p in all_parts if isinstance(p, str) and p.lower().startswith('color=')]
            if color_tokens:
                color_override = color_tokens[-1].split('=', 1)[1].strip()
                all_parts = [p for p in all_parts if p not in color_tokens]

        padding = 0
        # PyMOL may pass 'padding=5' either as a keyword argument or as a
        # positional token, depending on comma placement; handle both.
        if _kwargs.get('padding') is not None:
            try:
                padding = int(_kwargs.get('padding'))
            except (TypeError, ValueError):
                gui.logger.error(
                    f"Invalid padding value '{_kwargs.get('padding')}'. Use padding=<integer>, e.g. padding=5.")
                return
        padding_tokens = [p for p in all_parts if p.lower().startswith('padding=')]
        if padding_tokens:
            try:
                padding = int(padding_tokens[-1].split('=', 1)[1])
                all_parts = [p for p in all_parts if p not in padding_tokens]
            except ValueError:
                gui.logger.error(
                    f"Invalid padding value '{padding_tokens[-1]}'. Use padding=<integer>, e.g. padding=5.")
                return

        all_parts, alias, alias_err = _extract_trailing_alias(all_parts)
        if alias_err:
            gui.logger.error(alias_err)
            return

        def _is_view_target(name):
            if name in gui.query_groups:
                return True
            if any(name in table._rows for table in gui.annotation_tables.values()):
                return True
            # Individual merged row inside a saved group (e.g. group_SR_001).
            return any(
                row.motif_id == name
                for saved in gui.query_groups.values()
                for row in saved.get("rows", [])
            )

        # 'rmv_view <target> <color>' positional color, but only when the 2nd
        # token is not itself a motif ID / group and not 'hide'.
        if (len(all_parts) == 2 and not color_override
                and all_parts[1].lower() != 'hide'
                and _is_view_target(all_parts[0])
                and not _is_view_target(all_parts[1])):
            color_override = all_parts[1]
            all_parts = all_parts[:1]

        # Trailing 'hide' after one or more motif IDs / groups.
        hide_target = (
            len(all_parts) >= 2
            and all_parts[-1].lower() == 'hide'
            and all(_is_view_target(p) for p in all_parts[:-1])
        )
        if hide_target:
            all_parts = all_parts[:-1]

        # Route stable motif IDs and saved groups (one or many, mixed) to the
        # consolidated-table highlighter. When several targets are given and
        # some don't exist (e.g. a group that was never created because its
        # motif matched nothing), highlight the valid ones and warn about the
        # rest instead of failing the whole command.
        if all_parts and any(_is_view_target(p) for p in all_parts):
            valid_targets = [p for p in all_parts if _is_view_target(p)]
            missing_targets = [p for p in all_parts if not _is_view_target(p)]
            for name in missing_targets:
                gui.logger.warning(
                    f"'{name}' is not a saved group or motif ID; skipping it.")
            # Gray every involved base structure once so highlighting a later
            # target does not erase the colors of an earlier one.
            if not hide_target:
                involved = {
                    gui.loaded_structures.get(row.structure_id, row.structure_id)
                    for name in valid_targets
                    for row in gui._rows_for_target(name)
                }
                gui._gray_out_base_structures(involved)
            for index, name in enumerate(valid_targets):
                gui.view_annotation_results(
                    name,
                    color_override=color_override,
                    hide=hide_target,
                    padding=padding,
                    gray_base=False,
                )
            return

        # "rmv_view [MOTIF] db N [and/or/not db M ...]" - residue-based query syntax
        motif_filter, db_ops = _split_motif_and_db_expression(all_parts)
        if db_ops is not None or alias:
            gui.view_source_query(motif_filter, db_ops or [], alias=alias, color_override=color_override)
            return

        if not all_parts:
            print("\n  [rmv_view] Highlight motif residues on the loaded structure (no objects created)")
            print("\n  Usage:")
            print("    rmv_view <MOTIF_ID>")
            print("    rmv_view <GROUP>")
            print("    rmv_view <ID_OR_GROUP>[, <ID_OR_GROUP> ...]")
            print("    rmv_view <TARGET>, color=<name>")
            print("    rmv_view <TARGET>, padding=<n>")
            print("    rmv_view <TARGET> hide")
            print("\n  Examples:")
            print("    rmv_view 1S72_00016")
            print("    rmv_view group_SR")
            print("    rmv_view 1S72_00016, 1S72_00019, color=red")
            print("    rmv_view group_SR, color=red, padding=5")
            print("    rmv_view group_SR hide")
            print("\n  Related: rmv_list (motif IDs) | rmv_colors (color names) | rmv_create_object | rmv_super")
            if not gui.loaded_pdb_id:
                print("\n  Note: no structure loaded yet. Run: rmv_fetch <PDB_ID>  then  rmv_db <SOURCE>")
            elif not gui.annotation_tables:
                print("\n  Note: no annotations loaded yet. Run: rmv_db <SOURCE>")
            return

        # Handle 'rmv_view hide' - reset ALL view coloring
        if len(all_parts) == 1 and all_parts[0].upper() == 'HIDE':
            fpdb = gui.loaded_pdb_id or ''
            fsuf = gui._get_source_suffix()
            gui.viz_manager.reset_view_coloring(
                filter_pdb=fpdb, filter_suffix=fsuf)
            return

        # Handle 'rmv_view K-TURN hide' or 'rmv_view K-TURN 1 hide'
        # Also handle 'rmv_view all hide' as a full reset
        if len(all_parts) >= 2 and all_parts[-1].upper() == 'HIDE':
            motif_parts = [p for p in all_parts[:-1] if not p.isdigit()]
            raw_motif = " ".join(motif_parts)
            fpdb = gui.loaded_pdb_id or ''
            fsuf = gui._get_source_suffix()
            # 'rmv_view all hide' = reset everything
            if raw_motif.upper() in ('ALL', 'MOTIF'):
                gui.viz_manager.reset_view_coloring(
                    filter_pdb=fpdb, filter_suffix=fsuf)
            else:
                motif_arg, _ = _resolve_motif_type_and_instance(raw_motif, '')
                gui.viz_manager.reset_view_coloring(
                    motif_arg, filter_pdb=fpdb, filter_suffix=fsuf)
            return

        # Handle 'rmv_view all' - highlight all motif regions on structure
        if len(all_parts) == 1 and all_parts[0].upper() in ('ALL', 'MOTIF'):
            structure_name = None
            if gui.viz_manager and gui.viz_manager.structure_loader:
                structure_name = gui.viz_manager.structure_loader.get_current_structure()
            if not structure_name:
                gui.logger.error("No structure loaded. Use rmv_fetch first.")
                return
            gui._auto_color_motifs_on_structure(structure_name)
            return

        # Application 3: multiple motif families in one call (see show_motif
        # for the matching heuristic and why commas can't be relied on here).
        loaded_motifs_for_split = gui.viz_manager.motif_loader.get_loaded_motifs() if gui.viz_manager and gui.viz_manager.motif_loader else {}
        multi_types = _partition_into_known_types(all_parts, loaded_motifs_for_split, gui._resolve_loaded_motif_type)
        if multi_types:
            for one_type in multi_types:
                view_motif(one_type, color=color_override)
            return

        # Separate trailing numeric parts
        instance_nums = []
        while all_parts and all_parts[-1].isdigit():
            instance_nums.insert(0, int(all_parts.pop()))

        if not all_parts:
            print("\n  [rmv_view] Usage: rmv_view all | rmv_view <TYPE> [<NO>] | rmv_view hide")
            return

        raw_motif = " ".join(all_parts)
        motif_arg, inst_no = _resolve_motif_type_and_instance(raw_motif, '')
        if inst_no is not None and inst_no not in instance_nums:
            instance_nums.insert(0, inst_no)

        fpdb = gui.loaded_pdb_id or ''
        fsuf = gui._get_source_suffix()

        if instance_nums:
            for inum in instance_nums:
                gui.viz_manager.view_motif_instance(
                    motif_arg, inum,
                    filter_pdb=fpdb, filter_suffix=fsuf,
                    color_spec=color_override)
        else:
            gui.viz_manager.view_motif_type(
                motif_arg,
                filter_pdb=fpdb, filter_suffix=fsuf,
                color_spec=color_override)

    def combine_sets(name_list='', *extra_args, **_kwargs):
        """PyMOL command: Combine saved aliases and/or PyMOL objects into a new alias (Application 6).

        Usage:
            rmv_combine_groups novel_SR known_SR, as group_SR
            rmv_combine_groups obj_1 obj_2, as group_1

        Each name must be either a previously saved alias (rmv_summary/
        rmv_show/rmv_view '... as ALIAS') or an existing PyMOL object name.
        The new alias can then be superimposed directly: rmv_super group_SR
        """
        first_parts = str(name_list).split() if name_list else []
        rest_parts = []
        for a in extra_args:
            rest_parts.extend(str(a).split())
        all_parts = first_parts + rest_parts
        all_parts = [p.strip().strip(',') for p in all_parts if p.strip().strip(',')]

        all_parts, new_alias, alias_err = _extract_trailing_alias(all_parts)
        if alias_err or not new_alias or not all_parts:
            gui.command_error("Usage: rmv_combine_groups <group>, <group>[, ...], as <NEW_GROUP>")
            return

        # rmv_combine_groups over rmv_select groups: this is where cross-family
        # containment + Jaccard merging happens (professor's design). Fully
        # contained annotations from a different family are removed here, but
        # the original groups and the raw loaded annotations stay intact.
        if all(part in gui.query_groups for part in all_parts):
            if new_alias in gui.query_groups:
                gui.logger.error(f"Group '{new_alias}' already exists. Choose a different name.")
                return

            from collections import defaultdict

            # struct -> group_name -> [rows]
            per_struct_group: Dict[str, Dict[str, list]] = defaultdict(lambda: defaultdict(list))
            for part in all_parts:
                for row in gui._rows_for_target(part):
                    per_struct_group[row.structure_id][part].append(row)

            if not per_struct_group:
                gui.logger.error("Nothing to combine - the selected groups are empty.")
                return

            input_total = sum(
                len(rows) for groups in per_struct_group.values() for rows in groups.values()
            )

            # Merge across groups while preserving each database's own columns.
            # The surviving rows carry per-source labels (RNA3DMotifAtlas, Rfam,
            # FR3D, RNAMotifScanX), not the input group names; the group origin
            # is tracked separately for per-group coloring.
            merged_pairs = []  # (AnnotationRow, [contributing_group_names])
            for structure_id, group_rows in per_struct_group.items():
                ordered = [(part, group_rows[part]) for part in all_parts if part in group_rows]
                merged_pairs.extend(
                    gui._combine_group_rows(structure_id, ordered, new_alias)
                )

            if not merged_pairs:
                gui.logger.error("Nothing to combine - no residues resolved from the given groups.")
                return

            # Unique sequential IDs across structures. Remember which input group
            # each surviving instance came from (separately from the database
            # columns) so distinct per-group colors survive into the object.
            merged_rows = []
            member_colors = {}
            member_groups = {}
            for index, (row, contributing_groups) in enumerate(merged_pairs, 1):
                row.motif_id = f"{new_alias}_{index:03d}"
                merged_rows.append(row)
                ordered_groups = [part for part in all_parts if part in contributing_groups]
                member_groups[row.motif_id] = ordered_groups
                if ordered_groups:
                    member_colors[row.motif_id] = ordered_groups[0]

            # A clean Motif: line from each input group's own family name.
            families = []
            for part in all_parts:
                fam = gui.query_groups.get(part, {}).get("motif")
                if fam and fam not in families:
                    families.append(fam)
            motif_display = ", ".join(families) if families else "combined"

            gui.query_groups[new_alias] = {
                "rows": merged_rows,
                "motif_ids": [row.motif_id for row in merged_rows],
                "query": f"combine({', '.join(all_parts)})",
                "motif": motif_display,
                "structures": sorted({row.structure_id for row in merged_rows}),
                "member_groups": member_groups,
                "member_colors": member_colors,
            }
            colors.get_color(new_alias)  # reserve a stable color for the group
            gui.logger.success(
                f"Combined {len(all_parts)} group(s) into '{new_alias}': "
                f"{input_total} input row(s) -> {len(merged_rows)} instance(s) after residue merge.")
            gui.logger.info(
                "  Each database's original labels are kept as separate columns; "
                "overlapping instances retain every source's annotation.")
            gui.logger.info(
                "  The original groups and loaded annotations are unchanged.")
            return

        from .database.motif_hierarchy_cache import get_hierarchy_cache, residue_key
        cache = get_hierarchy_cache()
        if cache.alias_exists(new_alias):
            gui.logger.error(f"Alias '{new_alias}' is already in use. Please choose a different alias.")
            return

        combined_residue_keys: List[str] = []
        combined_pdb_id = None
        existing_objects = set(cmd.get_object_list())

        for name in all_parts:
            entry = cache.get_alias(name)
            if entry:
                if combined_pdb_id is None:
                    combined_pdb_id = entry['pdb_id']
                elif entry['pdb_id'] != combined_pdb_id:
                    gui.logger.error(
                        f"Cannot combine '{name}' ({entry['pdb_id']}) with a set from "
                        f"{combined_pdb_id}; rmv_combine_groups only supports members from the same PDB.")
                    return
                combined_residue_keys.extend(entry['residue_keys'])
                continue

            if name in existing_objects:
                if combined_pdb_id is None:
                    combined_pdb_id = gui.loaded_pdb_id or ''
                try:
                    pairs = {(atom.chain, int(atom.resi)) for atom in cmd.get_model(name).atom}
                except Exception as exc:
                    gui.logger.error(f"Could not read residues from object '{name}': {exc}")
                    return
                if pairs:
                    combined_residue_keys.append(residue_key(list(pairs)))
                continue

            gui.command_error(f"Unknown alias or PyMOL object: '{name}'")
            return

        if not combined_residue_keys:
            gui.command_error("Nothing to combine - no residues resolved from the given names.")
            return

        cache.create_alias(
            new_alias, combined_pdb_id or '', combined_residue_keys,
            motif_filter='', db_expression=f"combine({', '.join(all_parts)})",
            labels_snapshot=[], source_command='rmv_combine_groups',
        )
        colors.get_color(new_alias)  # force-assign a stable, unique color now

        # Best-effort visual grouping of whichever member objects/groups
        # already exist in the PyMOL scene; missing ones are skipped.
        for member in all_parts:
            if member in existing_objects or member in cmd.get_names('group_objects'):
                try:
                    cmd.group(new_alias, member)
                except Exception:
                    pass

        gui.logger.success(
            f"Alias '{new_alias}' created: combined {len(all_parts)} set(s) into "
            f"{len(combined_residue_keys)} instance(s) for {combined_pdb_id}.")
        gui.logger.info(f"  rmv_super {new_alias}     Superimpose the combined set")
        gui.logger.info(f"  rmv_align {new_alias}     Sequence-dependent superimposition")

    def load_user_annotations(tool='', pdb_id=''):
        """
        PyMOL command: Load motifs from user-uploaded annotation files.
        
        Supports: FR3D, RNAMotifScan
        
        Usage:
            rmv_user fr3d 1S72          Load FR3D annotations for 1S72
            rmv_user rnamotifscan 1A00  Load RNAMotifScan annotations
            rmv_user list               Show available user annotation files
        """
        # Handle PyMOL argument parsing - may get as single string or separate args
        tool_arg = str(tool).strip() if tool else ''
        pdb_arg = str(pdb_id).strip() if pdb_id else ''
        
        # If tool contains both tool name and pdb_id (space-separated)
        if tool_arg and not pdb_arg:
            parts = tool_arg.split()
            if len(parts) >= 2:
                tool_arg = parts[0]
                pdb_arg = parts[1]
        
        if not tool_arg:
            print("\n" + "="*60)
            print("User Annotation Loader")
            print("="*60)
            print("\nUsage: rmv_user <TOOL> <PDB_ID>")
            print("\nSupported tools:")
            print("  fr3d            FR3D output format")
            print("  rnamotifscan    RNAMotifScan output format")
            print("  rnamotifscanx   RNAMotifScanX output format")
            print("\nExamples:")
            print("  rmv_user fr3d 1S72")
            print("  rmv_user rnamotifscan 1A00")
            print("  rmv_user rnamotifscanx 1A00")
            print("  rmv_user list               Show available files")
            print("\nFile locations:")
            print("  FR3D files:        database/user_annotations/fr3d/")
            print("  RNAMotifScan:      database/user_annotations/rnamotifscan/")
            print("  RNAMotifScanX:     database/user_annotations/RNAMotifScanX/")
            print("\nFR3D (Source 5) commands (wraps external official BGSU fr3d-python):")
            print("  rmv_setup FR3D                  One-shot: install deps + register FR3D")
            print("  rmv_fr3d register <config>     Register external FR3D from a custom config")
            print("  rmv_fr3d status                 Show FR3D registration status")
            print("  rmv_load_motif                  Run FR3D search on loaded PDB")
            print("\nRNAMotifScanX wrapper commands:")
            print("  rmv_db 7                     Activate integrated Source-7 runtime")
            print("  rmv_rmsx_doctor             Validate Source-7 runtime installation")
            print("  rmv_rmsx setup              Attempt first-run runtime setup")
            print("  rmv_rmsx test")
            print("  rmv_rmsx run 1S72")
            print("  rmv_rmsx run 1S72")
            print("="*60 + "\n")
            return
        
        tool_arg = tool_arg.lower().strip()
        
        if tool_arg == 'list':
            gui._list_user_annotations()
            return
        
        if not pdb_arg:
            gui.command_error("Please specify PDB ID")
            print(f"  Usage: rmv_user {tool_arg} <PDB_ID>")
            return
        
        gui.load_user_annotations_action(tool_arg, pdb_arg)

    def fr3d_wrapper(action='', arg1='', *extra_args, **_kwargs):
        """PyMOL command: Manage the Source-5 FR3D integration.

        Source 5 wraps an external, user-installed official BGSU fr3d-python.

        Usage:
            rmv_fr3d status               Show FR3D registration / environment status
            rmv_fr3d setup [PYTHON]       Install FR3D deps (alias of rmv_setup FR3D)
            rmv_fr3d register <config>    Register FR3D from a JSON config file
            rmv_fr3d run [PDB_ID]         Run FR3D search on the loaded (or given) PDB

        Registration is normally done via `rmv_db FR3D` (bundled config) or
        `rmv_fr3d register <config>`.
        Running is normally done via a bare `rmv_load_motif`.
        """
        action_arg = str(action).strip() if action else ''
        arg1_str = str(arg1).strip() if arg1 else ''

        # Handle combined-string invocation from PyMOL
        if action_arg and not arg1_str:
            parts = action_arg.split()
            if len(parts) > 1:
                action_arg = parts[0]
                arg1_str = ' '.join(parts[1:])

        sub = action_arg.lower() if action_arg else 'status'

        if sub in ['', 'status', 'show']:
            gui.print_fr3d_status()
            return

        if sub == 'setup':
            gui.logger.info("Tip: 'rmv_setup FR3D' is the one-shot command that does this.")
            gui.auto_setup_fr3d(arg1_str)
            return

        if sub == 'register':
            if not arg1_str:
                gui.command_error("Usage: rmv_fr3d register /path/to/config.json")
                return
            gui.register_fr3d_source(arg1_str)
            return

        if sub == 'run':
            target_pdb = arg1_str or gui.loaded_pdb_id
            gui.run_fr3d_search(target_pdb)
            return

        gui.command_error(f"Unknown rmv_fr3d subcommand: {sub}")
        gui.logger.info("Use: rmv_fr3d status | setup [PYTHON] | register <config.json> | run [PDB_ID]")

    def setup_wrapper(source='', arg1='', *extra_args, **_kwargs):
        """PyMOL command: One-shot setup for an external source.

        Installs everything required and registers the source so the user only
        needs to paste the external software.

        Usage:
            rmv_setup FR3D [PYTHON]     Install deps + register FR3D (Source 5)
        """
        src = str(source or '').strip()
        arg1_str = str(arg1 or '').strip()
        if src and not arg1_str:
            parts = src.split()
            if len(parts) > 1:
                src = parts[0]
                arg1_str = ' '.join(parts[1:])
        key = src.lower()
        if key in ('', 'help'):
            gui.logger.info("Usage: rmv_setup FR3D [/abs/path/to/python]")
            gui.logger.info("  Prepares and registers an external source in one step.")
            return
        if key in ('fr3d', '5', 'source5'):
            gui.auto_setup_fr3d(arg1_str)
            return
        gui.command_error(f"rmv_setup: unknown source '{src}'")
        gui.logger.info("Supported: rmv_setup FR3D")

    def rmsx_wrapper(action='', arg1='', *extra_args, **_kwargs):
        """PyMOL command: Configure and run external RNAMotifScanX through RSMViewer.

        Usage:
            rmv_rmsx status
            rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH]
            rmv_rmsx args <ARG_TEMPLATE>
            rmv_rmsx doctor
            rmv_rmsx setup
            rmv_rmsx test
            rmv_rmsx run <PDB_ID> [EXTRA_ARGS]
            rmv_rmsx run_current [EXTRA_ARGS]
            rmv_rmsx scan_prepared <PDB_ID> [CHAINS] [compare]
            rmv_rmsx scan_cancel
        """
        action_arg = str(action).strip() if action else ''
        arg1_str = str(arg1).strip() if arg1 else ''

        if action_arg and not arg1_str:
            parts = action_arg.split()
            if len(parts) > 1:
                action_arg = parts[0]
                arg1_str = parts[1]

        sub = action_arg.lower() if action_arg else 'status'

        if sub in ['', 'status', 'show']:
            gui.print_rmsx_wrapper_status()
            return

        if sub == 'config':
            if not arg1_str:
                gui.command_error("Usage: rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH] [QUERY_FILE]")
                return
            extras = [str(x).strip() for x in extra_args if str(x).strip()]
            query_file = ''
            for idx, value in enumerate(list(extras)):
                lowered = value.lower()
                if lowered.endswith('.struct') or lowered.endswith('.txt') or lowered.startswith('query='):
                    query_file = value.split('=', 1)[1].strip() if '=' in value else value
                    extras.pop(idx)
                    break
            output_dir = extras[0] if len(extras) >= 1 else ''
            work_dir = extras[1] if len(extras) >= 2 else ''
            auto_on_fetch = extras[2] if len(extras) >= 3 else ''
            if not query_file and len(extras) >= 4:
                query_file = extras[3]
            gui.configure_rmsx_wrapper(arg1_str, output_dir, work_dir, auto_on_fetch, query_file)
            return

        if sub == 'args':
            parts = [arg1_str] if arg1_str else []
            parts.extend(str(x).strip() for x in extra_args if str(x).strip())
            template = ' '.join(parts).strip()
            gui.set_rmsx_args_template(template)
            return

        if sub == 'doctor':
            gui.rmsx_doctor(auto_setup=False)
            return

        if sub == 'setup':
            gui.rmsx_doctor(auto_setup=True)
            return

        if sub == 'test':
            gui.test_rmsx_wrapper()
            return

        if sub == 'run':
            if not arg1_str:
                gui.command_error("Usage: rmv_rmsx run <PDB_ID> [EXTRA_ARGS]")
                return
            extras = ' '.join(str(x).strip() for x in extra_args if str(x).strip())
            gui.run_rmsx_wrapper(arg1_str, extras, force_fresh=True)
            return

        if sub in ['scan_prepared', 'scan']:
            if not arg1_str:
                gui.command_error("Usage: rmv_rmsx scan_prepared <PDB_ID> [CHAINS] [compare]")
                return
            tokens = [str(x).strip() for x in extra_args if str(x).strip()]
            compare = False
            filtered = []
            for token in tokens:
                if token.lower() in ('compare', '--compare', 'validate'):
                    compare = True
                else:
                    filtered.append(token)
            gui.run_rmsx_scan_prepared(arg1_str, chains=' '.join(filtered).strip(), compare=compare)
            return

        if sub in ['scan_cancel', 'cancel']:
            gui.cancel_rmsx_scan()
            return

        if sub in ['run_current', 'current']:
            if not gui.loaded_pdb_id:
                gui.logger.error("No active structure. Use rmv_fetch <PDB_ID> first.")
                return
            extras = [arg1_str] if arg1_str else []
            extras.extend(str(x).strip() for x in extra_args if str(x).strip())
            gui.run_rmsx_wrapper(gui.loaded_pdb_id, ' '.join(extras).strip(), force_fresh=False)
            return

        gui.command_error(f"Unknown rmv_rmsx subcommand: {sub}")
        gui.logger.info("Use: rmv_rmsx status | config | args | doctor | setup | test | run | run_current | scan_prepared | scan_cancel")

    def rmsx_doctor_cmd(*_args, **_kwargs):
        """PyMOL command: Show integrated RMSX runtime diagnostics."""
        gui.rmsx_doctor(auto_setup=False)
    
    # Add commands to PyMOL
    cmd.extend('rmv_fetch', fetch_raw_pdb)
    cmd.extend('rmv_load', load_structure)
    cmd.extend('rmv_toggle', toggle_motif)
    cmd.extend('rmv_debug', set_debug)
    cmd.extend('rmv_hide', hide_annotation)
    cmd.extend('rmv_help', show_help)
    cmd.extend('rmv_bg_color', set_bg_color)
    cmd.extend('rmv_select', select_annotation)
    cmd.extend('rmv_list', list_annotations)
    cmd.extend('rmv_create_object', create_annotation_object)
    cmd.extend('rmv_db', select_database)
    cmd.extend('rmv_source', set_source)
    cmd.extend('rmv_refresh', refresh_motifs)
    cmd.extend('rmv_view', view_motif)
    cmd.extend('rmv_combine_groups', combine_sets)
    cmd.extend('rmv_fr3d', fr3d_wrapper)
    cmd.extend('rmv_setup', setup_wrapper)
    cmd.extend('rmv_rmsx', rmsx_wrapper)
    cmd.extend('rmv_rmsx_doctor', rmsx_doctor_cmd)
    
    def show_colors():
        """PyMOL command: Show color legend for all motif types."""
        from . import colors as color_module
        loaded = gui.viz_manager.motif_loader.get_loaded_motifs()
        if loaded:
            color_module.print_color_legend(loaded)
        else:
            color_module.print_color_legend()
    
    cmd.extend('rmv_colors', show_colors)
    
    def set_motif_color(motif_type='', color=''):
        """PyMOL command: Change color of a specific motif type.
        
        Usage:
            rmv_color HL red         Change HL to red
            rmv_color GNRA blue      Change GNRA to blue
            rmv_color IL 0.5 1.0 0.5 Change IL to RGB values
        
        Available colors: red, green, blue, yellow, cyan, magenta, orange,
                         pink, purple, teal, gold, coral, turquoise, etc.
        """
        if not motif_type:
            print("\nUsage: rmv_color <MOTIF_TYPE> <COLOR>")
            print("Examples:")
            print("  rmv_color HL red")
            print("  rmv_color GNRA blue")
            print("  rmv_color IL green")
            print("\nAvailable colors: red, green, blue, yellow, cyan, magenta,")
            print("                  orange, pink, purple, teal, gold, coral, etc.")
            return
        
        if not color:
            gui.command_error("Please specify a color")
            gui.logger.error("Example: rmv_color HL red")
            return
        
        from . import colors as color_module
        
        motif_arg = str(motif_type).strip().upper()
        color_arg = str(color).strip().lower()
        
        # Set the custom color
        result = color_module.set_custom_motif_color(motif_arg, color_arg)
        
        gui.logger.success(f"Changed {motif_arg} color to {color_arg}")
        
        # Re-apply color to currently loaded motifs if any
        loaded_motifs = gui.viz_manager.motif_loader.get_loaded_motifs()
        if motif_arg in loaded_motifs:
            info = loaded_motifs[motif_arg]
            structure_name = info.get('structure_name')
            main_selection = info.get('main_selection')
            
            # Re-color the motif residues in the structure
            if main_selection:
                try:
                    color_module.set_motif_color_in_pymol(cmd, main_selection, motif_arg)
                    gui.logger.info(f"Applied new color to {motif_arg} residues")
                except Exception as e:
                    gui.logger.debug(f"Could not apply color: {e}")
        
        print(f"\n  {motif_arg} is now colored {color_arg}")
        print(f"  Use 'rmv_show {motif_arg}' or 'rmv_show ALL' to see the change\n")
    
    cmd.extend('rmv_color', set_motif_color)

    def set_property(subject='', *extra_args, **_kwargs):
        """PyMOL command: Set a property on an alias/group or motif type.

        Currently supports 'color' (Application 6: distinct novel/known colors).

        Usage:
            rmv_set color, novel_SR, red
            rmv_set color, known_SR, blue
        """
        first_parts = str(subject).split() if subject else []
        rest_parts = [str(a) for a in extra_args]
        all_parts = [p.strip() for p in (first_parts + rest_parts) if p.strip()]

        if len(all_parts) < 3 or all_parts[0].lower() != 'color':
            gui.command_error("Usage: rmv_set color, <ALIAS_OR_TYPE>, <COLOR>")
            gui.logger.error("Example: rmv_set color, novel_SR, red")
            return

        target = all_parts[1]
        color_arg = all_parts[2].lower()

        from . import colors as color_module
        # Cache the preference under the target's name so future object
        # creation (rmv_show <alias>, rmv_super <alias>) also uses it.
        color_module.set_custom_motif_color(target, color_arg)

        applied = False
        # Stable-ID query group (from rmv_select): recolor any objects that
        # were already created for its members, and remember the preference
        # for future rmv_create_object/rmv_view calls.
        if target in gui.query_groups:
            member_objects = [
                f"motif_{mid}" for mid in gui.query_groups[target].get("motif_ids", [])
            ]
            live = set(cmd.get_object_list())
            for obj in member_objects:
                if obj in live:
                    try:
                        color_module.set_motif_color_in_pymol(cmd, obj, target)
                        applied = True
                    except Exception as exc:
                        gui.logger.error(f"Could not apply color to '{obj}': {exc}")
            if applied:
                gui.logger.success(f"Set color of group '{target}' to {color_arg}")
            else:
                gui.logger.success(
                    f"Color preference for group '{target}' set to {color_arg}; "
                    f"it will apply on rmv_create_object/rmv_view.")
            return

        existing_names = set(cmd.get_object_list()) | set(cmd.get_names('group_objects'))
        if target in existing_names:
            # Alias/group or plain PyMOL object - recolor directly.
            try:
                cmd.color(color_arg, target)
                applied = True
            except Exception as exc:
                gui.logger.error(f"Could not apply color to '{target}': {exc}")
        else:
            # Fall back to the regular loaded-motif-type recoloring path.
            loaded_motifs = gui.viz_manager.motif_loader.get_loaded_motifs()
            motif_key = target.upper()
            if motif_key in loaded_motifs:
                main_selection = loaded_motifs[motif_key].get('main_selection')
                if main_selection:
                    try:
                        color_module.set_motif_color_in_pymol(cmd, main_selection, motif_key)
                        applied = True
                    except Exception as exc:
                        gui.logger.error(f"Could not apply color to '{motif_key}': {exc}")

        if applied:
            gui.logger.success(f"Set color of '{target}' to {color_arg}")
        else:
            gui.logger.warning(
                f"Color preference for '{target}' saved, but no existing "
                f"PyMOL object/group/motif type named '{target}' was found to recolor now."
            )

    cmd.extend('rmv_set', set_property)

    def set_color_alias(*args, **_kwargs):
        """PyMOL command alias: rmv_set_color <GROUP_OR_TYPE>, <COLOR>.

        Convenience wrapper around 'rmv_set color, ...' so the documented
        Application 6 syntax works verbatim.
        """
        parts = [str(a).strip() for a in args if str(a).strip()]
        if parts and parts[0].lower() == 'color':
            parts = parts[1:]
        if len(parts) < 2:
            gui.command_error("Usage: rmv_set_color <GROUP_OR_TYPE>, <COLOR>")
            gui.logger.error("Example: rmv_set_color group_FP, red")
            return
        set_property('color', *parts)

    cmd.extend('rmv_set_color', set_color_alias)

    def save_motif_images(argument=''):
        """PyMOL command: Save a group/motif-ID as mmCIF, or save an image.

        Create a group first (rmv_select ... as <group>) or use a stable motif
        ID (e.g. 1S72_00016). Every save prints its output location.

        Usage (mmCIF structure export, original coordinates):
            rmv_save <group> cif              Export a saved group as mmCIF
            rmv_save <MOTIF_ID> cif           Export one stable motif ID as mmCIF
            rmv_save ALL cif                  Export every consolidated motif as mmCIF

        Usage (image / view save, PNG):
            rmv_save current                  Save the current PyMOL view
            rmv_save current my_view.png      Save the current view to a named file
            rmv_save ALL                      Save an image of every motif (cartoon)
            rmv_save ALL sticks               Save all motif images in a representation

        mmCIF export extracts ORIGINAL coordinates from the on-disk CIF file,
        NOT PyMOL's internal coordinates (which may be slightly modified).
        Output is a minimal coordinates-only mmCIF containing filtered
        _atom_site rows for motif residues.

        Available representations (for image save):
            - cartoon       (default) - Shows RNA backbone ribbon
            - sticks        - Shows all atoms as sticks
            - spheres       - Shows all atoms as spheres
            - ribbon        - Simplified backbone ribbon
            - lines         - Wire representation
            - licorice      - Thick bonds representation
            - surface       - Molecular surface
            - cartoon+sticks - Combination of cartoon and sticks

        Output folder structure:
            Images:     plugin_dir/motif_images/pdb_id/
            Structures: plugin_dir/motif_structures/pdb_id/
        """
        arguments = str(argument).strip().split()
        
        if not arguments:
            print("\nUsage: rmv_save <group | MOTIF_ID | ALL | current> [cif | representation]")
            print("\n  First create a group (rmv_select ... as <group>) or use a stable")
            print("  motif ID (e.g. 1S72_00016). Every save prints its output location.")
            print("\n  STRUCTURE EXPORT (mmCIF, original coordinates):")
            print("    rmv_save <group> cif         Export a saved group as mmCIF")
            print("    rmv_save <MOTIF_ID> cif      Export one stable motif ID (e.g. 1S72_00016)")
            print("    rmv_save ALL cif             Export every consolidated motif as mmCIF")
            print("\n  IMAGE / VIEW SAVE (PNG):")
            print("    rmv_save current             Save the current PyMOL view")
            print("    rmv_save current out.png     Save the current view to a named file")
            print("    rmv_save ALL                 Save an image of every motif (cartoon)")
            print("    rmv_save ALL sticks          Save all motif images in a representation")
            print("\n  Representations: cartoon (default), sticks, spheres, ribbon, lines, licorice, surface, cartoon+sticks")
            print("\n  Note: mmCIF export uses ORIGINAL coordinates from the on-disk CIF,")
            print("        not PyMOL's internal coordinates (coordinates-only _atom_site rows).")
            print("\n  Output locations (printed after each save):")
            print("    Structures: plugin_dir/motif_structures/<pdb_id>/")
            print("    Images:     plugin_dir/motif_images/<pdb_id>/")
            print("    View PNG:   the filename you provide (or the working directory)")
            print("\n  Example:")
            print("    rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas or RNAMotifScanX, as group_SR")
            print("    rmv_save group_SR cif")
            print("    rmv_save current group_SR.png")
            return
        
        pdb_id = gui.viz_manager.structure_loader.get_current_pdb_id()
        if not pdb_id:
            gui.logger.error("No structure loaded")
            return

        if arguments and len(arguments) >= 2 and arguments[-1].lower() in ('cif', 'mmcif'):
            stable_target = arguments[0]
            if stable_target in gui.query_groups or any(
                stable_target in table._rows for table in gui.annotation_tables.values()
            ):
                gui.export_annotation_rows(stable_target)
                return
        
        loaded_motifs = gui.viz_manager.motif_loader.get_loaded_motifs()
        if not loaded_motifs:
            gui.logger.error("No motifs loaded for this structure")
            return
        
        # Helper: check if a string is the 'cif' keyword
        def _is_cif(s):
            return s.lower() in ('cif', 'mmcif')
        
        # Helper: suggest closest match for a mistyped argument
        def _suggest(word, candidates):
            """Return the closest candidate using simple edit-distance heuristic."""
            word_up = word.upper()
            # Exact prefix match first
            prefix_hits = [c for c in candidates if c.startswith(word_up)]
            if prefix_hits:
                return prefix_hits[0]
            # Substring match
            sub_hits = [c for c in candidates if word_up in c or c in word_up]
            if sub_hits:
                return sub_hits[0]
            # Simple character-overlap score
            def _score(a, b):
                a, b = a.upper(), b.upper()
                if not a or not b:
                    return 0
                common = sum(1 for c in a if c in b)
                return common / max(len(a), len(b))
            best = max(candidates, key=lambda c: _score(word, c), default=None)
            if best and _score(word, best) > 0.4:
                return best
            return None
        
        representation = 'cartoon'  # Default
        
        if arguments[0].upper() == 'ALL':
            # rmv_save ALL [representation | cif]
            if len(arguments) > 1 and _is_cif(arguments[1]):
                gui.export_all_motif_structures_action()
            else:
                if len(arguments) > 1:
                    representation = arguments[1].lower()
                gui.save_all_motif_images_action(representation=representation)
        
        elif arguments[0].upper() == 'CURRENT':
            # Save current view: rmv_save current [filename]
            if len(arguments) > 1:
                filename = arguments[1]
            else:
                from datetime import datetime
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                filename = f"pymol_view_{timestamp}.png"
            gui.save_current_view_action(filename)
        
        else:
            # rmv_save <TYPE> [INSTANCE_ID] [representation | cif]
            motif_type = arguments[0].upper()
            
            if motif_type not in loaded_motifs:
                # A saved group / stable motif ID needs 'cif' to export, or use
                # 'rmv_save current <file>.png' to save an image of the view.
                raw = arguments[0]
                if raw in gui.query_groups or any(
                    raw in table._rows for table in gui.annotation_tables.values()
                ):
                    gui.logger.error(
                        f"'{raw}' is a saved group or motif ID; add 'cif' to export it, "
                        f"or save an image of the current view."
                    )
                    gui.logger.info(f"  rmv_save {raw} cif            Export as mmCIF")
                    gui.logger.info(f"  rmv_save current {raw}.png    Save the current view as PNG")
                    return
                # Check for possible typos against known keywords and motif types
                all_candidates = ['ALL', 'CURRENT'] + sorted(loaded_motifs.keys())
                suggestion = _suggest(arguments[0], all_candidates)
                gui.command_error(f"Unknown argument '{arguments[0]}'")
                if suggestion:
                    gui.logger.info(f"Did you mean: rmv_save {suggestion}?")
                gui.logger.info(f"Available: {', '.join(sorted(loaded_motifs.keys()))}")
                gui.logger.info("Other options: ALL, CURRENT")
                return
            
            if len(arguments) == 1:
                # rmv_save HL  ->  save all HL images (default cartoon)
                gui.save_motif_type_images_action(motif_type, representation=representation)
            
            elif len(arguments) == 2:
                arg2 = arguments[1]
                if _is_cif(arg2):
                    # rmv_save HL cif  ->  export all HL structures
                    gui.export_motif_type_structures_action(motif_type)
                else:
                    try:
                        instance_id = int(arg2)
                        # rmv_save HL 3  ->  save HL instance #3 image
                        gui.save_motif_instance_by_id_action(motif_type, instance_id,
                                                            representation=representation)
                    except ValueError:
                        # rmv_save HL sticks  ->  save all HL images as sticks
                        representation = arg2.lower()
                        gui.save_motif_type_images_action(motif_type, representation=representation)
            
            elif len(arguments) >= 3:
                arg2 = arguments[1]
                arg3 = arguments[2]
                try:
                    instance_id = int(arg2)
                    if _is_cif(arg3):
                        # rmv_save HL 3 cif  ->  export HL instance #3 as mmCIF
                        gui.export_motif_instance_by_id_action(motif_type, instance_id)
                    else:
                        # rmv_save HL 3 spheres  ->  save HL instance #3 as spheres
                        representation = arg3.lower()
                        gui.save_motif_instance_by_id_action(motif_type, instance_id,
                                                            representation=representation)
                except ValueError:
                    # arg2 is not an integer - treat as representation
                    representation = arg2.lower()
                    gui.save_motif_type_images_action(motif_type, representation=representation)
    
    cmd.extend('rmv_save', save_motif_images)
    
    def show_chain_diagnostics(structure_name=''):
        """PyMOL command: Show chain ID diagnostic information for a loaded structure.
        
        Usage:
            rmv_chains              Show chains for current structure
            rmv_chains 1s72         Show chains for specific structure
        """
        try:
            # Determine structure name
            if not structure_name:
                structure_name = gui.loaded_pdb if hasattr(gui, 'loaded_pdb') and gui.loaded_pdb else ''
            
            if not structure_name:
                print("\n  No structure specified. Usage: rmv_chains <structure_name>")
                return
            
            structure_name = structure_name.strip().lower()
            
            # Read current cif_use_auth from GUI state
            cif_auth_val = getattr(gui, 'cif_use_auth', 1)
            chain_mode = "auth_asym_id" if cif_auth_val == 1 else "label_asym_id"
            chain_label = "Auth chains" if cif_auth_val == 1 else "Label chains"
            
            # Get chains
            try:
                chains = cmd.get_chains(structure_name)
            except Exception as e:
                print(f"\n  ERROR: Could not get chains for '{structure_name}': {e}")
                return
            
            # Format chains in rows of 20
            print(f"\n  Structure: {structure_name.upper()}  |  cif_use_auth = {cif_auth_val} ({chain_mode})  |  Chains: {len(chains)}")
            print(f"  {chain_label}: ", end="")
            for i, ch in enumerate(chains):
                if i > 0 and i % 20 == 0:
                    print(f"\n               ", end="")
                print(f" {ch}", end="")
            print("\n")
            
        except Exception as e:
            print(f"\n  Error in chain diagnostics: {e}\n")
    
    cmd.extend('rmv_chains', show_chain_diagnostics)

    def show_loaded_tags(*args, **kwargs):
        """PyMOL command: Show all loaded PDB+source combination tags.

        Usage:
            rmv_loaded             Show currently loaded PDB_SRC tags
        """
        # Derive tags from both loaded_sources and loaded_motifs metadata
        tags_set = set()
        for pdb, suffix in gui.loaded_sources:
            tags_set.add(f"{pdb}{suffix}")
        loaded_motifs = gui.viz_manager.motif_loader.get_loaded_motifs()
        if loaded_motifs:
            for info in loaded_motifs.values():
                default_pdb = info.get('pdb_id', '').upper()
                default_sfx = info.get('source_suffix', '')
                for detail in info.get('motif_details', []):
                    d_pdb = detail.get('_pdb_id', default_pdb).upper()
                    d_sfx = detail.get('_source_suffix', default_sfx)
                    if d_pdb and d_sfx:
                        tags_set.add(f"{d_pdb}{d_sfx}")
        if not tags_set:
            print("\n  No PDB+source combinations loaded yet.")
            print("  Load data first:")
            print("    rmv_fetch 1S72")
            print("    rmv_db 7")
            print("    rmv_load_motif\n")
            return

        tags = sorted(tags_set)
        print(f"\n  Loaded PDB+source tags ({len(tags)}):")
        for t in tags:
            print(f"    {t}")
        print(f"\n  Use these tags with rmv_super / rmv_align:")
        print(f"    rmv_super MOTIF_TYPE, {', '.join(tags[:2])}")
        print()

    cmd.extend('rmv_loaded', show_loaded_tags)

    def reset_plugin():
        """PyMOL command: Reset everything - delete all objects and reset plugin to defaults.
        
        Usage:
            rmv_reset              Delete all PyMOL objects, reset plugin state
        """
        # Step 1: Delete all PyMOL objects
        try:
            cmd.delete('all')
            gui.logger.debug("Deleted all PyMOL objects")
        except Exception as e:
            gui.logger.debug(f"Could not delete objects: {e}")
        
        # Step 2: Reset all plugin state to defaults
        gui.loaded_pdb = None
        gui.loaded_pdb_id = None
        gui.loaded_structures = {}
        gui.annotation_tables = {}
        gui.query_groups = {}
        gui.motif_visibility = {}
        gui.current_source_mode = None
        gui.current_source_names = []
        gui.current_user_tool = None
        gui.current_local_source = None
        gui.current_web_source = None
        gui.combined_source_ids = []
        gui.current_source_id = None
        gui.user_rms_filtering_enabled = True
        gui.user_rmsx_filtering_enabled = True
        gui.user_rms_custom_pvalues = {}
        gui.user_rmsx_custom_pvalues = {}
        gui.cif_use_auth = 1
        gui.auth_to_label_map = {}
        gui.loaded_sources = set()
        gui.pdb_source_state = {}

        try:
            from .database.motif_hierarchy_cache import get_hierarchy_cache, close_hierarchy_cache
            cache = get_hierarchy_cache()
            cache.clear_all_hierarchy_data()
            cache.reset_all_aliases()
            cache_path = Path(cache.db_path)
            close_hierarchy_cache()
            for suffix in ("", "-wal", "-shm"):
                candidate = Path(str(cache_path) + suffix)
                if candidate.exists():
                    candidate.unlink()
        except Exception:
            pass

        # Step 2b: Clear the on-disk API response cache (~/.rsmviewer_cache/)
        # so the next rmv_load_motif hits the live API instead of returning
        # a file cached up to 30 days ago.
        try:
            from .database.cache_manager import get_cache_manager
            removed = get_cache_manager().clear_cache()
            if removed:
                gui.logger.debug(f"Cleared {removed} cached API response(s) from disk")
        except Exception:
            pass

        # Step 2c: Clear each provider's own in-process memory cache. Provider
        # instances are process-lifetime singletons (held by the source
        # selector), so a plain dict lookup like `_motif_cache[pdb_id]` would
        # keep returning the first-ever fetch for a PDB even after the disk
        # and SQLite caches above are wiped.
        try:
            from .database import get_source_selector
            source_selector = get_source_selector()
            if source_selector:
                for provider in source_selector.providers.values():
                    for attr in ('_motif_cache', '_pdb_motif_cache', '_annotation_cache', '_fetched_pdbs'):
                        cache_obj = getattr(provider, attr, None)
                        if cache_obj is not None:
                            cache_obj.clear()
        except Exception:
            pass

        # Step 2d: Clear the on-disk RMSX preannotated extraction cache so the
        # next load re-extracts fresh consensus logs from the archive/folder.
        try:
            import shutil
            rmsx_out = Path(getattr(gui, 'rmsx_output_path', '') or '')
            if rmsx_out:
                preannotated_cache = rmsx_out / '.preannotated_cache'
                if preannotated_cache.exists():
                    shutil.rmtree(preannotated_cache, ignore_errors=True)
                    gui.logger.debug(f"Cleared preannotated RMSX cache: {preannotated_cache}")
        except Exception:
            pass
        
        # Step 3: Reset chain ID convention to default
        try:
            cmd.set("cif_use_auth", 1)
        except:
            pass
        
        # Step 4: Clear motif loader data
        try:
            if gui.viz_manager and gui.viz_manager.motif_loader:
                gui.viz_manager.motif_loader.loaded_motifs = {}
        except:
            pass
        
        # Step 5: Reset colors
        try:
            from . import colors as color_module
            color_module.CUSTOM_COLORS.clear()
            color_module._dynamic_assigned.clear()
            color_module._dynamic_color_index = 0
        except:
            pass
        
        gui.logger.success("Plugin reset to defaults")
        print("\n  All objects deleted and plugin state and caches cleared.")
        print("  Ready for a fresh session.")
        print("\n  Quick Start:")
        print("     rmv_fetch 1S72                          # Load a PDB structure")
        print("     rmv_db RNA3DMotifAtlas,RNAMotifScanX    # Select and load sources")
        print("     rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR")
        print("     rmv_list group_SR                       # Inspect the motif IDs")
        print("     rmv_view group_SR                       # Highlight the family")
        print("     rmv_create_object group_SR              # Create selectable objects")
        print("     rmv_super group_SR                      # Medoid superimposition")
        print()
    
    cmd.extend('rmv_reset', reset_plugin)

    # Optional command modules. A failure to register any of these must NOT
    # abort plugin initialization or hamper the core pipeline (rmv_fetch/
    # rmv_db/rmv_select/rmv_view/rmv_create_object are already registered
    # above). Each is isolated so one bad import can't take down the rest.
    try:
        from .pair_visualizer import register_pair_commands
        register_pair_commands()
    except Exception as exc:
        gui.logger.warning(f"rmv_pair/rmv_pair_batch unavailable: {type(exc).__name__}: {exc}")

    try:
        from .alignment import register_alignment_commands
        register_alignment_commands()
    except Exception as exc:
        gui.logger.warning(f"rmv_super/rmv_align unavailable: {type(exc).__name__}: {exc}")

    gui.logger.success("RSMViewer GUI initialized")
