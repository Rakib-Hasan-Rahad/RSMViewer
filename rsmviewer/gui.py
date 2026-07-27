"""
RSMViewer - GUI Module
Provides PyMOL GUI interface for the plugin with multi-database support.

This module provides:
- MotifVisualizerGUI: Main GUI class for the plugin
- PyMOL command registration (24 commands)
- Database selection and switching functionality
- Multi-source combine mode with cascade merge

Author: CBB LAB @Rakib Hasan Rahad
Version: 1.0.0
Last Updated: 02 April 2026
"""

from typing import List, Optional, Dict, Tuple
import json
import os
import re
import shlex
import subprocess
import sys
from pymol import cmd
from .loader import VisualizationManager
from .utils import get_logger
from . import colors
from .database import get_registry
from .database.config import SOURCE_ID_MAP
from pathlib import Path
import platform


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
            normalized.setdefault(category, []).append(inst)
    return normalized


class MotifVisualizerGUI:
    """PyMOL GUI for RNA motif visualization with multi-database support."""

    def __init__(self):
        """Initialize GUI components."""
        self.logger = get_logger()

        # Get path to motif database
        plugin_dir = Path(__file__).parent
        self.database_dir = plugin_dir / 'motif_database'

        # Initialize visualization manager
        self.viz_manager = VisualizationManager(cmd, str(self.database_dir))

        # Track UI state
        self.motif_visibility = {}

        # Track currently loaded PDB
        self.loaded_pdb = None
        self.loaded_pdb_id = None

        # Track current source mode
        self.current_source_mode = None
        self.current_user_tool = None
        self.current_local_source = None      # 'atlas', 'rfam', or None (for both)
        self.current_web_source = None         # 'bgsu', 'rfam', or None (for auto)
        self.combined_source_ids = []          # List of source IDs for combining
        self.current_source_id = None          # Numeric source ID (1-8) for object tagging

        # Track filtering state for RMS, RMSX, and NoBIAS (for user annotations)
        self.user_rms_filtering_enabled = True   # Default: filters ON
        self.user_rmsx_filtering_enabled = True  # Default: filters ON
        self.user_nobias_filtering_enabled = True  # Default: filters ON

        # Track custom P-values for RMS, RMSX, and NoBIAS
        self.user_rms_custom_pvalues = {}        # Dict: {motif_name: p_value}
        self.user_rmsx_custom_pvalues = {}       # Dict: {motif_name: p_value}
        self.user_nobias_custom_pvalues = {}     # Dict: {motif_name: p_value}

        # Track chain ID convention: 1 = auth_asym_id (default), 0 = label_asym_id
        self.cif_use_auth = 1
        self.auth_to_label_map = {}   # auth_asym_id -> label_asym_id chain mapping

        # Track custom user annotation data paths per source (for rmv_db 5/6/7/8 with path)
        self.user_data_paths = {}  # Dict: {source_id_int: path_str}

        # Jaccard similarity threshold for cascade merge (0.0–1.0)
        self.jaccard_threshold = 0.60

        # Store within-source dedup stats for display in source attribution report
        # Dict: {source_id: (before_count, after_count)}
        self.dedup_stats = {}

        # Track loaded PDB+source combinations for cross-PDB superimposition
        # Each entry is a tuple (pdb_id_upper, source_suffix) e.g. ("1S72", "_S7")
        self.loaded_sources = set()

        # FR3D wrapper runtime settings (persisted in user_annotations/fr3d)
        fr3d_dir = plugin_dir / 'database' / 'user_annotations' / 'fr3d'
        self.fr3d_config_file = fr3d_dir / '.fr3d_wrapper.json'
        self.fr3d_root_path = ''
        self.fr3d_python_path = ''
        self.fr3d_input_path = ''
        self.fr3d_output_path = str(fr3d_dir.resolve())
        self.fr3d_auto_run_on_fetch = False
        self.fr3d_runtime_dir = str((plugin_dir / 'tools' / 'fr3d_runtime').resolve())
        self.fr3d_setup_attempted = False
        # Local FR3D pipeline config (set via rmv_db 5 /path/to/fr3d_pipeline_config.json)
        self.fr3d_pipeline_config_path = ''
        self.fr3d_pipeline_config = {}
        self._load_fr3d_wrapper_config()

        # RNAMotifScanX wrapper runtime settings (persisted in user_annotations/RNAMotifScanX)
        rmsx_dir = plugin_dir / 'database' / 'user_annotations' / 'RNAMotifScanX'
        self.rmsx_config_file = rmsx_dir / '.rmsx_wrapper.json'
        self.rmsx_executable_path = ''
        self.rmsx_output_path = str(rmsx_dir.resolve())
        self.rmsx_working_dir = ''
        self.rmsx_args_template = ''
        self.rmsx_auto_run_on_fetch = False
        # RMSX runtime config (integrated source-7 runtime)
        self.rmsx_pipeline_config_path = ''
        self.rmsx_pipeline_config = {}
        self.rmsx_runtime_dir = str((plugin_dir / 'tools' / 'rmsx_runtime').resolve())
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

        system = platform.system().lower()
        machine = platform.machine().lower()
        if system == 'darwin':
            platform_dir = 'macos-arm64' if machine in ('arm64', 'aarch64') else 'macos-x86_64'
        elif system == 'windows':
            platform_dir = 'windows-x86_64'
        else:
            platform_dir = 'linux-x86_64'

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
            runtime_dir.parent.parent / 'RNAMotifScanX' / 'PDB_prebuild.tgz',
            runtime_dir.parent / 'RNAMotifScanX' / 'PDB_prebuild.tgz',
        ]
        for candidate in archive_candidates:
            if candidate.exists() and candidate.is_file():
                prebuild_archive = str(candidate.resolve())
                break

        return {
            'rmsx_executable': rmsx_exe,
            'mc_annotate_executable': mc_exe,
            'rnaview_executable': rnaview_exe,
            'rnaview_dir': rnaview_dir,
            'pdb_prebuild_archive': prebuild_archive,
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
        """Get source suffix for PyMOL object naming (e.g., '_S3' for source 3).
        Returns empty string if no source is explicitly set."""
        if self.current_source_id is not None:
            cid = str(self.current_source_id)
            if '_' in cid:  # Combine mode (e.g., "8_7")
                return f"_S_{cid}"
            return f"_S{cid}"
        return ""

    def _load_fr3d_wrapper_config(self):
        """Load persisted FR3D wrapper settings if available."""
        try:
            cfg_path = self.fr3d_config_file
            if not cfg_path.exists():
                return

            with open(cfg_path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            self.fr3d_root_path = str(data.get('fr3d_root_path', '') or '').strip()
            self.fr3d_python_path = str(data.get('fr3d_python_path', '') or '').strip()
            self.fr3d_input_path = str(data.get('fr3d_input_path', '') or '').strip()
            output_path = str(data.get('fr3d_output_path', '') or '').strip()
            if output_path:
                self.fr3d_output_path = output_path
        except Exception as e:
            self.logger.warning(f"Could not read FR3D wrapper config from {self.fr3d_config_file}: {type(e).__name__}: {e}")

    def _save_fr3d_wrapper_config(self):
        """Persist FR3D wrapper settings."""
        try:
            self.fr3d_config_file.parent.mkdir(parents=True, exist_ok=True)
            payload = {
                'fr3d_root_path': self.fr3d_root_path,
                'fr3d_python_path': self.fr3d_python_path,
                'fr3d_input_path': self.fr3d_input_path,
                'fr3d_output_path': self.fr3d_output_path,
            }
            with open(self.fr3d_config_file, 'w', encoding='utf-8') as f:
                json.dump(payload, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Could not save FR3D wrapper config to {self.fr3d_config_file}: {type(e).__name__}: {e}")

    def print_fr3d_wrapper_status(self):
        """Print current FR3D wrapper settings and quick usage hints."""
        print("\n" + "=" * 70)
        print("FR3D Wrapper Status")
        print("=" * 70)
        print(f"Output path    : {self.fr3d_output_path}")
        print(f"Auto on fetch  : {'on' if self.fr3d_auto_run_on_fetch else 'off'}")

        # Show integrated runtime status.
        report = self._run_fr3d_runtime_setup(build=False)
        print(f"\nRuntime mode   : {report.get('runtime_mode', 'not_ready')}")
        print(f"FR3D tool path : {report.get('fr3d_python_dir', self.fr3d_root_path) or '(not found)'}")
        print(f"Python path    : {report.get('python_executable', self.fr3d_python_path) or '(auto-detect)'}")

        print("\nCommands:")
        print("  rmv_fr3d status           Show this status")
        print("  rmv_fr3d doctor           Diagnose Source-5 runtime readiness")
        print("  rmv_fr3d setup            Run one-time Source-5 runtime setup")
        print("  rmv_fr3d refresh [PDB_ID] Force fresh rerun")
        print("  rmv_fr3d run <PDB_ID>     Run local FR3D pipeline")
        print("  rmv_fr3d run_current      Run for the currently loaded PDB")
        print("  rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]  (advanced)")
        print("\nExecution policy:")
        print("  - Source 5 runs the local bundled FR3D pipeline fresh")
        print("  - No cached FR3D result reuse")
        print("  - No BGSU fallback in Source-5 flow")
        print("=" * 70 + "\n")

    def configure_fr3d_wrapper(self, fr3d_root: str, output_dir: str = '', input_dir: str = '', python_executable: str = ''):
        """Set FR3D wrapper paths using absolute paths and persist them."""
        if not fr3d_root:
            self.logger.error("FR3D root path is required")
            self.logger.info("Usage: rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [INPUT_DIR] [PYTHON]")
            return False

        expanded_root = os.path.abspath(os.path.expanduser(fr3d_root))
        if not os.path.isdir(expanded_root):
            self.logger.error(f"FR3D root directory not found: {expanded_root}")
            return False

        expected_pkg = os.path.join(expanded_root, 'fr3d', '__init__.py')
        if not os.path.isfile(expected_pkg):
            self.logger.error("Invalid FR3D root: expected fr3d package not found")
            self.logger.info("Expected file: <FR3D_ROOT>/fr3d/__init__.py")
            return False

        self.fr3d_root_path = expanded_root

        if output_dir:
            self.fr3d_output_path = os.path.abspath(os.path.expanduser(output_dir))
        if input_dir:
            self.fr3d_input_path = os.path.abspath(os.path.expanduser(input_dir))
        if python_executable:
            self.fr3d_python_path = os.path.abspath(os.path.expanduser(python_executable))
        elif not self.fr3d_python_path:
            detected_python = self._auto_detect_fr3d_python()
            if detected_python:
                self.fr3d_python_path = detected_python

        os.makedirs(self.fr3d_output_path, exist_ok=True)
        if self.fr3d_input_path and not os.path.isdir(self.fr3d_input_path):
            os.makedirs(self.fr3d_input_path, exist_ok=True)

        if self.fr3d_python_path and not os.path.isfile(self.fr3d_python_path):
            self.logger.error(f"Python executable not found: {self.fr3d_python_path}")
            return False

        self._save_fr3d_wrapper_config()

        self.logger.success("FR3D wrapper configuration updated")
        self.logger.info(f"FR3D root path: {self.fr3d_root_path}")
        self.logger.info(f"Python path: {self.fr3d_python_path or '(auto-detect)'}")
        self.logger.info(f"Output path: {self.fr3d_output_path}")
        if self.fr3d_input_path:
            self.logger.info(f"Input path: {self.fr3d_input_path}")
        else:
            self.logger.info("Input path: (empty, FR3D will auto-download when needed)")
        return True

    def _python_has_fr3d_deps(self, python_executable: str) -> bool:
        """Check whether a candidate Python executable can import FR3D dependencies."""
        if not python_executable or not os.path.isfile(python_executable):
            return False

        try:
            result = subprocess.run(
                [python_executable, '-c', 'import scipy, pdbx'],
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            return result.returncode == 0
        except Exception:
            return False

    def _auto_detect_fr3d_python(self) -> str:
        """Try to find a Python executable that can import scipy and pdbx."""
        candidates = []

        if self.fr3d_root_path:
            venv_python = Path(self.fr3d_root_path) / '.fr3d-venv' / 'bin' / 'python'
            candidates.append(str(venv_python))

        if self.fr3d_python_path:
            candidates.append(self.fr3d_python_path)
        candidates.extend([
            '/opt/homebrew/bin/python3',
            '/usr/local/bin/python3',
            sys.executable,
        ])

        seen = set()
        for candidate in candidates:
            candidate = os.path.abspath(os.path.expanduser(candidate))
            if candidate in seen:
                continue
            seen.add(candidate)
            if self._python_has_scipy(candidate):
                return candidate
        return ''
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

    def _build_internal_fr3d_config(self) -> Dict:
        """Build Source-5 FR3D runtime config from bundled runtime assets."""
        runtime_dir = Path(self.fr3d_runtime_dir)
        cache_dir = runtime_dir / 'cache'
        scripts_dir = runtime_dir / 'scripts'
        python_dir = runtime_dir / 'python'
        bin_dir = runtime_dir / 'bin'

        return {
            'runtime_dir': str(runtime_dir.resolve()),
            'fr3d_python_dir': self.fr3d_root_path,
            'python_executable': self.fr3d_python_path or sys.executable,
            'output_dir': self.fr3d_output_path,
            'cache_dir': str(cache_dir.resolve()),
            'scripts_dir': str(scripts_dir.resolve()),
            'python_dir': str(python_dir.resolve()),
            'bin_dir': str(bin_dir.resolve()),
            'cif_input_dir': self.fr3d_input_path,
            'auto_download_cif': True,
            'loop_types': ['HL', 'IL', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8'],
        }

    def _run_fr3d_runtime_setup(self, build: bool = False) -> Dict:
        """Run integrated Source-5 runtime doctor/setup script and return parsed report."""
        setup_script = Path(self.fr3d_runtime_dir) / 'setup_runtime.py'
        if not setup_script.exists():
            return {
                'ok': False,
                'runtime_mode': 'not_ready',
                'missing': ['setup_runtime.py'],
                'setup_message': f'Missing setup script: {setup_script}',
                'fr3d_python_dir': '',
                'python_executable': '',
                'output_dir': self.fr3d_output_path,
                'cache_dir': str((Path(self.fr3d_runtime_dir) / 'cache').resolve()),
            }

        cmdline = [
            sys.executable,
            str(setup_script),
            '--runtime-dir',
            self.fr3d_runtime_dir,
            '--json',
        ]
        if build:
            cmdline.append('--build')

        try:
            proc = subprocess.run(cmdline, capture_output=True, text=True, check=False)
            output = (proc.stdout or '').strip()
            if not output:
                return {
                    'ok': False,
                    'runtime_mode': 'not_ready',
                    'missing': ['runtime setup output'],
                    'setup_message': (proc.stderr or 'runtime setup produced no output').strip(),
                    'fr3d_python_dir': '',
                    'python_executable': '',
                    'output_dir': self.fr3d_output_path,
                    'cache_dir': str((Path(self.fr3d_runtime_dir) / 'cache').resolve()),
                }
            return json.loads(output)
        except Exception as e:
            return {
                'ok': False,
                'runtime_mode': 'not_ready',
                'missing': ['runtime setup execution'],
                'setup_message': f'Failed to run setup script: {type(e).__name__}: {e}',
                'fr3d_python_dir': '',
                'python_executable': '',
                'output_dir': self.fr3d_output_path,
                'cache_dir': str((Path(self.fr3d_runtime_dir) / 'cache').resolve()),
            }

    def ensure_fr3d_runtime_ready(self, auto_setup: bool = True) -> bool:
        """Ensure integrated Source-5 runtime is available for local pipeline execution."""
        report = self._run_fr3d_runtime_setup(build=False)
        if auto_setup and not self.fr3d_setup_attempted:
            self.fr3d_setup_attempted = True
            report = self._run_fr3d_runtime_setup(build=True)

        setup_message = str(report.get('setup_message', '') or '').strip()
        if setup_message:
            self.logger.info(setup_message)

        self.fr3d_root_path = str(report.get('fr3d_python_dir', '') or self.fr3d_root_path).strip()
        self.fr3d_python_path = str(report.get('python_executable', '') or self.fr3d_python_path).strip()
        output_dir = str(report.get('output_dir', '') or '').strip()
        if output_dir:
            self.fr3d_output_path = output_dir

        if report.get('ok'):
            return True

        missing = report.get('missing', []) or []
        if missing:
            self.logger.error(f"FR3D runtime is incomplete: {', '.join(str(x) for x in missing)}")
        self.logger.info('Run: rmv_fr3d doctor')
        return False

    def fr3d_doctor(self, auto_setup: bool = False):
        """Print integrated Source-5 runtime diagnostics."""
        report = self._run_fr3d_runtime_setup(build=auto_setup)

        print('\n' + '=' * 70)
        print('FR3D Doctor')
        print('=' * 70)
        print(f"Runtime dir     : {self.fr3d_runtime_dir}")
        print(f"Runtime mode    : {report.get('runtime_mode', 'not_ready')}")
        print(f"FR3D root       : {report.get('fr3d_python_dir', '(not found)') or '(not found)'}")
        print(f"Python          : {report.get('python_executable', '(not found)') or '(not found)'}")
        print(f"Output dir      : {report.get('output_dir', self.fr3d_output_path)}")
        print(f"Cache dir       : {report.get('cache_dir', '(not set)')}")
        print(f"Status          : {'LOCAL PIPELINE READY' if report.get('ok') else 'NOT READY'}")

        missing = report.get('missing', []) or []
        if missing:
            print('\nMissing local runtime items:')
            for item in missing:
                print(f"  - {item}")

        setup_message = str(report.get('setup_message', '') or '').strip()
        if setup_message:
            print(f"\nSetup message: {setup_message}")

        print('\nHints:')
        print('  - Normal workflow does not require manual config: rmv_db 5 -> rmv_fetch -> rmv_load_motif')
        print('  - Source 5 runs local FR3D pipeline only (no BGSU fallback)')
        print('  - Release packaging expects a bundled FR3D runtime under rsmviewer/tools/fr3d_runtime/python/fr3d-python')
        print('  - Python dependencies (numpy, scipy, mmcif-pdbx) must be available in the normal plugin environment')
        print('=' * 70 + '\n')

    def _load_fr3d_wrapper_config(self):
        """Load persisted FR3D wrapper settings if available."""
        try:
            cfg_path = self.fr3d_config_file
            if not cfg_path.exists():
                return

            with open(cfg_path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            self.fr3d_root_path = str(data.get('fr3d_root_path', '') or '').strip()
            self.fr3d_input_path = str(data.get('fr3d_input_path', '') or '').strip()
            self.fr3d_python_path = str(data.get('fr3d_python_path', '') or '').strip()
            self.fr3d_auto_run_on_fetch = bool(data.get('fr3d_auto_run_on_fetch', False))
            output_path = str(data.get('fr3d_output_path', '') or '').strip()
            if output_path:
                self.fr3d_output_path = output_path
        except Exception as e:
            self.logger.warning(f"Could not read FR3D wrapper config: {e}")

    def _save_fr3d_wrapper_config(self):
        """Persist FR3D wrapper settings."""
        try:
            self.fr3d_config_file.parent.mkdir(parents=True, exist_ok=True)
            payload = {
                'fr3d_root_path': self.fr3d_root_path,
                'fr3d_python_path': self.fr3d_python_path,
                'fr3d_input_path': self.fr3d_input_path,
                'fr3d_output_path': self.fr3d_output_path,
                'fr3d_auto_run_on_fetch': self.fr3d_auto_run_on_fetch,
            }
            with open(self.fr3d_config_file, 'w', encoding='utf-8') as f:
                json.dump(payload, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Could not save FR3D wrapper config: {e}")

    def print_fr3d_wrapper_status(self):
        """Print current FR3D wrapper settings and quick usage hints."""
        report = self._run_fr3d_runtime_setup(build=False)
        print("\n" + "=" * 70)
        print("FR3D Wrapper Status")
        print("=" * 70)
        print(f"Runtime dir    : {self.fr3d_runtime_dir}")
        print(f"Runtime mode   : {report.get('runtime_mode', 'not_ready')}")
        print(f"FR3D tool path : {report.get('fr3d_python_dir', self.fr3d_root_path) or '(auto-detect)'}")
        print(f"Python path    : {report.get('python_executable', self.fr3d_python_path) or '(auto-detect)'}")
        print(f"Output path    : {self.fr3d_output_path}")
        print(f"Auto on fetch  : {'on' if self.fr3d_auto_run_on_fetch else 'off'}")
        print("\nCommands:")
        print("  rmv_fr3d status               Show runtime status")
        print("  rmv_fr3d doctor               Diagnose runtime readiness")
        print("  rmv_fr3d setup                Run first-time runtime setup")
        print("  rmv_fr3d refresh [PDB_ID]     Force rerun for current/specified PDB")
        print("  rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]  (advanced)")
        print("  rmv_fr3d run <PDB_ID>         Run local FR3D pipeline")
        print("  rmv_fr3d run_current          Run for the currently loaded PDB")
        print("\nNormal workflow:")
        print("  rmv_db 5 -> rmv_fetch <PDB_ID> -> rmv_load_motif")
        print("  (always runs fresh pipeline; no cached FR3D result reuse)")
        print("=" * 70 + "\n")

    def configure_fr3d_wrapper(self, fr3d_root: str, output_dir: str = '', input_dir: str = '', python_executable: str = ''):
        """Set FR3D wrapper settings and persist them.

        Current motif pipeline uses BGSU loops download; FR3D root path is kept
        as an optional compatibility setting for users who want a configured
        local FR3D location in the plugin settings.
        """
        if fr3d_root:
            expanded_root = os.path.abspath(os.path.expanduser(fr3d_root))
            if not os.path.isdir(expanded_root):
                self.logger.error(f"FR3D root directory not found: {expanded_root}")
                return False
            self.fr3d_root_path = expanded_root

        if output_dir:
            self.fr3d_output_path = os.path.abspath(os.path.expanduser(output_dir))

        if input_dir:
            auto_val = input_dir.strip().lower()
            if auto_val in ('1', 'on', 'true', 'yes'):
                self.fr3d_auto_run_on_fetch = True
            elif auto_val in ('0', 'off', 'false', 'no'):
                self.fr3d_auto_run_on_fetch = False
            else:
                self.logger.error("AUTO_ON_FETCH must be one of: on/off, true/false, 1/0")
                return False

        try:
            os.makedirs(self.fr3d_output_path, exist_ok=True)
        except Exception as e:
            self.logger.error(f"Could not create FR3D output directory {self.fr3d_output_path}: {type(e).__name__}: {e}")
            return False

        self._save_fr3d_wrapper_config()

        self.logger.success("FR3D wrapper configuration updated")
        self.logger.info(f"FR3D tool path: {self.fr3d_root_path or '(not set)'}")
        self.logger.info(f"Output path: {self.fr3d_output_path}")
        self.logger.info(f"Auto run on rmv_fetch: {'on' if self.fr3d_auto_run_on_fetch else 'off'}")
        return True
    def _python_has_scipy(self, python_executable: str) -> bool:
        """Check whether a candidate Python executable can import scipy."""
        if not python_executable or not os.path.isfile(python_executable):
            return False

        try:
            result = subprocess.run(
                [python_executable, '-c', 'import scipy'],
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            return result.returncode == 0
        except Exception:
            return False

    def _auto_detect_fr3d_python(self) -> str:
        """Try to find a Python executable that can import scipy."""
        candidates = []

        if self.fr3d_root_path:
            candidates.append(str(Path(self.fr3d_root_path) / '.fr3d-venv' / 'bin' / 'python'))

        if self.fr3d_python_path:
            candidates.append(self.fr3d_python_path)
        candidates.extend([
            '/opt/homebrew/bin/python3',
            '/usr/local/bin/python3',
            sys.executable,
        ])

        seen = set()
        for candidate in candidates:
            candidate = os.path.abspath(os.path.expanduser(candidate))
            if candidate in seen:
                continue
            seen.add(candidate)
            self.logger.debug(f"Checking FR3D Python candidate: {candidate}")
            if self._python_has_fr3d_deps(candidate):
                self.logger.debug(f"Selected FR3D Python interpreter: {candidate}")
                return candidate

        self.logger.error("No FR3D Python interpreter with scipy and pdbx was found.")
        self.logger.info("Tried: local FR3D venv, /opt/homebrew/bin/python3, /usr/local/bin/python3, and the current interpreter.")
        return ''

    def _auto_detect_fr3d_root(self) -> str:
        """Try to auto-detect a sibling fr3d-python checkout."""
        here = Path(__file__).resolve()
        candidates = [
            here.parents[2] / 'fr3d-python',
            here.parents[1] / 'fr3d-python',
        ]
        for candidate in candidates:
            if (candidate / 'fr3d' / '__init__.py').exists():
                return str(candidate)
        return ''

    def run_fr3d_wrapper(self, pdb_id: str, category: str = 'motif', force_refresh: bool = False):
        """Run FR3D source-5 pipeline and load motifs.

                Source-5 policy:
                    - Always run local FR3D pipeline fresh.
                    - Never reuse cached FR3D loop outputs.
                    - Never fall back to BGSU loop download.
        """
        pdb_upper = str(pdb_id).strip().upper()
        if not pdb_upper:
            self.logger.error("PDB ID is required")
            self.logger.info("Usage: rmv_fr3d run <PDB_ID>")
            return False

        try:
            os.makedirs(self.fr3d_output_path, exist_ok=True)
        except Exception as e:
            self.logger.error(f"Could not create FR3D output directory {self.fr3d_output_path}: {type(e).__name__}: {e}")
            return False

        out_file = Path(self.fr3d_output_path) / f"{pdb_upper}_fr3d_loops.csv"

        # Source-5 always runs fresh pipeline and does not reuse cached outputs.
        if out_file.exists():
            try:
                out_file.unlink()
                self.logger.info(f"Removed previous FR3D output: {out_file.name}")
            except Exception as e:
                self.logger.warning(f"Could not remove previous FR3D output {out_file.name}: {e}")

        # -- Local FR3D pipeline --------------------------------------------
        if not self.ensure_fr3d_runtime_ready(auto_setup=True):
            return False
        runtime_report = self._run_fr3d_runtime_setup(build=False)

        pipeline_cfg = self._build_internal_fr3d_config()
        pipeline_cfg['fr3d_python_dir'] = str(runtime_report.get('fr3d_python_dir', '') or self.fr3d_root_path).strip()
        pipeline_cfg['python_executable'] = str(runtime_report.get('python_executable', '') or self.fr3d_python_path or sys.executable).strip()
        pipeline_cfg['cache_dir'] = str(runtime_report.get('cache_dir', pipeline_cfg.get('cache_dir', '')) or '').strip()
        pipeline_cfg['force_fresh'] = True

        fr3d_dir = str(pipeline_cfg.get('fr3d_python_dir', '') or '').strip()
        if not (fr3d_dir and os.path.isdir(fr3d_dir)):
            self.logger.error("FR3D local runtime not found. Source 5 does not use BGSU fallback.")
            self.logger.info("Use: rmv_fr3d doctor")
            return False

        self.logger.info(f"Running local FR3D pipeline fresh for {pdb_upper}...")
        csv_path = self._run_fr3d_local_pipeline(pdb_upper, pipeline_cfg)
        if not csv_path or not os.path.isfile(csv_path):
            self.logger.error("Local FR3D pipeline failed and fallback is disabled.")
            return False

        # Ensure source 5 points to this output directory and load directly.
        self.user_data_paths[5] = self.fr3d_output_path
        self._handle_source_by_id(5, self.fr3d_output_path)
        self.load_user_annotations_action('fr3d', pdb_upper, auto_pipeline=False)

        loaded_motifs = self.viz_manager.motif_loader.get_loaded_motifs()
        if loaded_motifs:
            total_instances = sum(len(info.get('motif_details', [])) for info in loaded_motifs.values())
            self.logger.success(f"Loaded FR3D motifs (fresh local pipeline): {len(loaded_motifs)} motif types, {total_instances} instances")
            return True

        self.logger.warning("FR3D pipeline generated output, but no motifs were loaded into RSMViewer.")
        return False

    def _run_fr3d_local_pipeline(self, pdb_id: str, config: dict) -> str:
        """Run the local FR3D loop-extraction pipeline via fr3d_loop_extractor.py.

        This implements the RNA 3D Motif Atlas pipeline locally:
          1. Run NA_pairwise_interactions.py to annotate canonical cWW basepairs
          2. Apply flankSS logic to extract hairpin (HL) and internal (IL) loops
          3. Write {PDB_ID}_fr3d_loops.csv in BGSU format to the configured output dir

        Args:
            pdb_id : PDB ID (uppercase)
            config : Parsed fr3d_pipeline_config.json dict

        Returns:
            Path to the output CSV file on success, empty string on failure.
        """
        try:
            # Import the pipeline module from rsmviewer/tools/
            tools_dir = str(Path(__file__).parent / 'tools')
            if tools_dir not in sys.path:
                sys.path.insert(0, tools_dir)

            from fr3d_loop_extractor import run_pipeline  # type: ignore

            # Resolve output directory (may have been set from config earlier)
            out_dir = str(config.get('output_dir', '') or '').strip()
            if not out_dir:
                out_dir = self.fr3d_output_path
            config_with_outdir = dict(config)
            config_with_outdir['output_dir'] = out_dir

            csv_path = run_pipeline(config_with_outdir, pdb_id)
            return csv_path or ''

        except ImportError as e:
            self.logger.error(f"Could not import FR3D pipeline module: {e}")
            return ''
        except Exception as e:
            self.logger.error(f"Local FR3D pipeline error: {type(e).__name__}: {e}")
            return ''

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
        
        Handles combine mode by loading from multiple sources and merging.
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
                            self.logger.warning(f"Mismatch: source_override={source_override}, config.specific_source={config.specific_source}")
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
            
            # Use pre-built auth->label chain mapping (from CIF parsing in rmv_fetch)
            # When cif_use_auth=0, motif data has auth_asym_id chains but PyMOL has label_asym_id
            auth_to_label = self.auth_to_label_map if self.cif_use_auth == 0 else {}
            
            # Count total motifs
            total_count = sum(len(instances) for instances in available_motifs.values())
            self.logger.success(f"Found {total_count} motifs in {pdb_id} (source: {source_name})")
            
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
                        # Rebuild combined selection
                        if main_motif_sel:
                            if existing['main_selection']:
                                existing['main_selection'] = f"{existing['main_selection']} or {main_motif_sel}"
                            else:
                                existing['main_selection'] = main_motif_sel
                        self.logger.success(f"Loaded {len(motif_details)} more {display_type_upper} motifs (total: {existing['count']})")
                    else:
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
                        }
                        self.logger.success(f"Loaded {len(motif_details)} {display_type_upper} motifs")
            
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
                tag = f"{pdb_id_upper}{source_suffix}" if source_suffix else pdb_id_upper
                self.logger.success(f"Loaded {len(motif_summary)} motif types from {pdb_id} (tag: {tag})")
                self.logger.info("")
                self.logger.info("Motif data ready (not rendered, no objects created)")
                self.logger.info("Next steps:")
                self.logger.info(f"  rmv_summary              Show all motifs")
                self.logger.info(f"  rmv_summary HL           Show HL instances")
                self.logger.info(f"  rmv_view all             Highlight all motifs on structure")
                self.logger.info(f"  rmv_show HL              Render hairpin loops")
                self.logger.info(f"  rmv_show HL 1            Zoom to specific instance")
                self.logger.info("")
            else:
                self.logger.warning(f"No valid motifs found for {pdb_id}")
                
        except Exception as e:
            self.logger.error(f"Failed to load motif data: {str(e)}")

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
        self.viz_manager.cmd.color('gray80', f"model {structure_name} and polymer.nucleic")

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
    
    def _load_combined_motifs(self, pdb_id: str, source_ids: List[int]):
        """
        Load, enrich, and cascade-merge motifs from multiple sources.
        
        Pipeline:
        1. Fetch raw motifs from each source
        2. Enrich generic sources (1,3,5) via homolog representative lookup
        3. Cascade merge (right-to-left, priority = source order)
        
        Args:
            pdb_id: PDB ID to fetch
            source_ids: List of source IDs in priority order (first = highest)
        
        Returns:
            Merged motif dictionary: {motif_type: [MotifInstance, ...]}
        """
        try:
            from .database.config import SOURCE_ID_MAP
            from .database.cascade_merger import CascadeMerger
            from .database.homolog_enricher import HomologEnricher
            from .database.representative_set import get_representative_loader
            
            pdb_id = pdb_id.upper()
            
            # Sources that use generic names and need enrichment
            GENERIC_SOURCES = {1, 3, 5}  # Atlas, BGSU, FR3D
            
            # --- Step 1: Fetch raw motifs from each source ---
            self.logger.info(f"Step 1/3: Fetching motifs from {len(source_ids)} sources...")
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
                    self.logger.info(f"  [{sid}] {label}: {total} motifs in {len(motifs)} categories")
                else:
                    raw_sources[sid] = {}
                    self.logger.warning(f"  [{sid}] {label}: no motifs found")
            
            if not any(raw_sources.values()):
                self.logger.error("No motifs found from any source")
                return {}
            
            # --- Step 2: Enrich generic sources via homolog ---
            generic_in_selection = [sid for sid in source_ids if sid in GENERIC_SOURCES and raw_sources.get(sid)]
            
            # Pre-filter: skip sources that already have specific (non-generic)
            # names for ALL their categories (e.g., source 3 after successful
            # HTML scraping).  Only enrich sources that still carry generic names.
            from .database.homolog_enricher import _is_generic_name
            sources_needing_enrichment = [
                sid for sid in generic_in_selection
                if any(_is_generic_name(cat) for cat in raw_sources[sid])
            ]
            skipped = set(generic_in_selection) - set(sources_needing_enrichment)
            if skipped:
                self.logger.info(
                    f"  Skipping enrichment for source(s) {skipped} "
                    f"(already have specific names)")

            if sources_needing_enrichment:
                self.logger.info(f"Step 2/3: Enriching {len(sources_needing_enrichment)} generic source(s) via homolog lookup...")
                try:
                    rep_loader = get_representative_loader()
                    
                    # Get or create BGSU provider for fetching representative annotations
                    from .database import get_source_selector
                    source_selector = get_source_selector()
                    bgsu_provider = None
                    if source_selector and hasattr(source_selector, 'providers'):
                        bgsu_provider = source_selector.providers.get('bgsu_api')
                    
                    if not bgsu_provider:
                        # Create a standalone BGSU provider
                        from .database.bgsu_api_provider import BGSUAPIProvider
                        from .database import get_cache_manager
                        bgsu_provider = BGSUAPIProvider(cache_manager=get_cache_manager())
                    
                    enricher = HomologEnricher(rep_loader, bgsu_provider)
                    
                    for sid in sources_needing_enrichment:
                        if raw_sources[sid]:
                            before_cats = len(raw_sources[sid])
                            raw_sources[sid] = enricher.enrich(pdb_id, raw_sources[sid])
                            after_cats = len(raw_sources[sid])
                            self.logger.info(f"  [{sid}] Enriched: {before_cats} -> {after_cats} categories")
                            
                except Exception as e:
                    self.logger.warning(f"Enrichment failed (using generic names): {e}")
                    import traceback
                    traceback.print_exc()
            else:
                self.logger.info("Step 2/3: No generic sources to enrich (all sources have specific names)")
            
            # --- Step 2.5: Stamp source origin on each instance ---
            for sid in source_ids:
                info = SOURCE_ID_MAP.get(sid, {})
                label = info.get('name', f'Source {sid}')
                for mtype, instances in raw_sources.get(sid, {}).items():
                    for inst in instances:
                        if inst.metadata is None:
                            inst.metadata = {}
                        inst.metadata['_source_id'] = sid
                        inst.metadata['_source_label'] = label
            
            # --- Step 2.7: Within-source deduplication ---
            # Some tools (e.g., RMSX) can produce duplicate entries with
            # identical residue sets. Remove exact duplicates within each
            # source before cascade merge (which only deduplicates ACROSS
            # sources).
            self.dedup_stats = {}
            for sid in source_ids:
                src_data = raw_sources.get(sid, {})
                if not src_data:
                    continue
                total_before = sum(len(v) for v in src_data.values())
                deduped = {}
                for mtype, instances in src_data.items():
                    seen = set()
                    unique = []
                    for inst in instances:
                        # Build a hashable key from (chain, residue_number) set
                        rset = frozenset(
                            (r.chain, r.residue_number) for r in inst.residues
                        )
                        if rset not in seen:
                            seen.add(rset)
                            unique.append(inst)
                    deduped[mtype] = unique
                total_after = sum(len(v) for v in deduped.values())
                self.dedup_stats[sid] = (total_before, total_after)
                if total_before != total_after:
                    self.logger.info(
                        f"  [{sid}] Within-source dedup: {total_before} -> {total_after} "
                        f"(removed {total_before - total_after} duplicates)")
                raw_sources[sid] = deduped
            
            # --- Step 3: Cascade merge ---
            self.logger.info(f"Step 3/3: Cascade merging {len(source_ids)} sources...")
            
            # Build ordered list matching source_ids order
            ordered_sources = [raw_sources.get(sid, {}) for sid in source_ids]
            ordered_labels = [SOURCE_ID_MAP.get(sid, {}).get('name', f'Source {sid}') for sid in source_ids]
            
            merger = CascadeMerger(jaccard_threshold=self.jaccard_threshold)
            merged = merger.merge_sources(ordered_sources, ordered_labels)
            
            if merged:
                total = sum(len(v) for v in merged.values())
                self.logger.success(f"Merge complete: {total} motifs in {len(merged)} categories")
            
            return merged
            
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
                    web_map = {'bgsu': 'bgsu_api', 'rfam_api': 'rfam_api'}
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
    
    def load_user_annotations_action(self, tool, pdb_id, auto_pipeline: bool = True, force_pipeline_refresh: bool = False):
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

            # -- FR3D auto-run pipeline (always fresh) ------------------------
            # Source 5 intentionally mirrors RMSX behavior:
            # - Always run the local pipeline on rmv_load_motif.
            # - Do not reuse cached FR3D outputs.
            # - Do not use BGSU fallback.
            if tool_lower == 'fr3d' and auto_pipeline:
                pdb_upper = pdb_id.strip().upper()
                self.logger.info(
                    f"Running FR3D pipeline fresh for {pdb_upper} (no cache, no fallback)..."
                )
                self.run_fr3d_wrapper(pdb_upper, force_refresh=True)
                # run_fr3d_wrapper already loads the motifs if it succeeds,
                # so we can return here to avoid double-loading.
                return

            # -- RMSX pipeline run (prebuilt/cache-aware by default) -----------
            # Source 7 reuses available/prebuilt results unless explicitly
            # forced to run fresh.
            elif tool_lower in ['rmsx', 'rnamotifscanx']:
                if not self.ensure_rmsx_runtime_ready(auto_setup=True):
                    return
                rmsx_cfg = getattr(self, 'rmsx_pipeline_config', {}) or self._build_internal_rmsx_config()
                self.rmsx_pipeline_config = dict(rmsx_cfg)
                if rmsx_cfg:
                    try:
                        tools_dir = str(Path(__file__).parent / 'tools')
                        if tools_dir not in sys.path:
                            sys.path.insert(0, tools_dir)
                        from rmsx_runner import run_pipeline as rmsx_run  # type: ignore
                        pdb_upper = pdb_id.strip().upper()
                        force_fresh = bool(force_pipeline_refresh)

                        run_cfg = dict(rmsx_cfg)
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
                        self.logger.info(
                            f"Running RMSX pipeline for {pdb_upper} ({mode_text}, no external download) in {out_dir}..."
                        )
                        existing = rmsx_run(run_cfg, pdb_upper, force_fresh=force_fresh)

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
                    self.logger.info(f"Using custom data path for {tool.upper()}: {_udp}")
            
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
                    self.logger.info(f"Checked FR3D annotation path: {_udp or user_annotations_dir / 'fr3d'}")
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
            
            # Process motifs (same as fetch_motif_data_action)
            motif_summary = {}
            from .utils.parser import SelectionParser
            
            total_count = sum(len(instances) for instances in available_motifs.values())
            self.logger.success(f"Found {total_count} motifs in {pdb_id} (source: {tool.upper()})")
            
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
                    self.logger.success(f"Loaded {len(motif_details)} {display_type_upper} motifs")
            
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
                self.logger.success(f"Loaded {len(motif_summary)} motif types from {tool.upper()}")
                self.logger.info("")
                self.logger.info("Motif data ready (not rendered)")
                self.logger.info("Next steps:")
                self.logger.info(f"  rmv_summary              Show all motifs")
                self.logger.info(f"  rmv_summary <TYPE>       Show specific motif type")
                self.logger.info(f"  rmv_view all             Highlight all motifs on structure")
                self.logger.info(f"  rmv_show <TYPE>          Render motif on structure")
                self.logger.info(f"  rmv_show <TYPE> <NO>     Zoom to specific instance")
                self.logger.info("")
            
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
                visible_str = "✓ visible" if data['visible'] else "✗ hidden"
                print(f"  {motif_type:20s} ({data['count']:2d} instances) {visible_str}")
        else:
            print("\nNo motifs loaded for this structure")
        
        print("="*60 + "\n")
    
    def print_sources(self):
        """Print available data sources with ID numbers - new format."""
        print("\n" + "="*80)
        print("  [DB]  AVAILABLE DATA SOURCES")
        print("="*80)
        
        try:
            from .database import get_config
            config = get_config()
            
            active_source = config.source_mode.value.upper()
            if active_source == 'BGSU':
                active_source += ' (BGSU RNA 3D Hub API)'
            elif active_source == 'AUTO':
                active_source += ' (Auto-select)'
            elif active_source == 'LOCAL':
                active_source += ' (Local databases)'
            
            print(f"\n  Currently Active: {active_source}\n")
            
            # Display sources with IDs grouped by category
            current_category = None
            for source_id in sorted(SOURCE_ID_MAP.keys()):
                info = SOURCE_ID_MAP[source_id]
                category = info.get('category')
                
                if category != current_category:
                    if current_category is not None:
                        print()
                    print(f"  {category}:")
                    print("  " + "-"*76)
                    current_category = category
                
                # Format: ID | Source Name | Command
                source_name = info.get('name', 'Unknown')
                cmd_str = f"rmv_db {source_id}"
                if info.get('supports_filtering'):
                    cmd_str += " [on|off|custom]" 
                
                print(f"    [{source_id}] {source_name:<24} {cmd_str}")
            
            print("\n" + "="*80)
            print("  COMMAND EXAMPLES")
            print("="*80)
            print("""
    Basic Usage:
      rmv_db 1                         Select RNA 3D Motif Atlas
      rmv_db 3                         Select BGSU RNA 3D Hub API
      rmv_db 6                         Select RNAMotifScan (RMS)
    
    Multi-Source Combine:
      rmv_db 1 3                       Combine Atlas + BGSU (Atlas = priority)
      rmv_db 2 5 3                     Combine 3 sources (Rfam = highest priority)
      -> Includes homolog enrichment & cascade merge
    
    With Filtering Control (RMS/RMSX/NoBIAS):
      rmv_db 6 off                     Show all motifs (no filtering)
      rmv_db 7 on                      Apply default P-value cutoffs
      rmv_db 6 C-LOOP 0.05 KINK-TURN 0.02   Custom P-values
    
    Combine with P-values (2-step workflow):
      Step 1: rmv_db 6 C-LOOP 0.05 KINK-TURN 0.02   Set P-values
      Step 2: rmv_db 1 6               Combine Atlas + RMS (P-values apply)
    
    After selecting a source:
      rmv_load_motif                   Fetch motif data
      rmv_summary                      Show motif types & counts
      rmv_summary HL                   Show hairpin loop instances
      rmv_show HL                      Render all hairpin loops
      rmv_show HL 1                    Zoom to specific instance
""")
            print("="*80)
            print("  MORE INFO")
            print("="*80)
            print("""
    Get detailed info about a source:
      rmv_source info 1                Show detailed info for source 1
    
    Get help:
      rmv_help                         List all available commands
""")
        except Exception as e:
            print(f"  Error loading sources: {e}")
        
        print("="*80 + "\n")
    
    def print_help(self):
        """Print all available commands in box format."""
        print("\n" + "+" + "-"*78 + "+")
        print("|" + "  RSMViewer v1.0.0 - COMMAND REFERENCE".center(78) + "|")
        print("|" + "  Last Updated: 26 July 2026".center(78) + "|")
        print("+" + "-"*78 + "+\n")
        
        print("+" + "-"*78 + "+")
        print("|  [DB] DATABASE SELECTION                                                  |")
        print("+" + "-"*78 + "+")
        print("|  rmv_sources               List all available data sources (1-8)        |")
        print("|  rmv_db <N>                Select motif data source by ID (1-8)          |")
        print("|  rmv_db <N> <N>            Combine multiple sources (e.g., 8 7)         |")
        print("|  rmv_db <N> /path/to/data  Use custom data path (sources 5-8)           |")
        print("|  rmv_db <N> <N>, jaccard_threshold=0.80  Set merge threshold            |")
        print("+" + "-"*78 + "+")
        print("|  [LD] LOADING & DATA                                                      |")
        print("+" + "-"*78 + "+")
        print("|  rmv_fetch <PDB_ID>        Load PDB structure (no motif data)            |")
        print("|  rmv_fetch /path/to/file   Load local PDB or mmCIF file                  |")
        print("|  rmv_fetch <ID> cif_use_auth=0  Load with label_asym_id chains          |")
        print("|  rmv_load_motif            Fetch motif data from selected source         |")
        print("|  rmv_load <PDB_ID>         Legacy: show workflow guide                  |")
        print("|  rmv_refresh               Force refresh cache and collect again       |")
        print("+" + "-"*78 + "+")
        print("|  [VZ] VISUALIZATION                                                       |")
        print("+" + "-"*78 + "+")
        print("|  rmv_show ALL              Show all motif types with objects            |")
        print("|  rmv_show <TYPE>           Highlight all instances of a motif type      |")
        print("|  rmv_show <TYPE> <NO>      Zoom to specific instance with details       |")
        print("|  rmv_show <TYPE> 1,3,5     Show multiple specific instances             |")
        print("|  rmv_show <TYPE>, padding=N  Expand residue ranges +/-N for visualization |")
        print("|  rmv_view all              Highlight all motif regions on structure     |")
        print("|  rmv_view <TYPE>           Zoom to motif regions (no objects created)   |")
        print("|  rmv_view <TYPE> <NO>      Zoom to instance & create selection          |")
        print("|  rmv_view hide             Reset all view coloring to gray              |")
        print("|  rmv_view <TYPE> hide      Reset coloring for a specific motif type     |")
        print("|  rmv_toggle <TYPE> on/off  Toggle motif visibility                      |")
        print("|  rmv_bg_color <COLOR>      Change background (non-motif) color          |")
        print("|  rmv_color <TYPE> <COLOR>  Change motif color                           |")
        print("|  rmv_colors                Show color legend for motif types             |")
        print("+" + "-"*78 + "+")
        print("|  [SV] SAVE & EXPORT COMMANDS                                              |")
        print("+" + "-"*78 + "+")
        print("|  rmv_save ALL [rep]        Save all motif instance images to disk        |")
        print("|  rmv_save <TYPE> [rep]     Save all instances of specific motif type    |")
        print("|  rmv_save <TYPE> <NO> [rep]Save specific motif instance by ID           |")
        print("|                            (Saves MOTIF ONLY - no background structure) |")
        print("|  rmv_save current          Save current PyMOL view (like png command)   |")
        print("|  rmv_save current <FILE>   Save current view to specific filename       |")
        print("|  [rep]: cartoon (default), sticks, spheres, ribbon, lines, etc.         |")
        print("|                            (Output: plugin_dir/motif_images/pdb_id/)    |")
        print("|                                                                          |")
        print("|  mmCIF Structure Export (original coordinates from disk):                |")
        print("|  rmv_save ALL cif          Export all motif structures as mmCIF          |")
        print("|  rmv_save <TYPE> cif       Export all instances of type as mmCIF         |")
        print("|  rmv_save <TYPE> <NO> cif  Export specific instance as mmCIF             |")
        print("|                            (Output: plugin_dir/motif_structures/pdb_id/)|")
        print("+" + "-"*78 + "+")
        print("|  [IF] INFORMATION COMMANDS                                                |")
        print("+" + "-"*78 + "+")
        print("|  rmv_summary               Show all motif types & counts                 |")
        print("|  rmv_summary <TYPE>        Show instances of specific type               |")
        print("|  rmv_summary <TYPE> <NO>   Show specific instance details                |")
        print("|  rmv_source info           Show currently selected source                 |")
        print("|  rmv_source info <N>       Show detailed info about source N              |")
        print("|  rmv_chains [structure]    Show chain ID diagnostics (auth/label mapping)|")
        print("|  rmv_loaded                Show loaded PDB+source combination tags       |")
        print("|  rmv_reset                 Reset plugin: delete all objects & clear state |")
        print("|  rmv_help                  Show this command reference                   |")
        print("+" + "-"*78 + "+")
        print("|  [AN] USER ANNOTATIONS (Sources 5-8)                                      |")
        print("+" + "-"*78 + "+")
        print("|  rmv_user <TOOL> <PDB_ID>  Load FR3D/RMS/RMSX/NoBIAS annotations       |")
        print("|  rmv_fr3d status           Show FR3D wrapper configuration              |")
        print("|  rmv_fr3d config <ROOT>    Configure local FR3D absolute path(s)        |")
        print("|  rmv_fr3d run <PDB_ID>     Download FR3D structural motifs (BGSU)       |")
        print("|  rmv_rmsx status           Show integrated RNAMotifScanX runtime status |")
        print("|  rmv_rmsx_doctor           Diagnose RMSX runtime & dependencies          |")
        print("|  rmv_rmsx setup            Attempt automatic RMSX runtime setup          |")
        print("|  rmv_rmsx run <PDB_ID>     Run integrated RNAMotifScanX pipeline (fresh)|")
        print("|  rmv_db <N> /path/to/data  Use custom data directory (any source 1-8)   |")
        print("|  rmv_db 7 off              Disable P-value filtering                    |")
        print("|  rmv_db 7 on               Enable P-value filtering                     |")
        print("|  rmv_db 7 MOTIF 0.01       Set custom P-value threshold for motif       |")
        print("+" + "-"*78 + "+")
        print("|  [SP] STRUCTURAL SUPERIMPOSITION                                          |")
        print("+" + "-"*78 + "+")
        print("|  rmv_super <TYPE>              Medoid superimposition (current PDB+src)   |")
        print("|  rmv_super <TYPE>, PDB1_SN, PDB2_SN  Cross-PDB superimposition          |")
        print("|  rmv_super <TYPE> 1,3,5        Superimpose specific instances only        |")
        print("|  rmv_super <TYPE>, padding=N   Expand residue ranges +/-N for objects      |")
        print("|  rmv_align <TYPE>              Same as rmv_super but sequence-dependent   |")
        print("+" + "-"*78 + "+")
        
        print("\n  QUICK EXAMPLES:")
        print("  ---------------")
        print("  1. Standard workflow (recommended):")
        print("     rmv_fetch 1S72            # Load PDB structure")
        print("     rmv_sources               # Check available data sources")
        print("     rmv_db 3                  # Select BGSU API")
        print("     rmv_load_motif            # Fetch motif data")
        print("     rmv_summary               # Show motif types & counts")
        print("     rmv_view all              # Highlight all motifs on structure")
        print("     rmv_summary HL            # Show hairpin loop instances")
        print("     rmv_show HL               # Render all hairpin loops")
        print("     rmv_show HL 1             # Zoom to specific instance")
        print()
        print("  2. Quick preview (no objects):")
        print("     rmv_view all              # Highlight all motifs (colored regions)")
        print("     rmv_view K-TURN           # Zoom to all K-TURN instances")
        print("     rmv_view K-TURN 1         # Zoom to instance + create selection")
        print("     rmv_view K-TURN hide      # Remove K-TURN coloring from structure")
        print("     rmv_view hide             # Reset all view coloring to gray")
        print()
        print("  3. Switch sources (no re-download):")
        print("     rmv_db 7                  # Switch to RMSX")
        print("     rmv_load_motif            # Fetch RMSX data (same PDB)")
        print("     rmv_show SARCIN-RICIN     # View RMSX sarcin-ricin motifs")
        print()
        print("  4. NoBIAS source with filtering:")
        print("     rmv_db 8                  # Select NoBIAS")
        print("     rmv_load_motif            # Fetch NoBIAS data")
        print("     rmv_db 8 off              # Disable P-value filtering")
        print("     rmv_db 8 7                # Combine NoBIAS + RMSX")
        print()
        print("  5. Multi-source combine:")
        print("     rmv_db 8 7                # Combine NoBIAS + RMSX")
        print("     rmv_load_motif            # Fetch, enrich, merge")
        print("     rmv_summary               # See merged motif counts")
        print("     rmv_show K-TURN           # Render (with source attribution)")
        print()
        print("  6. Cross-PDB superimposition:")
        print("     rmv_fetch 1S72; rmv_db 3; rmv_load_motif")
        print("     rmv_fetch 4V9F; rmv_load_motif")
        print("     rmv_loaded                # See all loaded tags")
        print("     rmv_super K-TURN, 1S72_S3, 4V9F_S3")
        print()
        print("  7. Export images + structures:")
        print("     rmv_save ALL               Save all motif images (PNG)")
        print("     rmv_save ALL cif           Export all structures (mmCIF)")
        print("     rmv_save HL 3 cif          Export HL instance #3")
        print("     rmv_save current           Save current view (high-res)")
        print()
        print("  8. Run FR3D locally, then visualize immediately:")
        print("     rmv_fr3d config /abs/path/to/fr3d-python")
        print("     rmv_fr3d run_current       # Uses currently loaded PDB")
        print("     rmv_summary                # Shows FR3D interaction families")
        print()


    def get_available_motifs(self):
        """Get list of available motif types for current PDB+source."""
        filtered = self._get_current_source_motifs()
        return list(filtered.keys()) if filtered else []
    
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
    
    def show_motif_summary_for_type(self, motif_type: str):
        """Print summary of a specific motif type without rendering.
        
        Args:
            motif_type (str): Motif type to show (e.g., 'HL', 'IL', 'GNRA')
        """
        motif_arg = motif_type.upper().strip()
        
        # Use filtered view (current PDB + source only)
        loaded_motifs = self._get_current_source_motifs()
        
        if not loaded_motifs:
            print("\nNo motifs loaded. Use 'rmv_fetch <PDB_ID>' first.\n")
            return
        
        motif_arg = self._resolve_loaded_motif_type(motif_arg, loaded_motifs)

        if motif_arg not in loaded_motifs:
            available = ', '.join(loaded_motifs.keys())
            print(f"\nMotif type '{motif_arg}' not loaded.")
            print(f"Available motifs: {available}\n")
            return
        
        # Get motif details from loaded_motifs (same structure as rmv_show uses)
        motif_info = loaded_motifs[motif_arg]
        motif_details = motif_info.get('motif_details', [])
        
        # Use the visualization manager's table printer (same as rmv_show)
        self.viz_manager._print_motif_instance_table(motif_arg, motif_details)
        
        # --- Source attribution report (combine mode only) ---
        self._print_source_attribution_report(motif_arg, motif_details)
        
        print("\n  Next steps:")
        print(f"    rmv_view {motif_arg}              Highlight {motif_arg} on structure")
        print(f"    rmv_show {motif_arg}              Render & create objects for {motif_arg}")
        # Source filter suggestions in combine mode
        if self.current_source_mode == 'combine' and self.combined_source_ids:
            from .database.config import SOURCE_ID_MAP
            for sid in self.combined_source_ids:
                info = SOURCE_ID_MAP.get(sid, {})
                tool = info.get('tool', '')
                subtype = info.get('subtype', '')
                name = info.get('name', f'Source {sid}')
                # Prefer tool > subtype > first word of name as the alias
                # shown in suggestions.  Must match aliases recognised by
                # ``_resolve_source_filter`` so the suggested command works.
                if tool:
                    alias = tool
                elif subtype:
                    alias = subtype
                else:
                    alias = name.split()[0].lower()
                print(f"    rmv_show {motif_arg} {alias:<16} Show {name}-only instances")
            print(f"    rmv_show {motif_arg} shared           Show shared instances")
        print(f"    rmv_summary {motif_arg} <NO>      Show details of specific instance")
        print(f"    rmv_super {motif_arg}             Superimpose all {motif_arg} instances")
        print("="*70)
        print()
    
    def _print_source_attribution_report(self, motif_type: str, motif_details: list):
        """Print unique merged-instance IDs grouped by source.

        Only prints when instances carry _source_label metadata (combine mode).
        Three categories: unique-to-source-A, unique-to-source-B, and shared
        (instances where the cascade merger detected overlap from both sources).
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
        """Print details of a specific motif instance (for rmv_summary MOTIF NO).
        
        Args:
            motif_type (str): Motif type (e.g., 'HL')
            instance_no (int): Instance number (1-indexed)
        """
        motif_arg = motif_type.upper().strip()
        
        # Use filtered view (current PDB + source only)
        loaded_motifs = self._get_current_source_motifs()
        
        if not loaded_motifs:
            print("\nNo motifs loaded. Use 'rmv_fetch <PDB_ID>' first.\n")
            return
        
        motif_arg = self._resolve_loaded_motif_type(motif_arg, loaded_motifs)

        if motif_arg not in loaded_motifs:
            available = ', '.join(loaded_motifs.keys())
            print(f"\nMotif type '{motif_arg}' not loaded.")
            print(f"Available motifs: {available}\n")
            return
        
        # Get motif details
        motif_info = loaded_motifs[motif_arg]
        motif_details = motif_info.get('motif_details', [])
        
        # Check instance number is valid
        if instance_no < 1 or instance_no > len(motif_details):
            print(f"\nInstance {instance_no} not found. Valid range: 1-{len(motif_details)}\n")
            return
        
        # Get the specific instance (1-indexed)
        detail = motif_details[instance_no - 1]
        
        # Use the visualization manager's instance detail printer
        # Don't show object name in summary (only shown in rmv_show)
        source_suffix = motif_info.get('source_suffix', '')
        pdb_id = motif_info.get('pdb_id', '')
        self.viz_manager._print_single_instance_info(
            motif_arg, instance_no, detail,
            source_suffix=source_suffix, pdb_id=pdb_id,
            show_object=False)
        
        # Next steps
        print("  Next steps:")
        print(f"    rmv_view {motif_arg} {instance_no}            Highlight on structure (no objects)")
        print(f"    rmv_show {motif_arg} {instance_no}            Render & zoom to this instance")
        if instance_no > 1:
            print(f"    rmv_show {motif_arg} {instance_no - 1}            View previous instance")
        if instance_no < len(motif_details):
            print(f"    rmv_show {motif_arg} {instance_no + 1}            View next instance")
        print(f"    rmv_show {motif_arg}               Show all {motif_arg} instances")
        print(f"    rmv_super {motif_arg}              Superimpose all {motif_arg} instances")
        print(f"    rmv_save {motif_arg} {instance_no}             Save image of this instance")
        print()
    
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
            self.logger.success(f"Motif source mode set to: {mode_display}")
            self._print_source_mode_info()
            
            # Print follow-up suggestions
            print("\n  Next steps:")
            if self.loaded_pdb_id:
                print(f"    rmv_load_motif             Fetch motif data for {self.loaded_pdb_id}")
            else:
                print(f"    rmv_fetch <PDB_ID>         Load PDB structure")
                print(f"    rmv_load_motif             Fetch motif data")
            print()
            
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
                SourceMode.AUTO: "Automatically select best available source (local first, then API)",
                SourceMode.LOCAL: "Use only local bundled databases (offline mode)",
                SourceMode.BGSU: "Use only BGSU RNA 3D Hub API (online, ~3000+ PDBs)",
                SourceMode.RFAM: "Use only Rfam API (online, named motifs)",
                SourceMode.ALL: "Combine all sources (comprehensive, may have duplicates)"
            }
            
            print(f"\nCurrent mode: {mode.value}")
            print(f"Description: {mode_descriptions.get(mode, 'Unknown')}")
            
        except ImportError:
            print("Source selector not available")
    
    def _handle_source_by_id(self, source_id: int, extra_args: str = None):
        """Handle source selection by ID number with support for custom P-values.
        
        Also detects multi-source mode when extra_args contains additional
        numeric source IDs (e.g., 'rmv_db 1 3' or 'rmv_db 2 5 3').
        
        Args:
            source_id (int): Source ID (1-8)
            extra_args (str): Optional arguments:
                - Additional source IDs for multi-source combine (e.g., "3" or "5 3")
                - "off" : disable filtering (for RMS/RMSX/NoBIAS)
                - "on"  : enable filtering (for RMS/RMSX/NoBIAS)
                - "MOTIF_NAME p_value MOTIF_NAME2 p_value2 ..." : custom P-values
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
        
        # Display what we're combining
        self.logger.success(f"Multi-source combine mode: {len(source_ids)} sources")
        for i, sid in enumerate(source_ids, 1):
            info = SOURCE_ID_MAP[sid]
            priority = "HIGHEST" if i == 1 else ("LOWEST" if i == len(source_ids) else f"#{i}")
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
            self.logger.info(f"  {i}. [{sid}] {info['name']} (priority: {priority}){pval_note}{path_note}")
        
        self.logger.info("")
        self.logger.info(f"Features: Homolog enrichment + Cascade merge (Jaccard threshold: {self.jaccard_threshold:.0%})")
        self.logger.info("Tip: Configure P-values first with 'rmv_db 6 MOTIF 0.05', then combine.")
        self.logger.info("Tip: Set Jaccard threshold with 'rmv_db 8 6, jaccard_threshold=0.80'")
        if self.loaded_pdb_id:
            self.logger.info(f"Use 'rmv_load_motif' to load, enrich, and merge data for {self.loaded_pdb_id}")
        else:
            self.logger.info("Use 'rmv_fetch <PDB_ID>' to load structure, then 'rmv_load_motif' for data")
    
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
        self.logger.success(f"Source {source_id}: {source_info['name']}")
        self.logger.info(f"Coverage: {source_info['coverage']}")
        self.logger.info(f"Type: {source_info['description']}")
        
        # VERIFICATION: Print source configuration
        self.logger.debug(f"SOURCE CONFIG VERIFICATION:")
        self.logger.debug(f"  - self.current_source_mode = {self.current_source_mode}")
        self.logger.debug(f"  - self.current_local_source = {self.current_local_source}")
        self.logger.debug(f"  - config.specific_source = {config.specific_source}")
        self.logger.debug(f"  - Expected to load from: {subtype} ONLY")
        
        self.logger.info("\nNext steps:")
        if self.loaded_pdb_id:
            self.logger.info(f"  rmv_load_motif             Fetch motif data for {self.loaded_pdb_id}")
        else:
            self.logger.info("  rmv_fetch <PDB_ID>         Load PDB structure")
            self.logger.info("  rmv_load_motif             Fetch motif data")
    
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
        subtype_to_provider = {'bgsu': 'bgsu_api', 'rfam_api': 'rfam_api'}
        provider_id = subtype_to_provider.get(subtype, subtype)
        config.specific_source = provider_id
        
        # Map subtype to SourceMode
        mode_map = {'bgsu': 'bgsu', 'rfam_api': 'rfam'}
        mode = mode_map.get(subtype, 'auto')
        self.set_source_mode(mode)
        self.logger.debug(f"Set config.specific_source = {provider_id} (from subtype={subtype})")
        
        self.logger.success(f"Source {source_id}: {source_info['name']}")
        self.logger.info(f"Coverage: {source_info['coverage']}")
        self.logger.info(f"Type: {source_info['description']}")
        
        # VERIFICATION: Print source configuration
        self.logger.debug(f"SOURCE CONFIG VERIFICATION:")
        self.logger.debug(f"  - self.current_source_mode = {self.current_source_mode}")
        self.logger.debug(f"  - self.current_web_source = {self.current_web_source}")
        self.logger.debug(f"  - config.specific_source = {config.specific_source}")
        self.logger.debug(f"  - Expected to load from: {provider_id} ONLY")
        
        self.logger.info("\nNext steps:")
        if self.loaded_pdb_id:
            self.logger.info(f"  rmv_load_motif             Fetch motif data for {self.loaded_pdb_id}")
        else:
            self.logger.info("  rmv_fetch <PDB_ID>         Load PDB structure")
            self.logger.info("  rmv_load_motif             Fetch motif data")
    
    def _handle_user_source_by_id(self, source_id: int, source_info: Dict, extra_args: str = None):
        """Handle user annotation source selection by ID.
        
        Supports:
        - rmv_db 6              (RMS with default filtering ON)
        - rmv_db 6 off          (RMS with filtering OFF)
        - rmv_db 6 on           (RMS with filtering ON - explicit)
        - rmv_db 6 C-LOOP 0.05 KINK-TURN 0.02  (RMS with custom P-values)
        - rmv_db 5 /path/to/data   (FR3D with custom data directory)
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

                # Special case: JSON pipeline config for FR3D (source 5) or RMSX (source 7)
                if (source_id == 5
                        and os.path.isfile(expanded_path)
                        and expanded_path.lower().endswith('.json')):
                    try:
                        with open(expanded_path, 'r', encoding='utf-8') as _f:
                            _cfg = json.load(_f)
                        _out = str(_cfg.get('output_dir', '') or '').strip()
                        if source_id == 5:
                            self.fr3d_pipeline_config_path = expanded_path
                            self.fr3d_pipeline_config = _cfg
                            if _out:
                                self.fr3d_output_path = os.path.abspath(os.path.expanduser(_out))
                                os.makedirs(self.fr3d_output_path, exist_ok=True)
                            self.logger.success(f"FR3D pipeline config loaded: {expanded_path}")
                            self.logger.info(f"  FR3D dir : {_cfg.get('fr3d_python_dir', '(not set)')}")
                            self.logger.info(f"  Output   : {self.fr3d_output_path}")
                            self.logger.info(f"  Python   : {_cfg.get('python_executable', 'python3')}")
                    except Exception as _e:
                        self.logger.error(f"Failed to parse pipeline config (source {source_id}): {_e}")
                elif (source_id == 7
                      and os.path.isfile(expanded_path)
                      and expanded_path.lower().endswith('.json')):
                    self.logger.warning(
                        "Source 7 now uses integrated runtime config; external JSON config is ignored."
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
                parts = extra_str.split()
                
                # Check if first arg is on/off
                if len(parts) > 0 and parts[0].lower() in ['on', 'off']:
                    filtering_enabled = (parts[0].lower() == 'on')
                    parts = parts[1:]  # Remove the on/off argument
                
                # Parse remaining arguments as custom P-values (MOTIF_NAME p_value pairs)
                if len(parts) > 0:
                    i = 0
                    while i < len(parts):
                        if i + 1 < len(parts):
                            try:
                                pvalue = float(parts[i + 1])
                                motif_name = self._normalize_filter_motif_name(parts[i])
                                custom_pvalues[motif_name] = pvalue
                                i += 2
                            except ValueError:
                                i += 1
                        else:
                            i += 1
        
        # Store filtering state and custom P-values
        if tool in ['rms', 'rnamotifscan']:
            self.user_rms_filtering_enabled = filtering_enabled
            self.user_rms_custom_pvalues = custom_pvalues
        elif tool in ['rmsx', 'rnamotifscanx']:
            self.user_rmsx_filtering_enabled = filtering_enabled
            self.user_rmsx_custom_pvalues = custom_pvalues
            if source_id == 7:
                self.rmsx_pipeline_config = self._build_internal_rmsx_config()
        elif tool in ['nobias']:
            self.user_nobias_filtering_enabled = filtering_enabled
            self.user_nobias_custom_pvalues = custom_pvalues
        
        # Build status message
        status_msg = f"Source {source_id}: {source_info['name']}"
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
        
        self.logger.info("\nNext steps:")
        if self.loaded_pdb_id:
            self.logger.info(f"  rmv_load_motif             Fetch motif data for {self.loaded_pdb_id}")
        else:
            self.logger.info("  rmv_fetch <PDB_ID>         Load PDB structure")
            self.logger.info("  rmv_load_motif             Fetch motif data")
    
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
            self.logger.error("Usage: rmv_source info [<ID>] or rmv_sources")
    
    def _print_active_source_info(self):
        """Print info about the currently active source only."""
        if self.current_source_id is None or self.current_source_mode is None:
            self.logger.info("No source is currently active.")
            self.logger.info("")
            self.logger.info("Select a source first:")
            self.logger.info("  rmv_sources      List all available sources")
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
                # FR3D - integrated runtime info
                report = self._run_fr3d_runtime_setup(build=False)
                print("\n--- Pipeline ---")
                print("  Extracts Hairpin Loops (HL) and Internal Loops (IL) from any PDB")
                print("  using the RNA 3D Motif Atlas flankSS algorithm (Petrov et al. 2013).")
                print("  Uses the bundled integrated FR3D runtime for local extraction.")
                print("\n--- Execution order (rmv_load_motif) ---")
                print("  1. Ensure bundled FR3D runtime is present")
                print("  2. Run fr3d_loop_extractor.py fresh for the active PDB")
                print("  3. Write {PDB}_fr3d_loops.csv and load motifs into RSMViewer")
                print("\n--- Current configuration ---")
                print("  Mode         : integrated runtime (no external config)")
                print(f"  Runtime mode : {report.get('runtime_mode', 'not_ready')}")
                print(f"  FR3D dir     : {report.get('fr3d_python_dir', '(auto-detect / not found)')}")
                print(f"  Python       : {report.get('python_executable', '(auto-detect / not found)')}")
                print(f"  Output dir   : {self.fr3d_output_path}")
                print("\n--- One-time setup ---")
                print("  1. Install normal Python dependencies required by the plugin")
                print("  2. Ensure the bundled FR3D runtime is shipped with the plugin package")
                print("  3. In PyMOL: rmv_fr3d doctor")
            elif source_id == 7:
                pipe_cfg = getattr(self, 'rmsx_pipeline_config', {})
                print("\n--- Pipeline ---")
                print("  Runs RNAMotifScanX per motif family, parses result_*.log files,")
                print("  and loads source 7 motifs directly into RSMViewer.")
                print("\n--- Current configuration ---")
                print("  Mode         : integrated runtime (no external config)")
                print(f"  Runtime dir  : {self.rmsx_runtime_dir}")
                print(f"  Executable   : {pipe_cfg.get('rmsx_executable', '(not found)')}")
                print(f"  Output dir   : {self.rmsx_output_path}")
                print(f"  Families     : {', '.join(pipe_cfg.get('motif_families', [])) or 'all defaults'}")
                print("\n--- Suggested workflow ---")
                print("  1. rmv_db 7")
                print("  2. rmv_fetch <PDB_ID>")
                print("  3. rmv_rmsx run_current     (or: rmv_rmsx run <PDB_ID>)")
                print("     rmv_rmsx run <PDB_ID>        # Full fresh rerun")
                print("  4. rmv_summary / rmv_show")
            else:
                print(f"\nAvailable motif types will be shown after loading a structure")
                print(f"with rmv_fetch <PDB_ID>")
        
        # Display sample commands
        print(f"\n--- Sample commands ---")
        if source_id == 5:
            print(f"  rmv_db 5                               Select integrated FR3D source")
            print(f"  rmv_fetch 1S72                          Load PDB structure")
            print(f"  rmv_load_motif                          Run pipeline if needed, then load")
            print(f"  rmv_summary                             Show HAIRPIN LOOP / INTERNAL LOOP counts")
            print(f"  rmv_show HAIRPIN LOOP                   Render all hairpin loops")
            print(f"  rmv_show INTERNAL LOOP 1                Zoom to internal loop instance 1")
            print(f"  rmv_fr3d run 1S72                       Force re-run for 1S72")
            print(f"  rmv_fr3d doctor                         Show runtime diagnostics")
            print(f"\n--- Combine FR3D with other sources ---")
            print(f"  rmv_db 5 7                              FR3D + RMSX combined")
            print(f"  rmv_db 5 3                              FR3D + BGSU API combined")
        else:
            print(f"  rmv_db {source_id}                 Select this source")
            print(f"  rmv_fetch 1S72                  Load structure")
            print(f"  rmv_load_motif                  Fetch motif data")
            print(f"  rmv_summary                     Show available motifs")
            print(f"  rmv_show HL                     Render motif")
        
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
        print("  rmv_sources                    List all sources (this display)")
        print("="*70 + "\n")
    
    def _handle_user_source(self, tool_name):
        """Handle user annotations source selection."""
        if not tool_name:
            self.logger.error("Usage: rmv_db user <tool_name> [on|off]")
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
        
        # Print next steps
        self.logger.info("")
        self.logger.info("Next steps:")
        if self.loaded_pdb_id:
            self.logger.info(f"  rmv_load_motif             Fetch motif data for {self.loaded_pdb_id}")
        else:
            self.logger.info("  rmv_fetch <PDB_ID>         Load PDB structure")
            self.logger.info("  rmv_load_motif             Fetch motif data")
        self.logger.info("")
    
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
            self.logger.error("Usage: rmv_db combine <ID1> <ID2> [<ID3> ...]")
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
        
        Uses the last loaded PDB and last selected source (or combined sources).
        Clears the cached data for that PDB, then re-fetches fresh motif data
        from the same source(s) that were last used.
        
        Args:
            pdb_id (str): PDB ID to refresh (uses currently loaded PDB if not specified)
        """
        try:
            # Determine PDB ID - use current if not specified
            if not pdb_id:
                pdb_id = self.loaded_pdb_id
            
            if not pdb_id:
                self.logger.error("No structure loaded. Use rmv_fetch <PDB_ID> first.")
                return
            
            pdb_id = pdb_id.upper()
            
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
            
            self.logger.info(f"Clearing cache and re-collecting motifs for {pdb_id} from {source_desc}...")
            
            # Clear cache for this PDB
            from .database import get_source_selector
            source_selector = get_source_selector()
            
            if source_selector and hasattr(source_selector, '_cache_manager'):
                try:
                    source_selector._cache_manager.clear_cache_for_pdb(pdb_id)
                    self.logger.debug(f"Cleared cache entries for {pdb_id}")
                except Exception:
                    pass  # Cache clearing is best-effort
            
            # Re-run the same fetch pipeline that rmv_load_motif uses
            self.fetch_motif_data_action(pdb_id)
            
            self.logger.success(f"Refresh complete for {pdb_id} from {source_desc}")
            self.logger.info(f"Next: rmv_summary | rmv_show <TYPE>")
                
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
            print(f"  Run: rmv_sources   (check available sources)")
            print(f"  Run: rmv_db <N>    (1-8)")
        
        # Always show workflow steps
        print("\n" + "-"*70)
        print("   WORKFLOW:")
        print("     Step 1: rmv_fetch <PDB_ID>       # Load PDB structure")
        print("     Step 2: rmv_sources               # Check available sources")
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
    def fetch_raw_pdb(pdb_id='', background_color='', cif_use_auth=''):
        """PyMOL command: Load raw PDB structure only (no motif data).
        
        Downloads and loads the PDB/mmCIF structure into PyMOL.
        Use rmv_db; rmv_load_motif after to select source and fetch motif data.
        Previously loaded PDB+source motif datasets are preserved for cross-PDB
        superimposition (rmv_super / rmv_align) until rmv_reset is called.
        
        Usage:
            rmv_fetch 1S72                           # Load PDB structure
            rmv_fetch 1S72, bg_color=lightgray       # With background color
            rmv_fetch 1S72, cif_use_auth=0           # Use label_asym_id chains
        
        Chain ID modes:
            cif_use_auth=1 (default)  - Use auth_asym_id (0, 9, A, B...)
            cif_use_auth=0            - Use label_asym_id (AA, BA, CA...)
        """
        if not pdb_id:
            gui.logger.error("Usage: rmv_fetch <PDB_ID> [, bg_color=gray80] [, cif_use_auth=0]")
            gui.logger.error("Examples:")
            gui.logger.error("  rmv_fetch 1S72")
            gui.logger.error("  rmv_fetch 1S72, bg_color=lightgray")
            gui.logger.error("  rmv_fetch 1S72, cif_use_auth=0    (use label_asym_id)")
            return
        
        pdb_arg = str(pdb_id).strip()
        bg_arg = str(background_color).strip() if background_color else None
        
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
        if cif_use_auth:
            cif_str = str(cif_use_auth).strip()
            if cif_str in ('0', 'off', 'false', 'label'):
                cif_auth_val = 0
            elif cif_str in ('1', 'on', 'true', 'auth'):
                cif_auth_val = 1
        
        # Detect whether the argument is a local file path or a PDB ID
        import os
        is_file_path = (
            os.sep in pdb_arg or
            pdb_arg.startswith('~') or
            pdb_arg.startswith('.') or
            pdb_arg.endswith(('.pdb', '.cif', '.mmcif', '.pdb.gz', '.cif.gz', '.ent', '.ent.gz'))
        )
        
        if is_file_path:
            expanded = os.path.expanduser(pdb_arg)
            if not os.path.isfile(expanded):
                gui.logger.error(f"File not found: {expanded}")
                return
            
            # Derive a structure name from the filename (strip extension)
            base = os.path.basename(expanded)
            structure_name = base.split('.')[0].lower()
            display_id = structure_name.upper()
        else:
            # Validate PDB ID
            if not pdb_arg or len(pdb_arg) != 4 or not pdb_arg.isalnum():
                gui.logger.error(f"Invalid PDB ID: '{pdb_arg}'")
                gui.logger.error("PDB ID must be exactly 4 alphanumeric characters (e.g., 1S72)")
                gui.logger.error("For local files: rmv_fetch /path/to/file.cif")
                return
            expanded = None
            structure_name = pdb_arg.lower()
            display_id = pdb_arg.upper()
        
        gui.cif_use_auth = cif_auth_val
        
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
                cmd.fetch(pdb_arg, structure_name)

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
            gui.logger.success(f"Loaded structure {display_id} as '{structure_name}'")
            gui.logger.info(f"Chain ID mode: {chain_mode}")
            
            # Get and display chains
            try:
                chains = cmd.get_chains(structure_name)
                if chains:
                    gui.logger.info(f"Chains found: {', '.join(chains[:20])}" +
                                   (f" ... (+{len(chains)-20} more)" if len(chains) > 20 else ""))
            except:
                pass
            
            gui.logger.info("")
            gui.logger.info("Next steps:")
            if gui.current_source_mode:
                gui.logger.info("  rmv_load_motif             Fetch motif data from current source")
                gui.logger.info("  rmv_view all               Highlight all motifs on structure")
            else:
                gui.logger.info("  rmv_db <N>                 Select a motif data source (1-8)")
            gui.logger.info("  rmv_sources                List all available sources")
            gui.logger.info("")

            # Optional automation: run FR3D structural motif download right after fetch.
            if gui.fr3d_auto_run_on_fetch:
                gui.logger.info("FR3D auto-run is ON: fetching FR3D structural motifs...")
                ok = gui.run_fr3d_wrapper(display_id)
                if ok:
                    gui.logger.success("FR3D motifs loaded automatically after rmv_fetch")
                else:
                    gui.logger.warning("FR3D auto-run after rmv_fetch failed; structure is still loaded")

            if gui.rmsx_auto_run_on_fetch:
                gui.logger.info("RNAMotifScanX auto-run is ON: running RNAMotifScanX wrapper...")
                ok = gui.run_rmsx_wrapper(display_id)
                if ok:
                    gui.logger.success("RNAMotifScanX motifs loaded automatically after rmv_fetch")
                else:
                    gui.logger.warning("RNAMotifScanX auto-run after rmv_fetch failed; structure is still loaded")
            
        except Exception as e:
            gui.logger.error(f"Failed to load {pdb_arg}: {str(e)}")
    
    def load_motif_data(argument=''):
        """PyMOL command: Fetch motif data from the selected source for the loaded PDB.
        
        Requires:
            1. A PDB structure must be loaded first (rmv_fetch)
            2. A source must be selected (rmv_db)
        
        Usage:
            rmv_load_motif               Fetch motifs from current source
        
        Workflow:
            rmv_fetch 1S72               # Step 1: Load PDB structure
            rmv_sources                  # Step 2: Check available sources
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
            gui.logger.info("  rmv_sources                (list all)")
            return
        
        pdb_id = gui.loaded_pdb_id
        
        # Dispatch to appropriate loader based on source mode
        if gui.current_source_mode == 'user' and gui.current_user_tool:
            gui.load_user_annotations_action(gui.current_user_tool, pdb_id)
        else:
            gui.fetch_motif_data_action(pdb_id, None)
    
    def load_structure(pdb_id_or_path='', background_color='', database=''):
        """PyMOL command: Load structure and automatically show all motifs.
        
        NOTE: This command is deprecated in the recommended workflow.
        Users should use rmv_fetch first, then rmv_db; rmv_load_motif.
        
        Usage:
            rmv_load <pdb_id_or_path>
        """
        if not pdb_id_or_path:
            gui.logger.error("Usage: rmv_load <PDB_ID>")
            return
        
        pdb_arg = str(pdb_id_or_path).strip().upper()
        
        # Instead of trying to fetch + load motifs all at once (which hangs PyMOL),
        # guide the user through the proper step-by-step workflow.
        print("\n" + "=" * 60)
        print("  RECOMMENDED WORKFLOW")
        print("=" * 60)
        print(f"\n  To visualize motifs for {pdb_arg}, follow these steps:\n")
        print(f"  Step 1:  rmv_fetch {pdb_arg}        # Fetch the PDB structure")
        print(f"  Step 2:  rmv_sources               # Check available data sources")
        print(f"  Step 3:  rmv_db <N>                # Select data source (1-8)")
        print(f"  Step 4:  rmv_load_motif             # Fetch motif data")
        print(f"  Step 5:  rmv_summary               # Show motif types & counts")
        print(f"  Step 6:  rmv_summary <TYPE>         # Show instances of a type")
        print(f"  Step 7:  rmv_show <TYPE>            # Render a motif type")
        print(f"  Step 8:  rmv_show <TYPE> <NO>       # Zoom to specific instance")
        print(f"\n  Type 'rmv_sources' to see all available data sources.")
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
                gui.logger.error(f"Usage: rmv_toggle MOTIF_TYPE on/off")
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
    
    def show_help():
        """PyMOL command: Show all available commands."""
        gui.print_help()
    
    def set_bg_color(color_name='gray80'):
        """PyMOL command: Change background color of non-motif residues."""
        color_arg = str(color_name).strip()
        if not color_arg:
            color_arg = 'gray80'
        gui.set_background_color(color_arg)
    
    def motif_summary(motif_type='', instance_no=''):
        """PyMOL command: Show motif summary table (console only, no rendering).
        
        Usage:
            rmv_summary              Show all motifs summary for loaded PDB
            rmv_summary HL           Show detailed instances of HL motif
            rmv_summary HL 1         Show specific HL instance #1
        """
        if not motif_type:
            # Show general motif summary
            gui.print_motif_summary()
        else:
            # Check if instance number is provided
            motif_arg = str(motif_type).strip().upper()
            
            # Handle both formats: "HL 1" and separate args
            if instance_no:
                try:
                    inst_no = int(instance_no)
                    gui.show_motif_instance_summary(motif_arg, inst_no)
                except ValueError:
                    gui.logger.error("Instance number must be an integer")
            else:
                # Check if the motif_type contains instance number at the end
                parts = motif_arg.split()
                if len(parts) >= 2 and parts[-1].isdigit():
                    # Last part is a number - treat as instance number
                    motif_name = ' '.join(parts[:-1])
                    inst_no = int(parts[-1])
                    gui.show_motif_instance_summary(motif_name, inst_no)
                else:
                    # Show all instances of the motif type
                    gui.show_motif_summary_for_type(motif_arg)
    
    def select_database(mode='', tool='', jaccard_threshold=''):
        """PyMOL command: Select motif data source by ID number.
        
        Usage:
            rmv_db 1                  - RNA 3D Motif Atlas (Local)
            rmv_db 2                  - Rfam (Local)
            rmv_db 3                  - BGSU RNA 3D Hub (Online)
            rmv_db 4                  - Rfam API (Online)
            rmv_db 5                  - FR3D Annotations (User)
            rmv_db 6                  - RNAMotifScan (RMS - User)
            rmv_db 7                  - RNAMotifScanX (RMSX - User)
        
        Multi-Source Combine (with enrichment + cascade merge):
            rmv_db 1 3               - Combine Atlas + BGSU (Atlas = priority)
            rmv_db 2 5 3             - Combine 3 sources (Rfam = highest priority)
        
        With optional Jaccard threshold for cascade merge:
            rmv_db 8 6, jaccard_threshold=0.80   - Combine with 80% threshold
            rmv_db 8 6, jaccard_threshold=80      - Same (values >1 treated as percentage)
            rmv_db 1 3, jaccard_threshold=0.70    - Combine Atlas+BGSU at 70%
            (default: 0.60 i.e. 60% if not specified)
        
        With optional parameters:
            rmv_db 6 off                       - RMS with filtering OFF
            rmv_db 6 on                        - RMS with filtering ON
            rmv_db 6 C-LOOP 0.05 KINK-TURN 0.02 - RMS with custom P-values
            rmv_db 7 C-LOOP_CONSENSUS 0.01    - RMSX with custom P-value
        
        With custom data path (sources 5-8):
            rmv_db 5 /path/to/fr3d/data       - FR3D with custom data directory
            rmv_db 6 /path/to/rms/data        - RMS with custom data directory
            rmv_db 7 /path/to/rmsx/data       - RMSX with custom data directory
            rmv_db 8 /path/to/nobias/data     - NoBIAS with custom data directory
        """
        if not mode:
            print("\n  [rmv_db] Data source selection")
            print("\n  Usage:")
            print("    rmv_db <SOURCE_ID> [options]")
            print("\n  Common examples:")
            print("    rmv_db 3                              Select BGSU API")
            print("    rmv_db 5                              Select FR3D (strict local pipeline)")
            print("    rmv_db 7                              Select RMSX")
            print("    rmv_db 8 7                            Combine NoBIAS + RMSX")
            print("    rmv_db 8 7, jaccard_threshold=0.80    Combine with custom merge threshold")
            print("    rmv_db 7 off                          Disable RMSX P-value filtering")
            print("    rmv_db 7 K-TURN 0.02                  Set custom P-value for one motif")
            print("\n  Path override examples (sources 5-8):")
            print("    rmv_db 5 /path/to/fr3d/data")
            print("    rmv_db 7 /path/to/rmsx/data")
            print("\n  Stage status:")
            print(f"    Current source ID: {gui.current_source_id if gui.current_source_id is not None else '(not selected)'}")
            print(f"    Loaded structure : {gui.loaded_pdb_id if gui.loaded_pdb_id else '(none)'}")
            print("\n  Related:")
            print("    rmv_sources      List all source IDs")
            print("    rmv_source info  Show selected source details")
            return
        
        # --- Parse optional jaccard_threshold kwarg ---
        jt_str = str(jaccard_threshold).strip() if jaccard_threshold else ''
        if jt_str:
            # Strip trailing '%' if present
            jt_clean = jt_str.rstrip('%')
            try:
                jt_val = float(jt_clean)
                # Values > 1 are treated as percentages (e.g., 80 -> 0.80)
                if jt_val > 1.0:
                    jt_val = jt_val / 100.0
                if not (0.0 < jt_val <= 1.0):
                    gui.logger.error(f"jaccard_threshold must be between 0 and 1 (or 1-100 as %). Got: {jt_str}")
                    return
                gui.jaccard_threshold = jt_val
                gui.logger.info(f"Jaccard similarity threshold set to {gui.jaccard_threshold:.0%}")
            except ValueError:
                gui.logger.error(f"Invalid jaccard_threshold value: '{jt_str}'. Use a number like 0.80 or 80")
                return
        
        mode_arg = str(mode).strip()
        tool_arg = str(tool).strip() if tool else None
        
        # Handle PyMOL passing arguments as combined string
        parts = mode_arg.split(None, 1)  # Split on first whitespace
        first_part = parts[0].lower()
        remaining_arg = parts[1] if len(parts) > 1 else tool_arg
        
        # Check if first argument is a number (source ID)
        try:
            source_id = int(first_part)
            gui._handle_source_by_id(source_id, remaining_arg)
            return
        except ValueError:
            pass  # Not a number, check for other commands
        
        # Legacy support for old syntax (backward compat during transition)
        if first_part == 'combine':
            gui._handle_combine_sources(remaining_arg)
        elif first_part == 'user':
            gui._handle_user_source(remaining_arg)
        elif first_part == 'local':
            gui._handle_local_source(remaining_arg)
        elif first_part == 'web':
            gui._handle_web_source(remaining_arg)
        else:
            # Old-style: auto, all, etc.
            gui.set_source_mode(first_part)
    
    def set_source(mode='', source_id=''):
        """PyMOL command: Show current source info or detailed info about a specific source.
        
        Usage:
            rmv_source info          - Show currently selected source info
            rmv_source info <N>      - Show detailed info about source N (1-8)
        """
        if not mode:
            gui.logger.error("Usage: rmv_source info [<ID>]")
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
        
        gui.logger.error(f"Unknown subcommand: {first_part}")
        gui.logger.error("Usage: rmv_source info [<ID>]")
        gui.logger.error("       rmv_db <ID>        - Select a source")
        gui.logger.error("  rmv_source info <N>  - N can be 1-8")
    
    def refresh_motifs(pdb_id=''):
        """PyMOL command: Force refresh cache and collect motif data again.
        
        Clears cached data for the currently loaded PDB and re-fetches
        motif information from the last selected source (or combined
        sources if combine mode was used).
        
        Usage:
            rmv_refresh        - Refresh current PDB from last selected source
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
        
        # Separate trailing numeric parts (instance numbers) from motif name.
        # Walk from the end: pure-digit parts are instance numbers.
        instance_nums = []
        while all_parts and all_parts[-1].isdigit():
            instance_nums.insert(0, int(all_parts.pop()))
        
        if not all_parts:
            gui.logger.error("Usage: rmv_show <MOTIF_TYPE> [<INSTANCE_NO>]")
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
            rmv_view all              Highlight ALL motif regions on structure
            rmv_view K-TURN           Zoom to all K-TURN instances
            rmv_view K-TURN 1         Zoom to instance #1 and create selection
            rmv_view hide             Reset all view coloring to gray
            rmv_view K-TURN hide      Reset only K-TURN view coloring
        """
        # Build parts list - split first arg by spaces for robustness
        # (PyMOL may pass 'K-TURN hide' as a single string)
        first_parts = str(motif_type).split() if motif_type else []
        rest_parts = [str(a) for a in extra_args]
        all_parts = first_parts + rest_parts
        all_parts = [p.strip() for p in all_parts if p.strip()]

        if not all_parts:
            print("\n  [rmv_view] Highlight motif regions (no motif objects)")
            print("\n  Usage:")
            print("    rmv_view all")
            print("    rmv_view <MOTIF_TYPE>")
            print("    rmv_view <MOTIF_TYPE> <INSTANCE_NO>")
            print("    rmv_view hide")
            print("    rmv_view <MOTIF_TYPE> hide")
            print("\n  Examples:")
            print("    rmv_view all")
            print("    rmv_view K-TURN")
            print("    rmv_view K-TURN 1")
            print("    rmv_view K-TURN hide")
            print("\n  Stage checks:")
            if not gui.loaded_pdb_id:
                print("    ERROR: No structure loaded. Run: rmv_fetch <PDB_ID>")
            loaded = gui.viz_manager.motif_loader.get_loaded_motifs() if gui.viz_manager and gui.viz_manager.motif_loader else {}
            if not loaded:
                print("    ERROR: No motifs loaded. Run: rmv_load_motif")
            else:
                print(f"    OK: {len(loaded)} motif types available for view operations")
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
                    filter_pdb=fpdb, filter_suffix=fsuf)
        else:
            gui.viz_manager.view_motif_type(
                motif_arg,
                filter_pdb=fpdb, filter_suffix=fsuf)

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
            print("\nFR3D wrapper commands (strict local FR3D pipeline, no fallback):")
            print("  rmv_fr3d status")
            print("  rmv_fr3d doctor")
            print("  rmv_fr3d setup")
            print("  rmv_fr3d run 1S72            Download hairpin/internal loop/junction motifs")
            print("  rmv_fr3d refresh [PDB_ID]    Force rerun and replace cached file")
            print("  rmv_fr3d run_current         Run for the currently loaded structure")
            print("  rmv_fr3d config /path/to/fr3d-python /abs/output on   (advanced)")
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
            gui.logger.error("Please specify PDB ID")
            print(f"  Usage: rmv_user {tool_arg} <PDB_ID>")
            return
        
        gui.load_user_annotations_action(tool_arg, pdb_arg)

    def fr3d_wrapper(action='', arg1='', *extra_args, **_kwargs):
        """PyMOL command: Run Source-5 FR3D local pipeline.

        Usage:
            rmv_fr3d status
            rmv_fr3d doctor
            rmv_fr3d setup
            rmv_fr3d refresh [PDB_ID]
            rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]
            rmv_fr3d run <PDB_ID>
            rmv_fr3d run_current

        'run' executes local FR3D extraction for Hairpin Loops (HL), Internal
        Loops (IL), and Junctions and loads them as structural RNA motifs in
        RSMViewer (source 5).
        """
        action_arg = str(action).strip() if action else ''
        arg1_str = str(arg1).strip() if arg1 else ''

        # Handle combined-string invocation from PyMOL
        if action_arg and not arg1_str:
            parts = action_arg.split()
            if len(parts) > 1:
                action_arg = parts[0]
                arg1_str = parts[1]

        sub = action_arg.lower() if action_arg else 'status'

        if sub in ['', 'status', 'show']:
            gui.print_fr3d_wrapper_status()
            return

        if sub == 'doctor':
            gui.fr3d_doctor(auto_setup=False)
            return

        if sub == 'setup':
            gui.fr3d_doctor(auto_setup=True)
            return

        if sub == 'refresh':
            target_pdb = arg1_str or gui.loaded_pdb_id
            if not target_pdb:
                gui.logger.error("No active structure. Use rmv_fetch <PDB_ID> first or pass PDB ID.")
                return
            gui.run_fr3d_wrapper(target_pdb, force_refresh=True)
            return

        if sub == 'config':
            if not arg1_str:
                gui.logger.error("Usage: rmv_fr3d config <FR3D_ROOT> [OUTPUT_DIR] [AUTO_ON_FETCH]")
                gui.logger.info("Example: rmv_fr3d config /path/to/fr3d-python /abs/fr3d_out on")
                return

            extras = [str(x).strip() for x in extra_args if str(x).strip()]
            output_dir = extras[0] if len(extras) >= 1 else ''
            auto_on_fetch = extras[1] if len(extras) >= 2 else ''
            gui.configure_fr3d_wrapper(arg1_str, output_dir, auto_on_fetch)
            return

        if sub == 'run':
            if not arg1_str:
                gui.logger.error("Usage: rmv_fr3d run <PDB_ID>")
                return
            gui.run_fr3d_wrapper(arg1_str)
            return

        if sub in ['run_current', 'current']:
            if not gui.loaded_pdb_id:
                gui.logger.error("No active structure. Use rmv_fetch <PDB_ID> first.")
                return
            gui.run_fr3d_wrapper(gui.loaded_pdb_id)
            return

        gui.logger.error(f"Unknown rmv_fr3d subcommand: {sub}")
        gui.logger.info("Use: rmv_fr3d status | doctor | setup | refresh [PDB_ID] | config | run <PDB_ID> | run_current")

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
                gui.logger.error("Usage: rmv_rmsx config <EXECUTABLE> [OUTPUT_DIR] [WORK_DIR] [AUTO_ON_FETCH] [QUERY_FILE]")
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
                gui.logger.error("Usage: rmv_rmsx run <PDB_ID> [EXTRA_ARGS]")
                return
            extras = ' '.join(str(x).strip() for x in extra_args if str(x).strip())
            gui.run_rmsx_wrapper(arg1_str, extras, force_fresh=True)
            return

        if sub in ['run_current', 'current']:
            if not gui.loaded_pdb_id:
                gui.logger.error("No active structure. Use rmv_fetch <PDB_ID> first.")
                return
            extras = [arg1_str] if arg1_str else []
            extras.extend(str(x).strip() for x in extra_args if str(x).strip())
            gui.run_rmsx_wrapper(gui.loaded_pdb_id, ' '.join(extras).strip(), force_fresh=False)
            return

        gui.logger.error(f"Unknown rmv_rmsx subcommand: {sub}")
        gui.logger.info("Use: rmv_rmsx status | config | args | doctor | setup | test | run | run_current")

    def rmsx_doctor_cmd(*_args, **_kwargs):
        """PyMOL command: Show integrated RMSX runtime diagnostics."""
        gui.rmsx_doctor(auto_setup=False)
    
    # Add commands to PyMOL
    cmd.extend('rmv_fetch', fetch_raw_pdb)
    cmd.extend('rmv_load_motif', load_motif_data)
    cmd.extend('rmv_load', load_structure)
    cmd.extend('rmv_toggle', toggle_motif)
    cmd.extend('rmv_sources', list_sources)
    cmd.extend('rmv_help', show_help)
    cmd.extend('rmv_bg_color', set_bg_color)
    cmd.extend('rmv_summary', motif_summary)
    cmd.extend('rmv_db', select_database)
    cmd.extend('rmv_source', set_source)
    cmd.extend('rmv_refresh', refresh_motifs)
    cmd.extend('rmv_show', show_motif)
    cmd.extend('rmv_view', view_motif)
    cmd.extend('rmv_user', load_user_annotations)
    cmd.extend('rmv_fr3d', fr3d_wrapper)
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
            gui.logger.error("Please specify a color")
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
    
    def save_motif_images(argument=''):
        """PyMOL command: Save motif instance images or extract mmCIF structures.
        
        Usage (images):
            rmv_save ALL                      Save all motif types and instances (default: cartoon)
            rmv_save ALL sticks               Save all motifs as sticks representation
            rmv_save HL                       Save all hairpin loop instances (default: cartoon)
            rmv_save HL sticks                Save all HL instances as sticks
            rmv_save HL 3                     Save 3rd HL instance (default: cartoon)
            rmv_save HL 3 spheres             Save 3rd HL instance as spheres
            rmv_save current                  Save current PyMOL view
            rmv_save current my_view.png      Save current view to file
        
        Usage (mmCIF structure export):
            rmv_save ALL cif                  Export all motif instances as mmCIF
            rmv_save HL cif                   Export all HL instances as mmCIF
            rmv_save HL 3 cif                 Export 3rd HL instance as mmCIF
        
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
            Images:     plugin_dir/motif_images/pdb_id/MOTIF_TYPE/<type>-<no>-<chain>-<residues>.png
            Structures: plugin_dir/motif_structures/pdb_id/MOTIF_TYPE/<type>-<no>-<chain>-<residues>.cif
        """
        arguments = str(argument).strip().split()
        
        if not arguments:
            print("\nUsage: rmv_save <ALL | MOTIF_TYPE | MOTIF_TYPE INSTANCE_ID | current> [representation | cif]")
            print("\n  IMAGE SAVE EXAMPLES:")
            print("  rmv_save ALL             Save all motif images (cartoon)")
            print("  rmv_save ALL sticks      Save all motif images (sticks)")
            print("  rmv_save HL              Save all hairpin loop images (cartoon)")
            print("  rmv_save HL sticks       Save all HL images (sticks)")
            print("  rmv_save HL 1            Save specific HL instance #1 (cartoon)")
            print("  rmv_save HL 1 spheres    Save specific HL instance #1 (spheres)")
            print("  rmv_save current         Save current PyMOL view")
            print("  rmv_save current out.png Save current view to file")
            print("\n  mmCIF STRUCTURE EXPORT:")
            print("  rmv_save ALL cif         Export ALL motif structures as mmCIF")
            print("  rmv_save HL cif          Export all HL instances as mmCIF")
            print("  rmv_save HL 3 cif        Export HL instance #3 as mmCIF")
            print("\n  Note: mmCIF export uses ORIGINAL coordinates from the on-disk CIF,")
            print("        not PyMOL's internal coordinates.")
            print("        Output is coordinates-only (_atom_site loop for motif residues).")
            print("\nRepresentations: cartoon, sticks, spheres, ribbon, lines, licorice, surface, cartoon+sticks")
            print("\nOutput goes to:")
            print("  Images:     plugin_dir/motif_images/pdb_id/MOTIF_TYPE/")
            print("  Structures: plugin_dir/motif_structures/pdb_id/MOTIF_TYPE/")
            return
        
        pdb_id = gui.viz_manager.structure_loader.get_current_pdb_id()
        if not pdb_id:
            gui.logger.error("No structure loaded")
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
                # Check for possible typos against known keywords and motif types
                all_candidates = ['ALL', 'CURRENT'] + sorted(loaded_motifs.keys())
                suggestion = _suggest(arguments[0], all_candidates)
                gui.logger.error(f"Unknown argument '{arguments[0]}'")
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
        gui.motif_visibility = {}
        gui.current_source_mode = None
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
        print("\n  All objects deleted and plugin state cleared.")
        print("  Ready for a fresh session.")
        print("\n  Quick Start:")
        print("     rmv_fetch <PDB_ID>       # Load a PDB structure")
        print("     rmv_sources               # Check available data sources")
        print("     rmv_db <N>                # Select data source (1-8)")
        print("     rmv_load_motif            # Fetch motif data")
        print("     rmv_summary               # Show motif types & counts")
        print("     rmv_view all              # Highlight all motifs on structure")
        print("     rmv_show HL               # Render hairpin loops")
        print()
    
    cmd.extend('rmv_reset', reset_plugin)
    
    # Register base-pair visualization commands (rmv_pair, rmv_pair_batch)
    from .pair_visualizer import register_pair_commands
    register_pair_commands()
    
    # Register medoid superimposition commands (rmv_super, rmv_align)
    from .alignment import register_alignment_commands
    register_alignment_commands()

    gui.logger.success("RSMViewer GUI initialized")
