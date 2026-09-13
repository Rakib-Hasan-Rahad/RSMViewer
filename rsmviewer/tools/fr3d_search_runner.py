#!/usr/bin/env python3
"""Minimal, non-invasive subprocess runner for official BGSU fr3d-python.

This script reproduces an official FR3D search WITHOUT editing, patching, or
vendoring the user's fr3d-python checkout. It works by:

  1. Putting BOTH the repo root and ``<repo>/fr3d/search`` on ``sys.path`` so
     that FR3D's mixed flat/package imports resolve.
  2. Importing ``fr3d_configuration`` and OVERRIDING its module attributes in
     ``sys.modules`` (SERVER=False and all paths) BEFORE importing FR3D. Because
     Python caches the module, every ``from fr3d_configuration import X`` inside
     FR3D binds the overridden values.
  3. Rewriting the query's ``searchFiles`` to point at a single staged target
     ``.cif`` (cif_local mode) so FR3D annotates the exact structure the user is
     viewing, using its own official pairwise annotator, with no network needed.

It never imports RSMViewer and can run under the user's own interpreter.

Contract (stdin JSON on argv[1] as a path, or --spec <path>):
    {
      "fr3d_root":      "/abs/path/to/fr3d-python",
      "query_file":     "/abs/path/to/query.json",   # a WebFR3D-style JSON query
      "target_cif":     "/abs/path/to/TARGET.cif",    # cif_local reference target
      "run_dir":        "/abs/path/to/run",           # all outputs go under here
      "query_name":     "sanitized_name",             # used for output naming
      "allow_network":  false,
    "data_mode":      "run_fr3d_pipeline"
    }

Emits a single JSON object on stdout:
    {"ok": true, "result": {"csv_path": ..., "candidate_count": N,
                             "query_name": ..., "executed_query": {...}}}
or
    {"ok": false, "error": "..."}
"""

import argparse
import types
import importlib.util
import json
import os
import pickle
import re
import sys
import time
import traceback


def _fail(message):
    print(json.dumps({"ok": False, "error": str(message)}))
    sys.exit(1)


def _autofix_empty_blocks(source, filename):
    """Insert a semantically-inert ``pass`` into compile-blocking empty suites.

    The official BGSU fr3d-python checkout contains at least one comment-only
    block (``if len(chains) == 0:`` in ``query_processing.py``) that makes the
    whole module fail to COMPILE under Python 3. Because Python compiles a
    module in full at import time, this blocks ``import FR3D`` even though the
    offending branch is never reached in cif_local mode.

    This function ONLY reacts to ``IndentationError: expected an indented
    block`` and ONLY ever inserts ``pass``. It never modifies existing
    statements and never touches any file on disk. Returns
    ``(fixed_source, applied_header_lines)``.
    """
    lines = source.split("\n")
    applied = []
    for _ in range(500):  # bounded safety
        try:
            compile("\n".join(lines), filename, "exec")
            return "\n".join(lines), applied
        except IndentationError as exc:
            msg = exc.msg or ""
            if "expected an indented block" not in msg:
                raise
            match = re.search(r"on line (\d+)", msg)
            if not match:
                raise
            header_ln = int(match.group(1))  # 1-based header line
            if header_ln < 1 or header_ln > len(lines):
                raise
            header = lines[header_ln - 1]
            indent = len(header) - len(header.lstrip(" "))
            lines.insert(header_ln, (" " * (indent + 4)) + "pass")
            applied.append(header_ln)
    return source, applied


def _preload_with_autofix(mod_name, file_path):
    """Load a source module into ``sys.modules``, auto-fixing empty blocks.

    Used only to get past compile-blocking stubs in the official checkout
    WITHOUT editing any file on disk. If no fix is needed, returns ``[]`` and
    leaves normal import to handle the module. Returns the list of fixed
    (1-based) header line numbers.
    """
    if not os.path.isfile(file_path):
        return []
    with open(file_path, "r", encoding="utf-8") as handle:
        source = handle.read()
    fixed_source, applied = _autofix_empty_blocks(source, file_path)
    if not applied:
        return []
    spec = importlib.util.spec_from_loader(mod_name, loader=None, origin=file_path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = file_path
    sys.modules[mod_name] = module
    exec(compile(fixed_source, file_path, "exec"), module.__dict__)
    return applied


def _load_spec():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", default="")
    parser.add_argument("spec_positional", nargs="?", default="")
    args = parser.parse_args()
    spec_path = args.spec or args.spec_positional
    if not spec_path or not os.path.isfile(spec_path):
        _fail("Runner spec file not found: %r" % spec_path)
    try:
        with open(spec_path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except Exception as exc:  # noqa: BLE001
        _fail("Failed to read runner spec: %s" % exc)


def _detect_fr3d_entrypoint(search_dir):
    """Detect which public FR3D entrypoint this checkout exposes."""
    fr3d_py = os.path.join(search_dir, "FR3D.py")
    try:
        with open(fr3d_py, "r", encoding="utf-8") as handle:
            source = handle.read()
    except Exception:
        # Conservative fallback for older branches.
        return "main"

    if re.search(r"^def\s+main\s*\(", source, flags=re.MULTILINE):
        return "main"
    if re.search(r"^def\s+fr3d_search_from_query_names\s*\(", source, flags=re.MULTILINE):
        return "from_query_names"
    return "main"


def _normalize_query_schema(query, entrypoint_mode):
    """Normalize query JSON fields for the detected FR3D entrypoint style."""
    changed = {}

    interaction_matrix = query.get("interactionMatrix")

    if entrypoint_mode == "main":
        # Legacy/main path expects lower-case numpositions and list-like matrix indexing.
        if "numpositions" not in query and "numPositions" in query:
            try:
                query["numpositions"] = int(query.get("numPositions"))
                changed["numpositions"] = {
                    "original": query.get("numPositions"),
                    "executed": query.get("numpositions"),
                }
            except Exception:
                pass

        if isinstance(interaction_matrix, dict):
            # Convert dict-form matrices back to 2D list for older FR3D branches.
            keys = [int(k) for k in interaction_matrix.keys() if str(k).isdigit()]
            n = (max(keys) + 1) if keys else int(query.get("numpositions", query.get("numPositions", 0)) or 0)
            converted = [["" for _ in range(n)] for _ in range(n)]
            for i_key, row in interaction_matrix.items():
                if not str(i_key).isdigit() or not isinstance(row, dict):
                    continue
                i = int(i_key)
                if i >= n:
                    continue
                for j_key, cell in row.items():
                    if not str(j_key).isdigit():
                        continue
                    j = int(j_key)
                    if j < n:
                        converted[i][j] = str(cell or "")
            query["interactionMatrix"] = converted
            changed["interactionMatrix_format"] = {
                "original": "dict",
                "executed": "list",
            }
    else:
        # FR3D latest path expects numPositions and dict-style matrix keys.
        if "numPositions" not in query and "numpositions" in query:
            try:
                query["numPositions"] = int(query.get("numpositions"))
                changed["numPositions"] = {
                    "original": query.get("numpositions"),
                    "executed": query.get("numPositions"),
                }
            except Exception:
                pass

        if isinstance(interaction_matrix, list):
            converted = {}
            for i, row in enumerate(interaction_matrix):
                converted[str(i)] = {}
                if isinstance(row, list):
                    for j, cell in enumerate(row):
                        converted[str(i)][str(j)] = str(cell or "")
                else:
                    converted[str(i)][str(i)] = str(row or "")
            query["interactionMatrix"] = converted
            changed["interactionMatrix_format"] = {
                "original": "list",
                "executed": "dict",
            }

    return changed


def _build_fr3d_configuration_module(module_name, search_dir, cif_dir, data_dir, raw_dir, json_dir):
    """Create a minimal FR3D configuration module for checkouts missing it."""
    module = types.ModuleType(module_name)
    module.SERVER = False
    module.CIFPATH = cif_dir
    module.DATAPATH = data_dir
    module.DATAPATHUNITS = os.path.join(data_dir, "units")
    module.DATAPATHPAIRS = os.path.join(data_dir, "pairs")
    module.OUTPUTPATH = raw_dir + os.sep
    module.JSONPATH = json_dir + os.sep
    module.TEMPLATEPATH = search_dir + os.sep
    module.JSLOCATION = search_dir + os.sep
    module.MAXTIME = float("inf")
    module.MAXCANDIDATES = 1000000
    module.MAXCANDIDATESHEATMAP = 300
    module.REFRESHTIME = float("inf")
    module.Q = {}
    return module


def main():
    spec = _load_spec()

    fr3d_root = os.path.abspath(os.path.expanduser(str(spec.get("fr3d_root", "")).strip()))
    query_file = os.path.abspath(os.path.expanduser(str(spec.get("query_file", "")).strip()))
    target_cif = os.path.abspath(os.path.expanduser(str(spec.get("target_cif", "")).strip()))
    run_dir = os.path.abspath(os.path.expanduser(str(spec.get("run_dir", "")).strip()))
    query_name = str(spec.get("query_name", "") or "").strip() or "query"
    allow_network = bool(spec.get("allow_network", False))
    data_mode = str(spec.get("data_mode", "run_fr3d_pipeline") or "run_fr3d_pipeline").strip()
    mode_aliases = {
        "run_fr3d_pipeline": "cif_local",
        "fr3d_local_data": "fr3d_native_data",
        "rna3dhub_web_interactions": "rna3dhub_interactions",
        "cif_local": "cif_local",
        "fr3d_native_data": "fr3d_native_data",
        "rna3dhub_interactions": "rna3dhub_interactions",
    }
    canonical_mode = mode_aliases.get(data_mode, data_mode)

    search_dir = os.path.join(fr3d_root, "fr3d", "search")
    entrypoint_mode = _detect_fr3d_entrypoint(search_dir)
    if not os.path.isfile(os.path.join(search_dir, "FR3D.py")):
        _fail("Not a valid fr3d-python checkout (missing fr3d/search/FR3D.py): %s" % fr3d_root)
    if not os.path.isfile(query_file):
        _fail("Query file not found: %s" % query_file)
    if canonical_mode == "cif_local" and not os.path.isfile(target_cif):
        _fail("Target CIF not found for run_fr3d_pipeline mode: %s" % target_cif)

    # --- run-dir layout -----------------------------------------------------
    raw_dir = os.path.join(run_dir, "raw")            # OUTPUTPATH (FR3D CSV/HTML)
    data_dir = os.path.join(run_dir, "data")          # DATAPATH (units/pairs pickles)
    json_dir = os.path.join(run_dir, "json")          # JSONPATH (executed query)
    cif_dir = os.path.dirname(target_cif) if target_cif else run_dir
    for path in (raw_dir, data_dir, json_dir,
                 os.path.join(data_dir, "units"), os.path.join(data_dir, "pairs")):
        os.makedirs(path, exist_ok=True)

    # --- import path: fr3d/search FIRST, then classifiers, then repo root ---
    # The official checkout uses "flat" sibling imports (e.g. file_reading,
    # query_processing live in fr3d/search; NA_pairwise_interactions,
    # NA_unit_annotation live in fr3d/classifiers). search_dir goes first so
    # modules that exist in both trees (e.g. discrepancy) resolve to the
    # search copy, matching how FR3D itself is normally run.
    classifiers_dir = os.path.join(fr3d_root, "fr3d", "classifiers")
    for path in (fr3d_root, classifiers_dir, search_dir):
        if path in sys.path:
            sys.path.remove(path)
        sys.path.insert(0, path)

    # --- import + override configuration BEFORE importing FR3D -------------
    # The official checkout loads its config under TWO module identities:
    #   * flat  "fr3d_configuration"            (query_processing, write_output, FR3D)
    #   * package "fr3d.search.fr3d_configuration" (file_reading)
    # Each consumer does `from ... import NAME`, binding the value at import
    # time, so we must override BOTH module objects up-front -- before FR3D or
    # any of its dependencies import them -- otherwise hard-coded server paths
    # like /var/www/html leak through.
    cfg_modules = []
    cfg_import_error = None
    try:
        import fr3d_configuration as cfg
        cfg_modules.append(cfg)
    except Exception as exc:  # noqa: BLE001
        cfg_import_error = exc

    try:
        import fr3d.search.fr3d_configuration as pkg_cfg
        if pkg_cfg not in cfg_modules:
            cfg_modules.append(pkg_cfg)
    except Exception:
        pkg_cfg = None

    if not cfg_modules:
        cfg = _build_fr3d_configuration_module(
            "fr3d_configuration", search_dir, cif_dir, data_dir, raw_dir, json_dir
        )
        pkg_cfg = _build_fr3d_configuration_module(
            "fr3d.search.fr3d_configuration", search_dir, cif_dir, data_dir, raw_dir, json_dir
        )
        sys.modules["fr3d_configuration"] = cfg
        sys.modules["fr3d.search.fr3d_configuration"] = pkg_cfg
        cfg_modules.extend([cfg, pkg_cfg])
    elif cfg_import_error:
        # If one identity imports and the other does not, mirror the loaded
        # module under both names so FR3D's mixed imports stay consistent.
        fallback = cfg_modules[0]
        sys.modules.setdefault("fr3d_configuration", fallback)
        sys.modules.setdefault("fr3d.search.fr3d_configuration", fallback)

    overrides = {
        "SERVER": False,
        "CIFPATH": cif_dir,
        "DATAPATH": data_dir,
        "DATAPATHUNITS": os.path.join(data_dir, "units"),
        "DATAPATHPAIRS": os.path.join(data_dir, "pairs"),
        "OUTPUTPATH": raw_dir + os.sep,
        "JSONPATH": json_dir + os.sep,
        "TEMPLATEPATH": search_dir + os.sep,
        "JSLOCATION": search_dir + os.sep,
        "MAXTIME": float("inf"),
        "MAXCANDIDATES": 1000000,
        "MAXCANDIDATESHEATMAP": 300,
        "REFRESHTIME": float("inf"),  # never emit intermediate HTML refreshes
    }
    for module in cfg_modules:
        for key, value in overrides.items():
            setattr(module, key, value)

    # When network is disallowed, pre-stage an empty NA_datafile so FR3D's
    # readPDBDatafile() sees a fresh cache and never attempts a download.
    if not allow_network:
        try:
            stamp_path = os.path.join(data_dir, "units", "NA_datafile.pickle")
            with open(stamp_path, "wb") as handle:
                pickle.dump({}, handle)
        except Exception:
            pass

    # --- build the executed query: force cif_local target ------------------
    try:
        with open(query_file, "r", encoding="utf-8") as handle:
            query = json.load(handle)
    except Exception as exc:  # noqa: BLE001
        _fail("Failed to parse query JSON: %s" % exc)

    original_query = json.loads(json.dumps(query))  # deep copy for provenance
    changed_fields = {}

    schema_changes = _normalize_query_schema(query, entrypoint_mode)
    if schema_changes:
        changed_fields.update(schema_changes)

    # Deterministic output name so we can locate the CSV afterward.
    if str(query.get("name", "")).strip() != query_name:
        changed_fields["name"] = {"original": query.get("name"), "executed": query_name}
    query["name"] = query_name

    if canonical_mode == "cif_local":
        changed_fields["searchFiles"] = {"original": query.get("searchFiles"), "executed": [target_cif]}
        query["searchFiles"] = [target_cif]
        # Seed an empty PDB metadata map. In cif_local mode FR3D reads the
        # structure straight from the .cif via processPDBFile and does not use
        # precomputed PDB annotations. The official calculateQueryConstraints
        # only initialises Q["PDB_data_file"] on the 4-char-PDB-id code path,
        # so providing an empty map here keeps the local-file branch from
        # raising KeyError('PDB_data_file') and avoids any network fetch.
        if not isinstance(query.get("PDB_data_file"), dict):
            changed_fields["PDB_data_file"] = {"original": query.get("PDB_data_file"), "executed": {}}
            query["PDB_data_file"] = {}

    executed_name = "executed_query.json"
    executed_path = os.path.join(json_dir, executed_name)
    try:
        with open(executed_path, "w", encoding="utf-8") as handle:
            json.dump(query, handle, indent=2)
        with open(os.path.join(run_dir, "original_query.json"), "w", encoding="utf-8") as handle:
            json.dump(original_query, handle, indent=2)
    except Exception as exc:  # noqa: BLE001
        _fail("Failed to stage executed query: %s" % exc)

    # --- run official FR3D --------------------------------------------------
    started = time.time()
    fr3d_error = ""
    shim_fixes = []
    # Pre-load compile-blocking modules from source with a `pass`-only autofix,
    # so `import FR3D` succeeds without editing the user's checkout on disk.
    try:
        shim_fixes = _preload_with_autofix(
            "query_processing", os.path.join(search_dir, "query_processing.py")
        )
    except Exception:  # noqa: BLE001
        shim_fixes = []  # fall back to normal import; the real error will surface

    # In-memory compatibility bindings for upstream scoping gaps. file_reading's
    # writeNAUnitAnnotations() calls annotate_bond_orientation(), but the name is
    # only imported locally inside processPDBFile(), so it is missing from the
    # module namespace and raises NameError on the local-cif annotation path.
    # We bind the EXACT upstream function into file_reading -- no behavioural
    # change, no result fabrication, no edit to the checkout on disk.
    shim_bindings = []
    try:
        import file_reading as _fr
        if not hasattr(_fr, "annotate_bond_orientation"):
            from fr3d.classifiers.NA_unit_annotation import annotate_bond_orientation as _abo
            _fr.annotate_bond_orientation = _abo
            shim_bindings.append("file_reading.annotate_bond_orientation")
    except Exception:  # noqa: BLE001
        pass  # let the real error surface if the binding cannot be made

    try:
        import FR3D  # noqa: N813  (official module name)
        if hasattr(FR3D, "main") and callable(FR3D.main):
            FR3D.main([executed_name])
        elif hasattr(FR3D, "fr3d_search_from_query_names") and callable(FR3D.fr3d_search_from_query_names):
            FR3D.fr3d_search_from_query_names([executed_name])
        else:
            raise AttributeError("FR3D module has neither 'main' nor 'fr3d_search_from_query_names'")
    except SystemExit:
        pass
    except Exception:  # noqa: BLE001
        fr3d_error = traceback.format_exc()

    # --- locate the CSV FR3D wrote to OUTPUTPATH ---------------------------
    expected_csv = os.path.join(raw_dir, query_name.replace(" ", "_") + ".csv")
    csv_path = expected_csv if os.path.isfile(expected_csv) else ""
    if not csv_path and os.path.isdir(raw_dir):
        csvs = sorted(
            (os.path.join(raw_dir, f) for f in os.listdir(raw_dir) if f.lower().endswith(".csv")),
            key=lambda p: os.path.getmtime(p),
            reverse=True,
        )
        if csvs:
            csv_path = csvs[0]

    if not csv_path:
        detail = fr3d_error or "FR3D produced no CSV output"
        _fail("FR3D search produced no candidate CSV.\n%s" % detail)

    candidate_count = 0
    try:
        with open(csv_path, "r", encoding="utf-8") as handle:
            candidate_count = max(0, sum(1 for _ in handle) - 1)
    except Exception:
        pass

    result = {
        "csv_path": csv_path,
        "candidate_count": candidate_count,
        "query_name": query_name,
        "run_dir": run_dir,
        "executed_query_path": executed_path,
        "changed_fields": changed_fields,
        "entrypoint_mode": entrypoint_mode,
        "elapsed_seconds": round(time.time() - started, 3),
        "fr3d_error": fr3d_error,
        "shim_fixed_lines": shim_fixes,
        "shim_bindings": shim_bindings,
    }
    print(json.dumps({"ok": True, "result": result}))


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        _fail("Unexpected runner failure: %s\n%s" % (exc, traceback.format_exc()))
