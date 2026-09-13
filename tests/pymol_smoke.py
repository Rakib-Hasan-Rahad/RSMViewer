import sys
import os
from pathlib import Path

workspace_root = Path(os.environ["RSMVIEWER_ROOT"]).resolve()
sys.path.insert(0, str(workspace_root))

from pymol import cmd

from rsmviewer.plugin import __init_plugin__


__init_plugin__(None)
registered = set(cmd.keyword)
required = {"rmv_fetch", "rmv_db", "rmv_select", "rmv_list", "rmv_view", "rmv_create_object", "rmv_super"}
removed = {"rmv_load_motif", "rmv_summary", "rmv_show", "rmv_user"}
missing = required - registered
unexpected = removed & registered
if missing:
    raise RuntimeError(f"Missing commands: {sorted(missing)}")
if unexpected:
    raise RuntimeError(f"Removed commands still registered: {sorted(unexpected)}")
print("RSMViewer PyMOL integration smoke: OK")
cmd.quit()
