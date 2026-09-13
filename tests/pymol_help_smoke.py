import os
import sys
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["RSMVIEWER_ROOT"]).resolve()))
from pymol import cmd
from rsmviewer.plugin import __init_plugin__

__init_plugin__(None)
output = StringIO()
with redirect_stdout(output):
    cmd.do("rmv_help")
text = output.getvalue()
for required in ("RSMViewer v2.0.0", "Updated: 10 September 2026", "rmv_select", "rmv_list", "rmv_create_object", "rmv_super"):
    if required not in text:
        raise RuntimeError(f"Missing help text: {required}")
for removed in ("rmv_load_motif", "rmv_summary", "rmv_show", "rmv_user", "rmv_db 3"):
    if removed in text:
        raise RuntimeError(f"Legacy help text found: {removed}")
print("RSMViewer help smoke: OK")
cmd.quit()
