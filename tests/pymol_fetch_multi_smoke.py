import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["RSMVIEWER_ROOT"]).resolve()))
from pymol import cmd
from rsmviewer.plugin import __init_plugin__

__init_plugin__(None)
cmd.do("rmv_fetch 1S72, 1FFK")
objects = set(cmd.get_object_list("all"))
if not {"1s72", "1ffk"}.issubset(objects):
    raise RuntimeError(f"Multiple fetch failed; objects={sorted(objects)}")
print("RSMViewer multiple-fetch smoke: OK")
cmd.quit()
