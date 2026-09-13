import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["RSMVIEWER_ROOT"]).resolve()))
from pymol import cmd
from rsmviewer.plugin import __init_plugin__

__init_plugin__(None)
cmd.do("rmv_fetch 1S72, 1FFK")
cmd.do("rmv_db RNA3DMotifAtlas,Rfam")
cmd.do("rmv_select SARCIN-RICIN, 1S72 and 1FFK, RNA3DMotifAtlas and Rfam, as group_SR")
if not cmd.get_names("all"):
    raise RuntimeError("No PyMOL objects after combined selection")
print("RSMViewer combined selection smoke: OK")
cmd.quit()
