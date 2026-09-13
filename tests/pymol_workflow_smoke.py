import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["RSMVIEWER_ROOT"]).resolve()))

from pymol import cmd
from rsmviewer.plugin import __init_plugin__


__init_plugin__(None)
cmd.do("rmv_fetch 1S72")
cmd.do("rmv_db RNA3DMotifAtlas")
cmd.do("rmv_list")
cmd.do("rmv_select SR, 1S72, RNA3DMotifAtlas, as live_SR")
cmd.do("rmv_list live_SR")
print("RSMViewer named workflow smoke: OK")
cmd.quit()
