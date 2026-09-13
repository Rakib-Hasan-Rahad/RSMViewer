import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["RSMVIEWER_ROOT"]).resolve()))

from pymol import cmd
from rsmviewer.plugin import __init_plugin__


__init_plugin__(None)
cmd.do("rmv_fetch 1S72")
cmd.do("rmv_fetch 1FFK")
cmd.do("rmv_db RNA3DMotifAtlas")
cmd.do("rmv_select SR, all, RNA3DMotifAtlas, as multi_SR")
cmd.do("rmv_view multi_SR")
cmd.do("rmv_create_object multi_SR")
cmd.do("rmv_super multi_SR")
cmd.do("rmv_save multi_SR cif")
print("RSMViewer multi-structure workflow smoke: OK")
cmd.quit()
