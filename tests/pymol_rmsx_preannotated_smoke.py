import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["RSMVIEWER_ROOT"]).resolve()))
from pymol import cmd
from rsmviewer.plugin import __init_plugin__

__init_plugin__(None)
cmd.do("rmv_fetch 1a51")
cmd.do("rmv_db RNAMotifScanX")
cmd.do("rmv_list")
print("RSMViewer RMSX preannotated workflow smoke: OK")
cmd.quit()
