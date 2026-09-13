import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(os.environ['RSMVIEWER_ROOT']).resolve()))
from pymol import cmd
from rsmviewer.plugin import __init_plugin__
from rsmviewer.gui import get_gui
__init_plugin__(None)
cmd.do('rmv_fetch 1S72')
cmd.do('rmv_db RNA3DMotifAtlas,Rfam')
gui=get_gui()
for sid, table in gui.annotation_tables.items():
    rows=[r for r in table.rows if any('SARCIN' in l.upper() for vals in r.source_hierarchy.values() for l in vals)]
    print(sid, len(table.rows), len(rows), sorted({s for r in table.rows for s in r.source_annotations}))
cmd.quit()
