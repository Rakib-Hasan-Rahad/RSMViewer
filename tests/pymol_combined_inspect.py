import os, sys, re
from pathlib import Path
sys.path.insert(0, str(Path(os.environ['RSMVIEWER_ROOT']).resolve()))
from pymol import cmd
from rsmviewer.plugin import __init_plugin__
from rsmviewer.gui import get_gui
__init_plugin__(None)
cmd.do('rmv_fetch 1S72, 1FFK')
cmd.do('rmv_db RNA3DMotifAtlas,Rfam')
gui=get_gui()
for sid, table in gui.annotation_tables.items():
    both=[]
    for row in table.rows:
        labels={source: ' | '.join(row.source_hierarchy.get(source,row.source_annotations.get(source,()))) for source in row.source_annotations}
        if len(labels)>=2 and any('SARCIN' in value.upper() for value in labels.values()):
            both.append((row.motif_id, labels))
    print('STRUCTURE', sid, 'ROWS', len(table.rows), 'SHARED_SARCIN', len(both))
    print(both[:10])
cmd.quit()
