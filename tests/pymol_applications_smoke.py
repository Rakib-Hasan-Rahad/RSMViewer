import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(os.environ['RSMVIEWER_ROOT']).resolve()))
from pymol import cmd
from rsmviewer.plugin import __init_plugin__

__init_plugin__(None)
# Applications 1-2
cmd.do('rmv_fetch 1S72')
cmd.do('rmv_db RNA3DMotifAtlas')
cmd.do('rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as app1_SR')
cmd.do('rmv_list app1_SR')
if not cmd.get_names('all'):
    raise RuntimeError('Application 1 did not load a structure')
cmd.do('rmv_view app1_SR')
# Applications 3-4
cmd.do('rmv_fetch 1FFK')
cmd.do('rmv_db RNA3DMotifAtlas,Rfam')
cmd.do('rmv_select SR, 1S72 and 1FFK, RNA3DMotifAtlas and Rfam, as app3_SR')
cmd.do('rmv_select KT, all, RNA3DMotifAtlas and Rfam, as app3_KT')
cmd.do('rmv_select CL, all, RNA3DMotifAtlas and Rfam, as app3_CL')
cmd.do('rmv_select EL, all, RNA3DMotifAtlas and Rfam, as app3_EL')
cmd.do('rmv_view app3_SR,app3_KT,app3_CL,app3_EL')
cmd.do('rmv_create_object app1_SR')
cmd.do('rmv_super app1_SR')
# Application 5-6 grammar against available annotations
cmd.do('rmv_select SR, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as app_FP')
cmd.do('rmv_select SR, 1S72, RNA3DMotifAtlas and not RNAMotifScanX, as app_FN')
print('RSMViewer Applications 1-6 command smoke: OK')
cmd.quit()
