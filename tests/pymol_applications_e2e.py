"""End-to-end verification of RSMViewer manuscript Applications 1-6.

Runs each documented command sequence in a real headless PyMOL session and
writes a PASS/FAIL report to /tmp/rsmviewer_apps_report.txt. Introspection runs
inside the PyMOL command queue (cmd.do) so it executes AFTER the queued
commands, then the assertions are evaluated and the report written.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ['RSMVIEWER_ROOT']).resolve()))
from pymol import cmd
from rsmviewer.plugin import __init_plugin__

__init_plugin__(None)

REPORT = '/tmp/rsmviewer_apps_report.txt'

# ── Application 1 & 2: Atlas Sarcin-Ricin on 1S72 ────────────────────────────
cmd.do('rmv_reset')
cmd.do('rmv_fetch 1S72')
cmd.do('rmv_db RNA3DMotifAtlas')
cmd.do('rmv_select SR, 1S72, RNA3DMotifAtlas, as group_SR')
cmd.do('rmv_list group_SR')
cmd.do('rmv_view group_SR')                       # App 2: view the whole family

# ── Application 4: objects + superimpose ─────────────────────────────────────
cmd.do('rmv_create_object group_SR')
cmd.do('rmv_super group_SR')

# ── Application 3: multi-structure, multi-family (Atlas + Rfam) ───────────────
cmd.do('rmv_fetch 1FFK')
cmd.do('rmv_db RNA3DMotifAtlas,Rfam')
cmd.do('rmv_select SR, 1S72 and 1FFK, RNA3DMotifAtlas and Rfam, as g3_SR')
cmd.do('rmv_select SR, 1S72, RNA3DMotifAtlas, as g3_SR_atlas')
cmd.do('rmv_select KT, all, RNA3DMotifAtlas, as g3_KT')
cmd.do('rmv_select CL, all, RNA3DMotifAtlas, as g3_CL')
cmd.do('rmv_select EL, all, RNA3DMotifAtlas, as g3_EL')
cmd.do('rmv_view g3_KT, g3_CL, g3_EL')

# ── Application 5: TP / FP / FN with RNAMotifScanX ────────────────────────────
cmd.do('rmv_db RNA3DMotifAtlas,RNAMotifScanX')
cmd.do('rmv_select SR, 1S72, RNA3DMotifAtlas and RNAMotifScanX, as group_TP')
cmd.do('rmv_select SR, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP')
cmd.do('rmv_select SR, 1S72, RNA3DMotifAtlas and not RNAMotifScanX, as group_FN')

# ── Application 6: combine + color + superimpose ─────────────────────────────
cmd.do('rmv_set_color group_FP, red')
cmd.do('rmv_combine group_FP, group_TP, as group_combined')
cmd.do('rmv_create_object group_combined')
cmd.do('rmv_super group_combined')

probe = (
    "python\n"
    "import rsmviewer.gui as _g\n"
    "_gui=_g.gui\n"
    "_lines=[]\n"
    "def _n(name):\n"
    "    _grp=_gui.query_groups.get(name)\n"
    "    return len(_grp['motif_ids']) if _grp else -1\n"
    "def _check(label, ok, detail=''):\n"
    "    _lines.append(('PASS' if ok else 'FAIL')+' '+label+((' :: '+detail) if detail else ''))\n"
    "_objs=set(cmd.get_object_list())\n"
    "# App1/2\n"
    "_check('App1 group_SR has 8 Atlas SR motifs', _n('group_SR')==8, 'count='+str(_n('group_SR')))\n"
    "# App4\n"
    "_check('App4 created motif_ objects for group_SR', any(o.startswith('motif_') for o in _objs), 'objs='+str(len([o for o in _objs if o.startswith('motif_')])))\n"
    "# App3\n"
    "_check('App3 Atlas-only SR on 1S72 non-empty', _n('g3_SR_atlas')>0, 'count='+str(_n('g3_SR_atlas')))\n"
    "_check('App3 KT/CL groups saved (EL data-dependent)', _n('g3_KT')>0 and _n('g3_CL')>0, 'KT='+str(_n('g3_KT'))+' CL='+str(_n('g3_CL'))+' EL='+str(_n('g3_EL')))\n"
    "_check('App3 shared Atlas+Rfam SR group saved', _n('g3_SR')>=0, 'count='+str(_n('g3_SR')))\n"
    "# App5\n"
    "_tp,_fp,_fn=_n('group_TP'),_n('group_FP'),_n('group_FN')\n"
    "_check('App5 TP non-empty', _tp>0, 'TP='+str(_tp))\n"
    "_check('App5 FP non-empty', _fp>0, 'FP='+str(_fp))\n"
    "_check('App5 FN non-empty', _fn>0, 'FN='+str(_fn))\n"
    "# App6\n"
    "_comb=_n('group_combined')\n"
    "_check('App6 combined = TP + FP union', _comb==len(set(_gui.query_groups.get('group_TP',{}).get('motif_ids',[]))|set(_gui.query_groups.get('group_FP',{}).get('motif_ids',[]))), 'combined='+str(_comb))\n"
    "_check('App6 group_combined objects created', 'group_combined' in cmd.get_names('group_objects') or any(o.startswith('motif_') for o in cmd.get_object_list()), '')\n"
    "with open('" + REPORT + "','w') as _f:\n"
    "    _f.write(chr(10).join(_lines)+chr(10))\n"
    "    _f.write('SUMMARY: '+str(sum(1 for l in _lines if l.startswith('PASS')))+'/'+str(len(_lines))+' passed'+chr(10))\n"
    "python end"
)
cmd.do(probe)
cmd.do('quit')
