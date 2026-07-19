#!/usr/bin/python2.6
#!python

# The script is used to generate input file for RNAMotifScanX
# The script accepts a PDB file and use MC-Annotate and RNAVIEW to parse
#	the file into base-interaction patterns 
# by Cuncong Zhong (send comments or bug reports to czhong@jvci.org)

import MapPDBSequence
import LoadSequence
import ParseStructureAnnotation
import FormattedOutput

import argparse
import os
import sys
import re


# Parse MC-Annotate and RNAVIEW files
#mca_interactions = ParseStructureAnnotation.GetMCAnnotateInteractions(result_file_mca)
rvw_interactions = ParseStructureAnnotation.GetRNAVIEWInteractions("/home/cczhong/Codes/RNAMotifScanX/PDB/1S72.pdb.out")
#merged_interactions = ParseStructureAnnotation.MergeInteractions(mca_interactions, rvw_interactions)


for single_interaction in rvw_interactions:
  print "interactions %s  %s  %s\n" %(single_interaction[0], single_interaction[1], single_interaction[2])
