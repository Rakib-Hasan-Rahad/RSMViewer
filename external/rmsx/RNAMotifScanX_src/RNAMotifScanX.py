#!/usr/bin/env python

import os, glob, sys

if len(sys.argv) < 4:
	print "Error in executing the script: Argument missing!"
	print "Usage: RNAMotifScanX query target_DIR" 
	sys.exit(0)

query_name = sys.argv[1]
query_name = os.path.abspath(query_name)
target_dir = sys.argv[2]
target_dir = os.path.abspath(target_dir)
output_file = sys.argv[3]
output_file = os.path.abspath(output_file)

for infile in glob.glob(os.path.join(target_dir, '*')):
	print "Aligning query against:	%s" % (infile)
	os.system("./RNAMotifAlignX -i %s -j %s -e -u 0.3 >> %s" % (query_name, infile, output_file))
