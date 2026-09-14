#!/usr/bin/perl -w
use strict;

my $archive_dir = shift;
foreach(<$archive_dir/*.pdb>)	{
	print "################### Parsing $_ #########################\n";
	system "./PrepareInput.py $_ pdb_seqres.na.txt";
}
