#!/usr/bin/env python

import dendropy
import copy
import sys
import os
import subprocess
from readTable import readTable
from set_path import set_path
def findQuartetTable(trees,origKeys,keyType,tmpPath,verbose):
	
	WS_LOC_SH = os.environ['WS_HOME']+'/DISTIQUE/src/shell/'
	out_path = src_fpath = os.path.expanduser(os.path.expandvars(tmpPath))
	taxa = set()
	if keyType==1:
		for key in origKeys:	
			l = key.split("/")
			taxa = taxa | set(l)
	else:
		taxa = origKeys
	i = 0	
	treeList = dendropy.TreeList()
	for tree in trees:
		t = copy.deepcopy(tree)
		t.retain_taxa_with_labels(taxa, update_bipartitions=False, suppress_unifurcations=True)
		treeList.append(t)
	outfile = out_path + "/partialQuartetTable.nwk"
	treeList.write(path = outfile,schema="newick")
	command = WS_LOC_SH+"compute.quartet.freq.sh -o "+out_path+" -w 1 -f "+outfile

	os.system(command)
	frq = 	readTable(out_path+'/quartet_tmp.q')
	if verbose:
		print "The quartet table has been computed"
	return frq
