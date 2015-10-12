#!/usr/bin/env python

import dendropy
import copy
import sys
import os
import subprocess
from readTable import readTable
def findQuartetTable(trees,origKeys,tmpPath):
	out_path = src_fpath = os.path.expanduser(os.path.expandvars(tmpPath))
	taxa = set()
	for key in origKeys:	
		l = key.split("/")
		taxa = taxa | set(l)
	i = 0	
	treeList = dendropy.TreeList()
	for tree in trees:
		t = copy.deepcopy(tree)
		t.retain_taxa_with_labels(taxa, update_bipartitions=False, suppress_unifurcations=True)
		treeList.append(t)
	outfile = out_path + "/partialQuartetTable.nwk"
	treeList.write(path = outfile,schema="newick")
	command = "/Users/Erfan/Documents/Research/DISTIQUE/src/shell/compute.quartet.freq.sh -o "+out_path+" -w 1 -f "+outfile

	os.system(command)
	frq = 	readTable(out_path+'/quartet_tmp.q')
	return frq
