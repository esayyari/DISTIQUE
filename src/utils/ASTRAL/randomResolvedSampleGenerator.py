#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import numpy as np
import itertools
import random
import subprocess
import dendropy

if __name__ == '__main__':
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option("-g","--gene",dest="gt",type="string",
			help="read genetrees from FILENAME")
	parser.add_option("-o","--output",dest="out",type="string",
			help="the PATH to write the generated files")
	parser.add_option("-t","--threshold",dest="thr",type="string",
			help="the minimum frequency that consensus will use. Default is 0.2",default=0.2)
	parser.add_option("-v","--verbose",dest="verbose",
			help="Verbose",default=1)
	parser.add_option("-n","--replicates",dest="num",type=int,
			help="The number of sampling, and resolving polytomies",default=10)
	(options,args) = parser.parse_args()
	gt = options.gt
	outpath = options.out
	thrt = options.thr
	thr = thrt.split(",")
	print thr
	verbose=options.verbose
	num=options.num
	if ( not options.gt  or not options.out):
		sys.exit("Please enter genetrees file, and output folder location")

	src_fpath = os.path.expanduser(os.path.expandvars(gt))

	trees = dendropy.TreeList.get_from_path(src_fpath, 'newick',rooting='force-unrooted')
	
	treelist = dendropy.TreeList()
	for th in thr:
		print th
		con_tree = trees.consensus(min_freq=float(th))
		for i in range(0,num):
			tree_tmp = con_tree.clone(2)
			tree_tmp.resolve_polytomies(limit=2, update_bipartitions=False, rng=random.seed())
			treelist.append(tree_tmp)


	treelist.write(path=outpath, schema="newick")
