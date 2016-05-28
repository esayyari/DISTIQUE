#!/usr/bin/env python

import dendropy
import sys
import os
from optparse import OptionParser
import numpy as np
import itertools
import random
import subprocess
import tempfile
def expandName(filename):
        src_fpath = os.path.expanduser(os.path.expandvars(filename))
        if not os.path.exists(src_fpath):
                sys.stderr.write('Not found: "%s"' % src_fpath)
        return src_fpath

if __name__ == "__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option("-s","--species",dest="sp",type="string",
			help="read species tree from FILENAME")
	parser.add_option("-o","--output",dest="out",type="string",
			help="the PATH to write the generated files")
	(options,args) = parser.parse_args()
	sp = options.sp
	out = options.out

	sp = expandName(sp)

	tree = dendropy.Tree.get(path=sp,schema="newick",store_tree_weights=True)
	tree.encode_bipartitions()
	taxon_namespace = tree.taxon_namespace
	f = open(out,'w')
	for edge in tree.postorder_internal_edge_iter():
		if edge.length is not None:
			to_print = edge.bipartition.leafset_as_newick_string(taxon_namespace)
			to_print = to_print.replace(";",":")
			to_print = to_print.replace(" ","")
			to_print = to_print.replace("),(","}|{")
			to_print = to_print.replace("(","{")
			to_print = to_print.replace(")","}")
			to_print = to_print.replace("{{","{")
			to_print = to_print.replace("}}","}")
			print >>f, to_print, edge.length
	
