#!/usr/bin/env python

import dendropy
import sys
import os
import itertools
from optparse import OptionParser


def expandname(filename):
	src_fpath = os.path.expanduser(os.path.expandvars(filename))
	if not os.path.exists(src_fpath):
		sys.stderr.write('Not found: "%s"' % src_fpath)
	return src_fpath


def findInternalBranches(Tree):
	taxon_namespace = Tree.taxon_namespace
	for edge in Tree.postorder_internal_edge_iter(exclude_seed_edge=True):
		to_print = edge.bipartition.leafset_as_newick_string(taxon_namespace)
		to_print = to_print.replace(";", "")
		to_print = to_print.replace(" ", "")
		to_print = to_print.replace("),(", "|")
		to_print = to_print.replace("(", "")
		to_print = to_print.replace(")", "")
		g = to_print.split("|")
		g[0] = ",".join(sorted(g[0].split(",")))
		g[1] = ",".join(sorted(g[1].split(",")))
		to_print = "|".join(sorted([g[1], g[0]]))
		tmpList.append(to_print)
	return tmpList


if __name__ == "__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option("-s", "--species", dest="sp", type="string",
					  help="read species tree from FILENAME")
	parser.add_option("-o", "--output", dest="out", type="string",
					  help="the PATH to write the generated files")
	(options, args) = parser.parse_args()
	sp = options.sp
	out = options.out

	sp = expandname(sp)

	tree = dendropy.Tree.get(path=sp, schema="newick", store_tree_weights=True, rooting='force-unrooted')
	tree.encode_bipartitions()
	tmpList = list()
	f = open(out, 'w')
	tmpList = findInternalBranches(tree)
	for c in tmpList:
		print >> f,c
	f.close()
