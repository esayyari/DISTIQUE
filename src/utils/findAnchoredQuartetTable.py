#!/usr/bin/env python
import dendropy
import os
def findAnchoredQuartetTable(anch,gt):
	src_fpath = os.path.expanduser(os.path.expandvars(gt))
	trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')


