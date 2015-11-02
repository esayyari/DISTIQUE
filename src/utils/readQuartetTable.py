#!/usr/bin/env python

import os
import dendropy
from readTable import readTable
def readQuartetTable(gt):
	src_fpath = os.path.expanduser(os.path.expandvars(gt))
	frq = readTable(src_fpath)
	for k in frq:
		sz = sum(frq[k].values())
		for k2 in frq[k]:
			frq[k][k2] /= sz
	return frq
