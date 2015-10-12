#!/usr/bin/env python

import os
import dendropy
from readTable import readTable
def readQuartetTable(gt):
	src_fpath = os.path.expanduser(os.path.expandvars(gt))
	frq = readTable(src_fpath)
	return frq
