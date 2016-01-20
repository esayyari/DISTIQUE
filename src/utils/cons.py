#!/usr/bin/env python
import dendropy
import sys
import os
from optparse import OptionParser
import numpy as np
import itertools
import random
import subprocess
WS_LOC_SHELL= os.environ['WS_HOME']+'/DISTIQUE/src/shell'
WS_LOC_FM = os.environ['WS_HOME']+'/fastme-2.1.4/src'
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-g","--gene",dest="gt",type="string",
		help="read genetrees from FILENAME")
parser.add_option("-o","--output",dest="out",type="string",
		help="the PATH to write the generated files")
parser.add_option("-t","--threshold",dest="thr",type=float,
		help="the minimum frequency that consensus will use. Default is 0.5",default=0.5)

(options,args) = parser.parse_args()
gt = options.gt
thr=options.thr
outpath=options.out
if ( not options.gt  or not options.out):
	sys.exit("Please enter genetrees file, and output folder location")

src_fpath = os.path.expanduser(os.path.expandvars(gt))

trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')

con_tree = trees.consensus(min_freq=thr)   



print "writing the resulting tree as: "+outpath+"/distance.d_cons_tree.nwk"
ftmp=tempfile.mkstemp(suffix='.nwk', prefix="distance.d_cons_tree", dir=outpath, text=None)
con_tree.write(path=ftmp[1],schema="newick",suppress_rooting=True)

                
	


