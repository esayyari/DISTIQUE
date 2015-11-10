#!/usr/bin/env python
import dendropy
import sys
import os
import numpy as np
import itertools
import subprocess
import printTools as pr
import tableManipulationTools as tbs
import anchoredTableTools as atbs
import toolsTreeTaxa as tstt
import timer as tm
from optparse import OptionParser
WS_LOC_SHELL= os.environ['WS_HOME']+'/DISTIQUE/src/shell'
WS_LOC_FM = os.environ['WS_HOME']+'/fastme-2.1.4/src'

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--file", dest="filename", type="string",
              help="read quartet table from FILENAME")
parser.add_option("-g","--gene",dest="gt",type="string",
		help="read genetrees from FILENAME")
parser.add_option("-o","--output",dest="out",type="string",
		help="the PATH to write the generated files")
parser.add_option("-t","--threshold",dest="thr",type=float,
		help="the minimum frequency that consensus will use. Default is 0.5",default=0.5)
parser.add_option("-v","--verbose",dest="verbose",
		help="Verbose",default=1)
parser.add_option("-a","--achs",dest="a",type="string",
		help="Anchors that the table will be based on")
parser.add_option("-s","--sp",dest="sp",
		help="species tree")
parser.add_option("-n","--numStep",dest="num",type="int",
		help="The number of anchoring steps, default is 2",default=2)

(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
outpath = options.out
thr = options.thr
sp = options.sp
num = options.num
if (options.a):
	ac = sorted(options.a.split(','))
	randomSample=False
else:
	randomSample=True
if options.filename:
	readFromFile = True
else:
	readFromFile = False
print readFromFile
verbose=options.verbose
if ( not options.gt  or not options.out):
	sys.exit("Please enter genetrees file, and output folder location")

src_fpath = os.path.expanduser(os.path.expandvars(gt))

tm.tic()
print "reading trees"
trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')
tm.toc()

tm.tic()
print "majority consensus"
con_tree = trees.consensus(min_freq=thr)   
tm.toc()

taxa = list()
for e in con_tree.leaf_nodes():
	taxa.append(e.taxon.label)

print taxa
n = len(con_tree.leaf_nodes())
if verbose:
	print "Number of taxa is: " + str(n)
if readFromFile:
	print "computing the distance table, reading from file"
	tm.tic()
	atbs.anchoredDistance(achs=ac,qfile=filename,outfile=outpath+'/distancet.d')
	tm.toc()
else:
	print "computing the distance table, anchoring seperately"
	tm.tic()
	atbs.anchoredDistance(achs=ac,gt=trees,wrkPath=outpath,outfile=outpath+'/distancet.d',taxa=taxa)
	tm.toc()
subprocess.call([WS_LOC_FM+"/fastme", "-i",outpath+"/distancet.d","-w","none","-o",outpath+"/distance.d_fastme_tree.nwk"])
if verbose:
	
	print "writing the resulting tree as: "+outpath+"/distance.d_fastme_tree.nwk"
res=tstt.compareAnchoredRes(outpath+'/distance.d_fastme_tree.nwk',taxa,ac,sp,outpath)
f1 = open(outpath+'/res.txt','a')
print >>f1, res
f1.close()
