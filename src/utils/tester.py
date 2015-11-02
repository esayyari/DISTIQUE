#!/usr/bin/env python
import dendropy
import sys
import os
import numpy as np
import itertools
from readTable import readTable
from anchoredDistance import anchoredDistance
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
parser.add_option("-a","--anchors",dest=a,
		help="Anchors that the table will be based on")
(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
outpath = options.out
thr = options.thr
achs = sorted("".split(options.a))
if options.filename:
	readFromFile = True
else:
	readFromFile = False
verbose=options.verbose
if ( not options.gt  or not options.out):
	sys.exit("Please enter genetrees file, and output folder location")

src_fpath = os.path.expanduser(os.path.expandvars(gt))

trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')

con_tree = trees.consensus(min_freq=thr)   

(to_resolve,maxPolyOrder) = findPolytomies(con_tree)
taxa = list()
for e in con_tree.leaf_nodes():
	taxa.append(e.taxon.label)
n = len(con_tree.leaf_nodes())
if verbose:
	print "Number of taxa is: " + str(n)
	print "the number of polytomies is: "+str(len(to_resolve))
	print "the maximum order of polytomies is: "+str(maxPolyOrder)
	print "the maxPossiblePoly is: "+str(maxPossiblePoly)
if verbose:
	print "computing the total quartet table"
if readFromFile:
	frq = readTable(filename)
else:
	frq = findQuartetTable(trees,taxa,0,outpath,verbose)
for e in con_tree.postorder_node_iter():
	if e in to_resolve:
		val = to_resolve[e]
		(taxa_list,taxa_inv) =  getTaxaList(to_resolve[e])
		if verbose:
			print "computing the partial quartet table"
		
		quartTable = findTrueAverageTable(frq,taxa_list,av,po)
		if verbose:
			print "computing distance table using the method: "+str(method)
		distanceTable(quartTable,method,outpath+"/distancet.d")
		subprocess.call([WS_LOC_FM+"/fastme", "-i",outpath+"/distancet.d","-w","none","-o",outpath+"/distancet.d_fastme_tree.nwk"])
		if verbose:
			print "starting to resolve polytomy"	
		res= resolvePolytomy(outpath+"/distancet.d_fastme_tree.nwk",e,con_tree,verbose)	
		if verbose:
			print res
if verbose:
	print "writing the resulting tree as: "+outpath+"/distance.d_distique_tree.nwk"
con_tree.write(path=outpath+"/distance.d_distique_tree.nwk",schema="newick")
