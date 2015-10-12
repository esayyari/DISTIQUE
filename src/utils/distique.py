#!/usr/bin/env python
from averageQuartetTables import averageQuartetTables
import dendropy
import sys
import os
from optparse import OptionParser
from  prodDistance import prodDistance
from printDistanceTable import printDistanceTable
from minDistance import minDistance
import numpy as np
import itertools
import random
from getQuartetKeys import generateKey
from findPartialQuartetFreq import partialQuartetTable
from addQuartetTables import addQuartetTables
from resolvePolytomies import resolvePolytomy
import subprocess
from getDistanceTable import distanceTable
from  generateTaxaList import getTaxaList
from labelNodes import labelNodes
from findPolytomies import findPolytomies
from readTable import readTable
WS_LOC_SHELL= os.environ['WS_HOME']+'/DISTIQUE/src/shell'
WS_LOC_FM = os.environ['WS_HOME']+'/fastme-2.1.4/src'
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--file", dest="filename", type="string",
              help="read data from FILENAME")
parser.add_option("-g","--gene",dest="gt",type="string",
		help="read genetrees from FILENAME")
parser.add_option("-o","--output",dest="out",type="string",
		help="the PATH to write the generated files")
parser.add_option("-t","--threshold",dest="thr",
		help="the minimum frequency that consensus will use. Default is 0.5",default=0.5)
parser.add_option("-n","--numToStop",dest="numToStop",
		help="The number of steps that convergence criteria should be meet until the successful convergence. Default is 10",default=10)
parser.add_option("-m","--method",dest="method",
		help="The method to compute the distance of taxa. The default is prod.",default="prod")
parser.add_option("-x","--numMax",dest="numMax",
		help="The maximum number of steps for computing the average quartet table. The default is 100", default=100)
parser.add_option("-e","--epsilon",dest="epsilon",
		help="The threshold of convergence for computing the average quartet table using KL divergence. The default is 0.01 ",default=0.01)
parser.add_option("-v","--verbose",dest="verbose",
		help="Verbose",default=1)
parser.add_option("-r","--readFromFile",dest="readFromFile",
		help="Set it to 1 if you already computed the quartet table. If you want the code to compute the quartet tables set it to 0. Default is 0",default=0)
parser.add_option("-p","--maxPoly",dest="maxPoly",
		help="This indicates the maximum order of polytomy in the tree that you want the code to compute it using partial quartet tables. If you want the code to compute the whole quartet table your could set it to 0. Default is 20",default=20)

(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
outpath = options.out
thr = options.thr
thr=options.thr
numToStop = options.numToStop
numMax = options.numMax
eps = options.epsilon
verbose=options.verbose
if options.readFromFile == 1:
	readFromFile = True
else:
	readFromFile = False
maxPossiblePolyOrder = options.maxPoly
method = options.method

if ( not options.gt or not options.filename or not options.out):
	sys.exit("Please enter genetrees file and quartetTable file location")

frq = readTable(filename)

src_fpath = os.path.expanduser(os.path.expandvars(gt))

trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')

con_tree = trees.consensus(min_freq=thr)   

labelNodes(con_tree)

con_tree.write(path="consensusTree.nwk",schema="newick") 

(to_resolve,maxPolyOrder) = findPolytomies(con_tree)

for e in con_tree.postorder_node_iter():
	if e in to_resolve:
		val = to_resolve[e]
		(taxa_list,taxa_inv) =  getTaxaList(to_resolve[e])
		if maxPolyOrder > maxPossiblePolyOrder and readFromFile:
			quartTable = averageQuartetTables(limit=eps,NumToStop = numToStop, NumMax = numMax,ListTaxa=taxa_list,QTable=frq,QtablePath=filename,QtableReady=True,Inv=taxa_inv,V=verbose)
		else:
			quartTable = averageQuartetTables(limit=eps,NumToStop = numToStop, NumMax = numMax,ListTaxa=taxa_list,QTable=frq,QtablePath=filename,QtableReady=False,workingPath = outpath,Inv=taxa_inv,V=verbose,treeList=trees)
		distanceTable(quartTable,method,outpath+"/distancet.d")
		subprocess.call([WS_LOC_FM+"/fastme", "-i",outpath+"/distancet.d","-w","none","-o",outpath+"/distancet.d_fastme_tree.nwk"])	
		res= resolvePolytomy(outpath+"/distancet.d_fastme_tree.nwk",e,con_tree)	
		if verbose:
			print res

con_tree.write(path=outpath+"distance.d_fastme_tree.nwk",schema="newick")
tns = dendropy.TaxonNamespace()
tree1 = dendropy.Tree.get_from_path("/Users/Erfan/Documents/Research/data/mammalian/mammalian-model-species.tre","newick",taxon_namespace=tns,rooting="force-unrooted")
tree2 = dendropy.Tree.get_from_path("trees1.nwk","newick",taxon_namespace=tns,rooting="force-unrooted")
res1 = dendropy.calculate.treecompare.false_positives_and_negatives(tree1,tree2)
print res1

