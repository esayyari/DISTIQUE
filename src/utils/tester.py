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
from greedy_cons import greedy_cons
import subprocess
from getDistanceTable import distanceTable
from  generateTaxaList import getTaxaList
from labelNodes import labelNodes
from findPolytomies import findPolytomies
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--file", dest="filename", type="string",
              help="read data from FILENAME")
parser.add_option("-g","--gene",dest="gt",type="string",
		help="read genetrees from FILENAME")
parser.add_option("-o","--output",dest="out",type="string",
		help="the PATH to write the generated files")
(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
outpath = options.out


thr=0.9
if ( not options.gt or not options.filename or not options.out):
	sys.exit("Please enter genetrees file and quartetTable file location")
	
f = open(filename, 'r')
frq = dict()
for line in f:
    k=line.split()
    v = dict()
    d = k[0].split('/')
    v[d[1]] = float(k[1])
    v[d[2]] = float(k[2])
    v[d[3]] = float(k[3])
    frq[k[0]] = v
keyDict = sorted(np.unique(("/".join(frq.keys())).split("/")));
src_fpath = os.path.expanduser(os.path.expandvars(gt))
trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')
con_tree = trees.consensus(min_freq=thr)   
labelNodes(con_tree)
con_tree.write(path="consensusTree.nwk",schema="newick") 
#con_tree.print_plot() 
numToStop = 10
numMax = 100
eps = 0.01
verbose=1
readFromFile = True
maxPossiblePolyOrder = 20
(to_resolve,maxPolyOrder) = findPolytomies(con_tree)
#if len(to_resolve)!= 1:
for e in con_tree.postorder_node_iter():
	if e in to_resolve:
		val = to_resolve[e]
		(taxa_list,taxa_inv) =  getTaxaList(to_resolve[e])
		if maxPolyOrder > maxPossiblePolyOrder and readFromFile:
			quartTable = averageQuartetTables(limit=eps,NumToStop = numToStop, NumMax = numMax,ListTaxa=taxa_list,QTable=frq,QtablePath=filename,QtableReady=True,Inv=taxa_inv,V=verbose)
		else:
			quartTable = averageQuartetTables(limit=eps,NumToStop = numToStop, NumMax = numMax,ListTaxa=taxa_list,QTable=frq,QtablePath=filename,QtableReady=False,workingPath = outpath,Inv=taxa_inv,V=verbose,treeList=trees)
		distanceTable(quartTable,"prod",outpath+"/distancet.d")
		subprocess.call(["/Users/Erfan/Documents/Research/fastme-2.1.4/src/fastme", "-i",outpath+"/distancet.d","-w","none","-o",outpath+"/distancet.d_fastme_tree.nwk"])	
		res= resolvePolytomy(outpath+"/distancet.d_fastme_tree.nwk",e,con_tree)	
		print res
#else:
#	e = to_resolve.keys()
#	e = e[0]
#	distanceTable(frq,"prod",outpath+"/distancet.d")
#	subprocess.call(["/Users/Erfan/Documents/Research/fastme-2.1.4/src/fastme", "-i",outpath+"/distancet.d","-w","none","-o",outpath+"/distancet.d_fastme_tree.nwk"])
#	res= resolvePolytomy(outpath+"/distancet.d_fastme_tree.nwk",e,con_tree)
#	print res

con_tree.write(path="trees1.nwk",schema="newick")
tns = dendropy.TaxonNamespace()
tree1 = dendropy.Tree.get_from_path("/Users/Erfan/Documents/Research/data/mammalian/mammalian-model-species.tre","newick",taxon_namespace=tns,rooting="force-unrooted")
tree2 = dendropy.Tree.get_from_path("trees1.nwk","newick",taxon_namespace=tns,rooting="force-unrooted")
res1 = dendropy.calculate.treecompare.false_positives_and_negatives(tree1,tree2)
print res1

