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
from findQuartetTable import findQuartetTable
from findTrueAverageTable import findTrueAverageTable
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
parser.add_option("-n","--numToStop",dest="numToStop",type=int,
		help="The number of steps that convergence criteria should be meet until the successful convergence. Default is 10",default=10)
parser.add_option("-m","--method",dest="method",type=str,
		help="The method to compute the distance of taxa. The default is prod.",default="prod")
parser.add_option("-x","--numMax",dest="numMax",type=int,
		help="The maximum number of steps for computing the average quartet table. The default is 100", default=100)
parser.add_option("-e","--epsilon",dest="epsilon",type=float,
		help="The threshold of convergence for computing the average quartet table using KL divergence. The default is 0.01 ",default=0.01)
parser.add_option("-v","--verbose",dest="verbose",
		help="Verbose",default=1)
parser.add_option("-r","--readFromFile",dest="readFromFile",
		help="Set it to 1 if you already computed the quartet table. If you want the code to compute the quartet tables set it to 0. Default is 0",default=0)

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
elif options.filename:
	readFromFile = True
else:
	readFromFile = False
method = options.method

if ( not options.gt  or not options.out):
	sys.exit("Please enter genetrees file, and output folder location")
if readFromFile:
	frq = readTable(filename)

src_fpath = os.path.expanduser(os.path.expandvars(gt))

trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')

con_tree = trees.consensus(min_freq=thr)   

labelNodes(con_tree)

con_tree.write(path="consensusTree.nwk",schema="newick") 

(to_resolve,maxPolyOrder) = findPolytomies(con_tree)
taxa = list()
for e in con_tree.leaf_nodes():
	taxa.append(e.taxon.label)
maxPossiblePoly = pow(len(to_resolve),0.25)*2.5*maxPolyOrder
n = len(con_tree.leaf_nodes())
if verbose:
	print "Number of taxa is: " + str(n)
	print "the number of polytomies is: "+str(len(to_resolve))
	print "the maximum order of polytomies is: "+str(maxPolyOrder)
	print "the maxPossiblePoly is: "+str(maxPossiblePoly)
if n<maxPossiblePoly:
	if verbose:
		print "computing the total quartet table"
<<<<<<< HEAD
	if readFromFile:
		frq = readTable(filename)
	else:
		frq = findQuartetTable(trees,taxa,0,outpath,verbose)
if maxPolyOrder< n-4:
=======
	frq = findQuartetTable(trees,taxa,0,outpath,verbose)
if len(to_resolve)<n:
>>>>>>> fc8361a5592fdbd6b0903caeabfc15dec790fd07
	for e in con_tree.postorder_node_iter():
		if e in to_resolve:
			val = to_resolve[e]
			(taxa_list,taxa_inv) =  getTaxaList(to_resolve[e])
			if readFromFile:
				if verbose:
					print "reading the quartet table from file: "+filename
				quartTable = averageQuartetTables(limit=eps,NumToStop = numToStop, NumMax = numMax,ListTaxa=taxa_list,QTable=frq,QtablePath=filename,QtableReady=True,Inv=taxa_inv,V=verbose,KeyType = 1)
			elif n<maxPossiblePoly:
				if verbose:
					print "using precomputed quartet table"
				quartTable = averageQuartetTables(limit=eps,NumToStop = numToStop, NumMax = numMax,ListTaxa=taxa_list,QTable=frq,QtableReady=False,workingPath = outpath,Inv=taxa_inv,V=verbose,treeList=trees,KeyType = 1)
			else:
				if verbose:
					print "computing the partial quartet table in each step"
				quartTable = averageQuartetTables(limit=eps,NumToStop = numToStop, NumMax = numMax,ListTaxa=taxa_list,QtableReady=False,workingPath = outpath,Inv=taxa_inv,V=verbose,treeList=trees,KeyType=1)
			
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
else:
	if verbose:
		print "computing distance matrix using the method: "+method
	distanceTable(frq,method,outpath+"/distancet.d")
	if verbose:
                print "writing the resulting tree as: "+outpath+"/distance.d_distique_tree.nwk"
	subprocess.call([WS_LOC_FM+"/fastme", "-i",outpath+"/distancet.d","-w","none","-o",outpath+"/distance.d_distique_tree.nwk"])	
		


