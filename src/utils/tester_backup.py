#!/usr/bin/env python

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
	
to_resolve=dict()
tmp = list()
taxa = set(con_tree.leaf_nodes())
p = 0
j = 0
samplingNum = 100
for e in con_tree.postorder_node_iter():
	n = len(e.child_nodes())
	tmp_set = set()
	if n>3 or ( n==3 and e.parent_node>0 ):
		e.label = "poly"+str(j)
		j += 1
		v = dict()
		i = 0
		for c in e.child_node_iter():
			c.label = "child"+str(i)
			if len(tmp)>0:
				t = tmp.pop()
				v["child"+str(i)] = t
				tmp_set = tmp_set | t
				i += 1
		tmp.append(tmp_set);
		if len(taxa-tmp_set)>0:
			e.parent_node.label = str("child"+str(n))
			t = taxa-tmp_set
			v["child"+str(n)] = t
		to_resolve[e] = v
	else:
		for i in range(0,n):
			if len(tmp)>0:
				tmp_set = tmp_set | tmp.pop()
		if e.is_leaf():
			tmp_set.add(e)
		tmp.append(tmp_set)

for e, val in to_resolve.iteritems():
        (taxa_list,taxa_inv) =  getTaxaList(val)
	chn = e.child_nodes()
	for ch in chn:
		print ch.label
	if e.parent_node>0:
		print e.parent_node.label
	for a in range(0,samplingNum):
		origKeys = generateKey(val,taxa_list)
		partialTable1	     = partialQuartetTable(frq,origKeys,taxa_inv)
		if a>0:
			quartTable=addQuartetTables(partialTable1,quartTable)
		else:
			quartTable = partialTable1
	distanceTable(partialTable1,"prod",outpath+"/distancet.d")
	subprocess.call(["/Users/Erfan/Documents/Research/fastme-2.1.4/src/fastme", "-i",outpath+"/distancet.d","-w","none","-o",outpath+"/distancet.d_fastme_tree.nwk"])	
	resolvePolytomy(outpath+"/distancet.d_fastme_tree.nwk",e)
