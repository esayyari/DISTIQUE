#!/usr/bin/env python
import dendropy
import itertools
import sys
import os
from readQuartetTable import readQuartetTable
from printTools import printDistanceTableToFile
import printTools as pr
from labelNodes import labelNodes

from findAnchoredQuartetTable import findAnchoredQuartetTable
import numpy as np
def anchoredDistance(**kwargs):
	readFromTable=False
	for k,v in kwargs.iteritems():
		if k == "qfile":
			qfile=v
			readFromTable = True
		elif k == "achs":
			achs = sorted(v)
		elif k == "gt":
			gt = v
		elif k == "outfile":
			outfile = v
		elif k == "wkrPath":
			out = v
	if(readFromTable):
		frq=readQuartetTable(qfile)	
	else:
		frq=findAnchoredQuartetTable(achs,gt)	
	D = anchoredDistanceFromFrq(frq,achs)
	keyDict = sorted(np.unique((" ".join(D.keys())).split(" ")));              
	print "print distance table to file"	
	printDistanceTableToFile(D,keyDict,outfile)
	return 
def anchoredDistanceFromFrq(frq,achs):
	D = dict()
	for k,v in frq.iteritems():
		kt = sorted(k.split("/"))
		if ((achs[0] in kt ) and( achs[1] in kt)):
			s = sorted(list(set(kt)-{achs[0],achs[1]}))
			key1 = s[0]+" "+s[1]
			key2 = s[1]+" "+s[0]
			if s[0]<achs[0]:
				D[key1]=-np.log(frq[k][s[1]])
				D[key2]=D[key1]	
			elif s[0]>achs[0]:
				D[key1]=-np.log(frq[k][achs[1]])	
				D[key2]=D[key1]
	return D		 
def findAnchoredQuartets(anch,gt,out):
	src_fpath = os.path.expanduser(os.path.expandvars(gt))
	trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')
	Q = list()
	anch = sorted(anch)
	print anch
	for tree in trees:
		labelNodes(tree)
		filter = lambda taxon: True if taxon.label==anch[0] else False
		root = tree.find_node_with_taxon(filter)

		tree.reroot_at_node(root,update_bipartitions=True, suppress_unifurcations=False)
		filter = lambda taxon: True if taxon.label==anch[1] else False
		node = tree.find_node_with_taxon(filter)
		print node.label
		while node.parent_node != root:
			nodePre = node
			node = node.parent_node

			ch = node.child_nodes()
			
			if len(ch)>2:
				for ct in ch:
					if ct == nodePre:
						continue
					leafset = ct.leaf_nodes()
					if len(leafset)<2: 
						continue
					else:
						Q += listQ(leafset,anch)
						
				
			else:
				for ct in ch:
					if ct == nodePre:
						continue
					leafset = ct.leaf_nodes()
					if len(leafset)<2:
						continue
					else:
						Q += listQ(leafset,anch)
	pr.printQuartetsToFile(Q,out+'/tmp.q')
	return 
def findAnchoredQuartetTables(outfile,anch,
def listQ(leafset,anch):
	Q = list()
	for l in list(itertools.combinations(leafset,2)):
		lt = sorted([l[0].label,l[1].label])
		Q.append(anch[0]+' '+anch[1]+'|'+lt[0]+' '+lt[1])
	return Q	
				
