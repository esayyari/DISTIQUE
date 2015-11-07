#!/usr/bin/env python
import dendropy
import itertools
import sys
import os
from readQuartetTable import readQuartetTable
import printTools as pr
from labelNodes import labelNodes

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
		elif k == "wrkPath":
			out = v
		elif k == "taxa":
			taxa = v
	if(readFromTable):
		frq=readQuartetTable(qfile)	
	else:
		frq=findAnchoredDistanceTable(achs,gt,taxa,out)	
	D = anchoredDistanceFromFrq(frq,achs)
	keyDict = sorted(np.unique((" ".join(D.keys())).split(" ")));              
	print "print distance table to file"	
	pr.printDistanceTableToFile(D,keyDict,outfile)
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
def findAnchoredQuartets(anch,gt,taxa,out):
	src_fpath = os.path.expanduser(os.path.expandvars(gt))
	trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')
	Q = list()
	anch = sorted(anch)
	n = len(trees)
	for tree in trees:
		labelNodes(tree)
		filter = lambda taxon: True if taxon.label==anch[0] else False
		root = tree.find_node_with_taxon(filter)

		tree.reroot_at_node(root,update_bipartitions=True, suppress_unifurcations=False)
		filter = lambda taxon: True if taxon.label==anch[1] else False
		Q = buildEmtyQuartetTable(anch,taxa,n) 
		node = tree.find_node_with_taxon(filter)
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
						Q  =listQ(Q,leafset,anch)
						
				
			else:
				for ct in ch:
					if ct == nodePre:
						continue
					leafset = ct.leaf_nodes()
					if len(leafset)<2:
						continue
					else:
						Q = listQ(Q,leafset,anch)
	pr.printQuartetsToFile(Q,out+'/tmp.q')
	return Q 
def findAnchoredQuartetTables(D,Q):
	for q in Q:
		q = q.replace("|","")
		key = "/".join(sorted(q.split(" ")))
		D[key] += 1
	for key in D:
		
	return D
		
def buildEmtyQuartetTable(anch,taxa,n):
	taxa = set(taxa)-set(anch)
	Q = dict()
	taxaList = list(itertools.combinations(taxa,2))
	for taxaPair in taxaList:
		taxaPair = sorted(list(taxaPair))
		l = anch[0]+" "+anch[1]+" | "+taxaPair[0]+" "+taxaPair[1]
		Q[l] = [0.5, n]
	return Q
		
		
		
def listQ(leafset,anch):
	for l in list(itertools.combinations(leafset,2)):
		lt = sorted([l[0].label,l[1].label])
		Q[anch[0]+' '+anch[1]+' | '+lt[0]+' '+lt[1]][0] += 1
	return Q	
def findAnchoredDistanceTable(achs,gt,taxa,out):
	Q = findAnchoredQuartets(achs,gt,taxa,out)
	n = len(taxa)
	Q = buildEmtyQuartetTable(achs,taxa,n)
	D = anchoredDistanceFromFrq(Q,achs)
	return D
