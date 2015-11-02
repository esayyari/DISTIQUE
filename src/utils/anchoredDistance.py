#!/usr/bin/env python
import dendropy
import itertools
import sys
import os
from readQuartetTable import readQuartetTable
from printDistanceTableToFile import printDistanceTableToFile 


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
			
