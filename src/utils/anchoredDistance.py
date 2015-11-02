#!/usr/bin/env python
import dendropy
import itertools
import sys
import os
from readQuartetTable import readQuartetTable
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
	if(readFromTable):
		frq=readQuartetTable(qfile)	
	else:
		frq=findAnchoredQuartetTable(achs,gt)	
	D = anchoredDistanceFromFrq(frq,achs)
	return D
def anchoredDistanceFromFrq(frq,achs):
	D = dict()
	for k,v in frq.iteritems():
		kt = sorted(k.split("/"))
		if ((achs[0] in kt ) and( achs[1] in kt)):
			s = sorted(list(set(kt)-{achs[0],achs[1]}))
			if s[0]<achs[0]:
				D[s[0]+' '+s[1]]=-np.log(frq[kt][s[1]])
				D[s[1]+' '+s[0]]=D[s[0]+' '+s[1]]
			elif s[0]>achs[0]:
				D[s[0]+' '+s[1]]=-np.log(frq[kt][achs[1]])	
				D[s[1]+' '+s[0]]=D[s[0]+' '+s[1]]
	return D		 
			
