#!/usr/bin/env python
import dendropy
import sys
import os

from getQuartetKeys import generateKey
import itertools

def partialQuartetTable(quartTable,origKeys,inv_taxa):	
	pQuartTable = dict()
	inv_keys = dict()
	for key in origKeys:
		tmp = key.split("/")
		inv_tmp = list()
		#print key,quartTable[key][tmp[1]],quartTable[key][tmp[2]], quartTable[key][tmp[3]]
		for k in tmp:
		#	print k+":"+inv_taxa[k]
			it = inv_taxa[k]
			inv_tmp.append(it)
			inv_keys[it] = k
		inv_tmp_sorted = sorted(inv_tmp)
		dummyKey = "/".join(inv_tmp_sorted)
		t = dict()
		for i in range(1,4):
			q1 = sorted([inv_keys[inv_tmp_sorted[0]], inv_keys[inv_tmp_sorted[i]]] )
			s  = {1,2,3}-{i}
			s = list(s)
			q2 = sorted([inv_keys[inv_tmp_sorted[s[0]]], inv_keys[inv_tmp_sorted[s[1]]]])
			if q2[0]<q1[0]:
				t[inv_tmp_sorted[i]] = quartTable[key][q2[1]]
		#		print q2[1]+":"+inv_tmp_sorted[0]+inv_tmp_sorted[i]
			else:
		#		print q1[1]+":"+inv_tmp_sorted[0]+inv_tmp_sorted[i]
				t[inv_tmp_sorted[i]] = quartTable[key][q1[1]]
		pQuartTable[dummyKey] = t
#	for k,v in pQuartTable.iteritems():
		#print k,
#		a = sorted(k.split("/"))
#		for i in range(1,4):
		#	print v[a[i]]/sum(v.values()),
		#print	
	return pQuartTable
