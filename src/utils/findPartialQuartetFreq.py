#!/usr/bin/env python
import dendropy
import sys
import os

from getQuartetKeys import generateKey
import itertools
from expandDict import expandDict 
def partialQuartetTable(quartTable,origKeys,inv_taxa):	
	expandedQTable = expandDict(quartTable)
	pQuartTable = dict()
	inv_keys = dict()
	for key in origKeys:
		tmp = key.split("/")
		inv_tmp = list()
		for k in tmp:
#			print k+":"+inv_taxa[k]
			it = inv_taxa[k]
			inv_tmp.append(it)
			inv_keys[it] = k
		inv_tmp_sorted = sorted(inv_tmp)
		dummyKey = "/".join(inv_tmp_sorted)
		t = dict()
		for i in range(1,4):
			t[inv_tmp_sorted[i]] = expandedQTable[key][inv_keys[inv_tmp_sorted[0]]+inv_keys[inv_tmp_sorted[i]]]	
		pQuartTable[dummyKey] = t
#	for k,v in pQuartTable.iteritems():
#		print k,
#		a = sorted(k.split("/"))
#		for i in range(1,4):
#			print v[a[i]]/sum(v.values()),
#		print	
	return pQuartTable
