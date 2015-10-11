#!/usr/bin/env python

import sys
import itertools
import dendropy
import random
def generateKey(taxa_list):
	chosen = list()
	for k, v in taxa_list.iteritems():
		rt = random.sample(v,1)
		chosen.append(rt[0])
	allQuartetComb = itertools.combinations(chosen,4)
	origKeys = ['/'.join(sorted(q)) for q in allQuartetComb ]
#	print len(origKeys)
	return origKeys
		
		
		
