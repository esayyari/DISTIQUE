#!/usr/bin/env python

import sys
import os
import numpy as np
import itertools
import dendropy
from collections import defaultdict
import random
def generateKey(taxaDict,taxa_list):
	chosen = list()
	for k, v in taxa_list.iteritems():
		rt = random.sample(v,1)
		chosen.append(rt[0])
	allQuartetComb = itertools.combinations(chosen,4)
	origKeys = ['/'.join(sorted(q)) for q in allQuartetComb ]
	return origKeys
		
		
		
