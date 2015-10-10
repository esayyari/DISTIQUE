#!/usr/bin/env python

import sys
import os
import numpy as np
import itertools
import dendropy
from collections import defaultdict
import random
def getTaxaList(taxaDict):
        taxa = dict()
        inv_taxa=dict()
        taxa_list=dict()
        n = len(taxaDict)
        for key1,vals1 in taxaDict.iteritems():
                v = list()
                for t in vals1:
                        k = t.taxon.__str__()
                        kv = k.split("'")
                        v.append(kv[1])
                taxa_list[key1] = v
	for k, v in taxa_list.iteritems():
                for v2 in v:
                        inv_taxa[v2] = k
	return (taxa_list,inv_taxa)
