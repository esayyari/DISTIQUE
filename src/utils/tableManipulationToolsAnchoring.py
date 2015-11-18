#! /usr/bin/env python
import sys
import csv
import numpy as np
from optparse import OptionParser
import dendropy
import itertools
import os
import printTools as pr
from scipy import stats
from numpy import mean, sqrt, square, arange
from  prodDistance import prodDistance
from minDistance import minDistance
import copy
import toolsTreeTaxa as tstt
import numpy as np

	
def generateKey(taxa_list):
	chosen = list()
	for k, v in taxa_list.iteritems():
		rt = random.sample(v,1)
		chosen.append(rt[0])
	allQuartetComb = itertools.combinations(chosen,4)
	origKeys = ['/'.join(sorted(q)) for q in allQuartetComb ]
	return origKeys
	
		
def findTrueAverageTableAnchoring(frq,anch,list_taxa,method):
	anch = sorted(list(anch))
	lst_taxa = list(list_taxa.keys())
	TotalKey = dict()
	n = len(lst_taxa)
	for i in range(0,n):
		for j in range(i+1,n):
			for taxon_i in list_taxa[lst_taxa[i]]:
				for taxon_j in list_taxa[lst_taxa[j]]:
					keyt = "/".join(sorted([taxon_i,taxon_j,anch[0],anch[1]]))
					lab_taxon_i = taxon_i
					lab_taxon_j = taxon_j
					lab_taxon_k = anch[0]
					lab_taxon_z = anch[1]
					key_orig = "/".join(sorted([lab_taxon_i,lab_taxon_j,lab_taxon_k,lab_taxon_z]))
					
					l = sorted([lst_taxa[i],lst_taxa[j],anch[0],anch[1]])
					key_inv = "/".join(l)
					v = frq[key_orig]
					v_inv = v[0]/v[1]	
					if key_inv in TotalKey:
						vt = TotalKey[key_inv] 
						vt.append(v_inv)
					else:
						vt = list()
						vt.append(v_inv)
					TotalKey[key_inv] = vt
									
	TotalKeyf = dict()
	for q,v in TotalKey.iteritems():
		if method == "gmean":
			vtt = stats.gmean(v)
		elif method == "mean":
			vtt = mean(v)
		else:
			vtt = sqrt(mean(square(v)))
		TotalKeyf[q] = vtt
	return TotalKeyf
