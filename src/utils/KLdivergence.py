#!/usr/bin/env python
import numpy as np
import itertools
def KLdivergence(qTable1,qTable2):
	kld = dict()
	for key in qTable1:
		v1 = qTable1[key]
		sz1 = sum(v1.values())
		v2 = qTable2[key]
		sz2 = sum(v2.values())
		kld[key] = 0
		for key2 in v1:
			kld[key] += float(v1[key2])/float(sz1)*np.log((float(v1[key2])*sz2)/(sz1*float(v2[key2])))
	totalKL = sum(kld.values())
	return totalKL
