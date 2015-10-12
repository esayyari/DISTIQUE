#!/usr/bin/env python

from KLdivergence import KLdivergence 
import numpy as np

def convergencedQuartTable(qTable1,qTable2,eps,verbose):
	kl = KLdivergence(qTable1,qTable2)
	if kl<eps:
		if verbose==1:
			print kl
		return True
	else:
		return False
	

