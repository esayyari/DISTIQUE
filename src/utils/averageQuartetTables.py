#!/usr/bin/env python
import dendropy
from getQuartetKeys import generateKey
from findPartialQuartetFreq import partialQuartetTable
from addQuartetTables import addQuartetTables
from checkConvergenceQuartTable import convergencedQuartTable
def averageQuartetTables(eps,numToStop,numMax,taxa_list,frq,taxa_inv,verbose):
	num = 0	
	for a in range(0,numMax):
		origKeys = generateKey(taxa_list)
		partialTable1= partialQuartetTable(frq,origKeys,taxa_inv)
		if a>0:
			partialTable1=addQuartetTables(partialTable1,quartTable)
			if convergencedQuartTable(partialTable1,quartTable,eps,verbose):
				quartTable = partialTable1
				if num == numToStop:
					print "Quartet Table converged with " + str(a)+" steps"
					return quartTable
				else:
					num += 1
			else:
				num = 0
		quartTable = partialTable1
	print "Warning: quit averaging Quartet Tables without convergence"
	return quartTable
