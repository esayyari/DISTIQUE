#!/usr/bin/env python
import dendropy
from getQuartetKeys import generateKey
from findPartialQuartetFreq import partialQuartetTable
from addQuartetTables import addQuartetTables
from checkConvergenceQuartTable import convergencedQuartTable
import itertools
from readQuartetTable import readQuartetTable
from findQuartetTable import findQuartetTable
from findTrueAverageTable import findTrueAverageTable
def averageQuartetTables( **kwargs):
	for k,v in kwargs.iteritems():
		if k == 'limit':
			eps = v
		elif k == 'NumToStop':
			numToStop = v
		elif k == 'NumMax':
			numMax = v
		if k == 'ListTaxa':
			taxa_list = v
		elif k == 'QTable':
			frq = v
		elif k == 'QtablePath':
			filename = v
		elif k == 'QtableReady':
			availTable = v
		elif k == 'Inv':
			taxa_inv = v
		elif k == 'V':
			verbose = v
		elif k == 'workingPath':
			wrkPath = v
		elif k == 'treeList':
			trees = v
		elif k == 'KeyType':
			keyType = v
	if availTable and ('QTable' not in kwargs.keys()):
		frq = readQuartetTable(QtablePath)
			
	
	num = 0	
	for a in range(0,numMax):
		origKeys = generateKey(taxa_list)
		if availTable:
			partialTable1= partialQuartetTable(frq,origKeys,taxa_inv)
		else:
			frq = findQuartetTable(trees,origKeys,keyType,wrkPath,verbose)
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
	print "Partial quartet table computed"
	return quartTable
