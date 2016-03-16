import numpy as np
from numpy import Inf

def prodDistance(frq,met):
	mapDict = dict()
	v = set([0, 1, 2, 3])
	ps_count = 5
	for k in sorted(frq.keys()):
		d = sorted(k.split('/'))
		
		sz = sum(frq[k].values())
		for i in range(1,4):
			distKey = d[0]+' '+d[i]
				
			if (met == "freq" or met == "log"):
				if distKey not in set(mapDict.keys()) and frq[k][d[i]]!=sz:
					if frq[k][d[i]]/sz <= 0:
							print k,frq[k],"It is zero!"
					elif frq[k][d[i]] == Inf:
							print k, frq[k][d[i]],"It is infinity"
					mapDict[distKey] = -np.log(float(frq[k][d[i]])/sz)
				elif distKey not in set(mapDict.keys()) and frq[k][d[i]]==sz:
					if frq[k][d[i]]/sz <= 0:
							print k,frq[k],"It is zero!"
					elif frq[k][d[i]] == Inf:
							print k, frq[k][d[i]],"It is infinity"
					mapDict[distKey] = -np.log(float(frq[k][d[i]]-(10.**(-ps_count)))/sz)
				elif frq[k][d[i]]==sz:
					if frq[k][d[i]]/sz <= 0:
							print k,frq[k],"It is zero!"
					elif frq[k][d[i]] == Inf:
							print k, frq[k][d[i]],"It is infinity"
					mapDict[distKey] -= np.log(float(frq[k][d[i]]-(10.**(-ps_count)))/sz)
				else:
					if frq[k][d[i]]/sz <= 0:
							print k,frq[k],"It is zero!"
					elif frq[k][d[i]] == Inf:
							print k, frq[k][d[i]],"It is infinity"
					mapDict[distKey] -=  np.log(float(frq[k][d[i]])/sz)
				s = set([0,i])
				g = sorted(list(v - s))
				distKey = d[g[0]]+ ' ' +d[g[1]]

				if distKey not in set(mapDict.keys())  and frq[k][d[i]]!=sz:
					if frq[k][d[i]]/sz <= 0:
							print k,frq[k],"It is zero!"
					elif frq[k][d[i]] == Inf:
							print k, frq[k][d[i]],"It is infinity"
					mapDict[distKey] = -np.log(float(frq[k][d[i]])/sz)
				elif distKey not in set(mapDict.keys()) and frq[k][d[i]]==sz:
					if frq[k][d[i]]/sz <= 0:
							print k,frq[k],"It is zero!"
					elif frq[k][d[i]] == Inf:
							print k, frq[k][d[i]],"It is infinity"
					mapDict[distKey] = -np.log(float(frq[k][d[i]]-(10.**(-ps_count)))/sz)
				elif frq[k][d[i]]==sz:
					if frq[k][d[i]]/sz <= 0:
							print k,frq[k],"It is zero!"
					elif frq[k][d[i]] == Inf:
							print k, frq[k][d[i]],"It is infinity"
					mapDict[distKey] -= np.log(float(frq[k][d[i]]-(10.**(-ps_count)))/sz)
				else:
					if frq[k][d[i]]/sz <= 0:
							print k,frq[k],"It is zero!"
					elif frq[k][d[i]] == Inf:
							print k, frq[k][d[i]],"It is infinity"
					mapDict[distKey] -=  np.log(float(frq[k][d[i]])/sz)
	return mapDict
