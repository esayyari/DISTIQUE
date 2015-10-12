import numpy as np
import sys

def prodDistance(frq):
	mapDict = dict()
	v = set([0, 1, 2, 3])
	for k in sorted(frq.keys()):
		d = sorted(k.split('/'))
		
		sz = sum(frq[k].values())
		for i in range(1,4):
			distKey = d[0]+' '+d[i]

			if distKey not in set(mapDict.keys()) and frq[k][d[i]]!=sz:
				mapDict[distKey] = -np.log(float(frq[k][d[i]])/sz)
			elif distKey not in set(mapDict.keys()) and frq[k][d[i]]==sz:
				mapDict[distKey] = -np.log(float(frq[k][d[i]]-(10.**(-ps_count)))/sz)
			elif frq[k][d[i]]==sz:
				mapDict[distKey] -= np.log(float(frq[k][d[i]]-(10.**(-ps_count)))/sz)
			else:
				mapDict[distKey] -=  np.log(float(frq[k][d[i]])/sz)
			s = set([0,i])
			g = sorted(list(v - s))
			distKey = d[g[0]]+ ' ' +d[g[1]]

			if distKey not in set(mapDict.keys())  and frq[k][d[i]]!=sz:
				mapDict[distKey] = -np.log(float(frq[k][d[i]])/sz)
			elif distKey not in set(mapDict.keys()) and frq[k][d[i]]==sz:
				mapDict[distKey] = -np.log(float(frq[k][d[i]]-(10.**(-ps_count)))/sz)
			elif frq[k][d[i]]==sz:
				mapDict[distKey] -= np.log(float(frq[k][d[i]]-(10.**(-ps_count)))/sz)
			else:
				mapDict[distKey] -=  np.log(float(frq[k][d[i]])/sz)
	return mapDict
