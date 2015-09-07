import numpy as np
import sys

def prodDistance(frq,ps_count):
	mapDict = dict()
	v = set([0, 1, 2, 3])
	sz = sum(frq[frq.keys()[0]].values())
	for k in sorted(frq.keys()):
		d = sorted(k.split('/'))
		for i in range(1,4):
			distKey = d[0]+' '+d[i]

			if distKey not in mapDict.keys() and frq[k][d[i]]>0 and frq[k][d[i]]!=sz:
				mapDict[distKey] = -np.log(float(frq[k][d[i]])/sz)
			elif distKey not in mapDict.keys() and frq[k][d[i]]<1:
				mapDict[distKey] = -np.log(ps_count)
			elif distKey not in mapDict.keys() and frq[k][d[i]]==sz:
				mapDict[distKey] = -np.log(float(frq[k][d[i]]-ps_count)/sz)
			elif frq[k][d[i]] < 1:
				mapDict[distKey] -= np.log(ps_count)
			elif frq[k][d[i]]==sz:
				mapDict[distKey] -= np.log(float(frq[k][d[i]]-ps_count)/sz)
			else:
				mapDict[distKey] -=  np.log(float(frq[k][d[i]])/sz)
			s = set([0,i])
			g = sorted(list(v - s))
			distKey = d[g[0]]+ ' ' +d[g[1]]

			if distKey not in mapDict.keys() and frq[k][d[i]]>0 and frq[k][d[i]]!=sz:
				mapDict[distKey] = -np.log(float(frq[k][d[i]])/sz)
			elif distKey not in mapDict.keys() and frq[k][d[i]]<1:
				mapDict[distKey] = -np.log(ps_count)
			elif distKey not in mapDict.keys() and frq[k][d[i]]==sz:
				mapDict[distKey] = -np.log(float(frq[k][d[i]]-ps_count)/sz)
			elif frq[k][d[i]] < 1:
				mapDict[distKey] -= np.log(ps_count)
			elif frq[k][d[i]]==sz:
				mapDict[distKey] -= np.log(float(frq[k][d[i]]-ps_count)/sz)
			else:
				mapDict[distKey] -=  np.log(float(frq[k][d[i]])/sz)
	return mapDict
