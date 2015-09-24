import numpy as np
import sys

def minDistance(frq):
	mapDict = dict()
        v = set([0, 1, 2, 3])
        for k in sorted(frq.keys()):
                d = sorted(k.split('/'))
        	sz = sum(frq[k].values())
			
                for i in range(1,4):
                        distKey = d[0]+' '+d[i]

                        if distKey not in mapDict.keys()  :
				mapDict[distKey] = list()
				p = np.log(min(1.,3.*frq[k][d[i]]/(1.*sz)))
                                mapDict[distKey].append(p)
                        else:
				p = np.log(min(1.,3.*frq[k][d[i]]/(1.*sz)))
                                mapDict[distKey].append(p)
                        s = set([0,i])
                        g = sorted(list(v - s))
                        distKey = d[g[0]]+ ' ' +d[g[1]]

                        if distKey not in mapDict.keys():
				mapDict[distKey] = list()
				p = np.log(min(1.,3.*frq[k][d[i]]/(1.*sz)))
                                mapDict[distKey].append(p)
                        else:
				p = np.log(min(1.,3.*frq[k][d[i]]/(1.*sz)))
                                mapDict[distKey].append(p)
			retDict = dict()

	for key in mapDict:
		m = min(mapDict[key])
		if m==np.log(1.):
			m = 0
		else: 	
			m = min(mapDict[key])
		retDict[key] = -m
        return retDict
	
