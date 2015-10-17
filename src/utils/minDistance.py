import numpy as np
import sys

def minDistance(frq):
	mapDict = dict()
	numPc = dict()
        v = set([0, 1, 2, 3])
        for k in sorted(frq.keys()):
                d = sorted(k.split('/'))
        	sz = sum(frq[k].values())
			
                for i in range(1,4):
                        distKey = d[0]+' '+d[i]

                        if distKey not in set(mapDict.keys())  :
				mapDict[distKey] = list()
				numPc[distKey] = 0
				if np.abs(frq[k][d[i]]) < 1: 
					numPc[distKey] += 1
				p = np.log(min(1.,3.*frq[k][d[i]]/(1.*sz)))
                                mapDict[distKey].append(p)
                        else:
				if np.abs(frq[k][d[i]]) < 1:
					numPc[distKey] += 1
				p = np.log(min(1.,3.*frq[k][d[i]]/(1.*sz)))
                                mapDict[distKey].append(p)
                        s = set([0,i])
                        g = sorted(list(v - s))
                        distKey = d[g[0]]+ ' ' +d[g[1]]

                        if distKey not in set(mapDict.keys()):
				numPc[distKey] = 0
                                if np.abs(frq[k][d[i]]) < 1: 
                                        numPc[distKey] += 1
				mapDict[distKey] = list()
				p = np.log(min(1.,3.*frq[k][d[i]]/(1.*sz)))
                                mapDict[distKey].append(p)
                        else:
				p = np.log(min(1.,3.*frq[k][d[i]]/(1.*sz)))
				if np.abs(frq[k][d[i]]) < 1:
                                        numPc[distKey] += 1
                                mapDict[distKey].append(p)
			retDict = dict()

	for key in mapDict:
#		ml = sorted(mapDict[key])
#		m = np.mean( ml[0:int(sz/100.*2)+1] )
		m = min(mapDict[key])
		pc = numPc[key]
		if m==np.log(1):
			m = 0
		elif pc>0:
			m = pc*np.log(0.5/sz)
		else: 	
			m = min(mapDict[key])
		retDict[key] = -m
        return retDict
	
