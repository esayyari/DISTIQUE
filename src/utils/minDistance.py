import numpy as np
import sys

def minDistance(frq,ps_count):
	mapDict = dict()
        v = set([0, 1, 2, 3])
        for k in sorted(frq.keys()):
                d = sorted(k.split('/'))
		
        	sz = sum(frq[k].values())
                for i in range(1,4):
                        distKey = d[0]+' '+d[i]

                        if distKey not in mapDict.keys() and frq[k][d[i]]>0 :
				mapDict[distKey] = list()
				p = min(1.5,3.*frq[k][d[i]]/float(sz))
                                mapDict[distKey].append(p)
                        elif distKey not in mapDict.keys() and frq[k][d[i]]<1:
				mapDict[distKey] = list()
                                mapDict[distKey].append(-ps_count)
                        elif frq[k][d[i]] < 1:
                                mapDict[distKey].append(-ps_count)
                        else:
				p = np.log(min(1.5,3.*frq[k][d[i]]/float(sz)))
                                mapDict[distKey].append(p)
                        s = set([0,i])
                        g = sorted(list(v - s))
                        distKey = d[g[0]]+ ' ' +d[g[1]]


                        if distKey not in mapDict.keys() and frq[k][d[i]]>0:
				mapDict[distKey] = list()
				p = min(1.5,3.*np.log(frq[k][d[i]])/float(sz))
                                mapDict[distKey].append(p)
                        elif distKey not in mapDict.keys() and frq[k][d[i]]<1:
				mapDict[distKey] = list()
                                mapDict[distKey].append(-ps_count)
                        elif frq[k][d[i]] < 1:
                                mapDict[distKey].append(-ps_count)
                        else:
				p = np.log(min(1.5,3.*frq[k][d[i]]/float(sz)))
                                mapDict[distKey].append(p)
	retDict = dict()

	for key in mapDict:
		m = min(mapDict[key])
#		print key,
#		for i in mapDict[key]:
#			print i,
#		print
		if m==np.log(1.5):
			m = 0
		else: 	
			m = min(mapDict[key])
		retDict[key] = -m
        return retDict
	
