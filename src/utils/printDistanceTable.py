def printDistanceTable (mapDict , keyDict):
	l = len(keyDict)
	print l
	for i in range(0,l):
        	sp = keyDict[i]
        	print sp,
        	for j in range(0,l):
                	if i==j:
                        	print '%0.6f' % 0,
                        	continue
                	k = sorted([keyDict[j],keyDict[i]])

                	print '%0.6f' % mapDict[k[0]+' '+k[1]],
        	print
