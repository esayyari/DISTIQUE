def printDistanceTableToFile (mapDict , keyDict, outfile):
	f1 = open(outfile, 'w+')
	l = len(keyDict)
	norm = 1.
	ps = 20.
	print >> f1, l
	for i in range(0,l):
        	sp = keyDict[i]
        	print >> f1, sp,
        	for j in range(0,l):
                	if i==j:
				b = 0.
                        	print >> f1, '%0.6f' % b,
                        	continue
                	k = sorted([keyDict[j],keyDict[i]])
			a=mapDict[k[0]+' '+k[1]]
			a += ps
			a /=norm
                	print >> f1,'%0.6f' % a,
        	print >> f1 
