import tempfile
import sys
def printDistanceTable (mapDict , keyDict):
	l = len(keyDict)
	norm = 1.
	ps = 20.
	print l
	for i in range(0,l):
        	sp = keyDict[i]
        	print sp,
        	for j in range(0,l):
                	if i==j:
				b = 0.
                        	print '%0.6f' % b,
                        	continue
                	k = sorted([keyDict[j],keyDict[i]])
			a=mapDict[k[0]+' '+k[1]]
			a += ps
			a /=norm
                	print '%0.6f' % a,
        	print
def printDistanceTableToFile (mapDict , keyDict, outfile):
		
	f1 = open(outfile, 'w')
	keyDict = list(keyDict)
	l = len(keyDict)
	norm = 1.
	ps = 0.
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
	f1.close()
def printQuartetTable(frq):
	for key in frq:
		print key,
		keyt = key.split("/")
		sz = sum(frq[key].values())
		for i in range(1,4):
			print frq[key][keyt[i]]/float(sz),
		print
	return
def printQuartetTableToFile(frq,filename):
	orig_stdout = sys.stdout
        f = open(filename, 'w+')
	sys.stdout = f
	for key in frq:
		print key,
		keyt = key.split("/")
		sz = sum(frq[key].values())
		for i in range(1,4):
			print frq[key][keyt[i]]/float(sz),
		print
	
	sys.stdout = orig_stdout
	f.close()
	return
def printQuartetsToFile(q,filename):
	orig_stdout = sys.stdout
        f = open(filename, 'w+')
	sys.stdout = f
	for l in q:
		print q
	sys.stdout = orig_stdout
	f.close()
	return
