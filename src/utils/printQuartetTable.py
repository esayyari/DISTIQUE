#!/usr/bin/env python
import sys
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
	f = file(filename, 'w')
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
