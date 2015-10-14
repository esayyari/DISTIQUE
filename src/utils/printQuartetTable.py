#!/usr/bin/env python

def printQuartetTable(frq):
	for key in frq:
		print key,
		keyt = key.split("/")
		sz = sum(frq[key].values())
		for i in range(1,4):
			print frq[key][keyt[i]]/float(sz),
		print
	return
