#!/usr/bin/env python

def printQuartetTable(frq):
	for key in frq:
		print key,
		keyt = key.split("/")
		for i in range(1,4):
			print frq[key][keyt[i]],
		print
	return
