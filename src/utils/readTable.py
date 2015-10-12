#!/usr/bin/env python

import sys
import os

def readTable(tmpPath):
	frq = dict()
	out_path = src_fpath = os.path.expanduser(os.path.expandvars(tmpPath))
	f = open(out_path, 'r')
	frq = dict()
	for line in f:
		k=line.split()
		v = dict()
	        d = k[0].split('/')
	        v[d[1]] = float(k[1])
	        v[d[2]] = float(k[2])
	        v[d[3]] = float(k[3])
	   	frq[k[0]] = v
	return frq

