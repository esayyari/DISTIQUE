#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import numpy as np
import itertools
import random
import subprocess
import re
def uniquePoolOfBranches(inpt):
	poolBranches = dict()
        listPoolBranches = list()
        f = open(inpt, 'r')
        for x in f:
                y = re.search('^\{',x)
                if y:
                        qtInfo=re.sub('\[.*','',x)
                        qtInfo=re.sub(' ','',qtInfo)
                        qtInfo= re.sub("\n","",qtInfo)
                        if qtInfo in poolBranches:
                                continue
                        else:
                                bpInfo=re.sub('.*\[','[',x)
                                bpInfo=re.sub('[\[\]{}]','',bpInfo)
                                poolBranches[qtInfo] = 1
                                y = x.replace("\n","")
                                listPoolBranches.append(y)
        f.close()
	return listPoolBranches
if __name__ == '__main__':
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option("-i","--input",dest="inpt",type="string",
			help="read input branches from FILENAME")
	parser.add_option("-o","--output",dest="out",type="string",
			help="the PATH to write the generated files")
	parser.add_option("-v","--verbose",dest="verbose",
			help="Verbose",default=1)
	(options,args) = parser.parse_args()
	inpt = options.inpt
	outpath = options.out
	if ( not options.inpt  or not options.out):
		sys.exit("Please enter pool of brnaches file, and output folder location")

	poolBranches = dict()
	listPoolBranches = list()
	f = open(inpt, 'r')
	n = 0
	m = 0
	c = 0
	for x in f:
		y = re.search('^\{',x)
		if y:
			n += 1
			qtInfo=re.sub('\[.*','',x)
			qtInfo=re.sub(' ','',qtInfo)
			qtInfo= re.sub("\n","",qtInfo)
			if qtInfo in poolBranches:
				
				m +=1
				continue
			else:
				c += 1
				bpInfo=re.sub('.*\[','[',x)
				bpInfo=re.sub('[\[\]{}]','',bpInfo)
				poolBranches[qtInfo] = 1
				y = x.replace("\n","")
				listPoolBranches.append(y)
	f.close()
	g = open(outpath+'/poolOfBranches.txt','w')
	for item in listPoolBranches:
		print>>g, item.replace("\n","")
	g.close()
	print m
	print c
	print n
