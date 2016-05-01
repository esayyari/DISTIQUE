#!/usr/bin/env python

from checkBranchIfAvail import checkBranchIfAvail;
from uniquePoolOfBranches import uniquePoolOfBranches
import sys;
import os 
from optparse import OptionParser
import itertools
import random
import re
if __name__ == '__main__':
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option("-i","--input",dest="inpt",type="string",
			help="read pool of branches from FILENAME")
	parser.add_option("-s","--species",dest="sp",type="string",
			help="read species info from FILENAME")
	parser.add_option("-o","--output",dest="out",type="string",
			help="the PATH to write the generated files")
	parser.add_option("-v","--verbose",dest="verbose",
			help="Verbose",default=1)
	(options,args) = parser.parse_args()
      	inpt = options.inpt
	out = options.out
	sp = options.sp
	if ( not options.inpt  or not options.out or not options.sp):
		sys.exit("Please enter pool of brnaches file, species info file, and output folder location")
	pool=uniquePoolOfBranches(inpt);
	ppInfo =checkBranchIfAvail(pool,sp)
	f=open(out,'w')
	for x in ppInfo:
		print>> f,x.replace("\n","")
	f.close()
