#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import itertools
import random
import re
def checkBranchIfAvail(g,sp):
	
	poolSpeciesBranches = dict()

	branches = dict()
	f = open(sp, 'r')
	for x in f:
		y = re.search('^\{',x)
		if y:
			bpInfo=re.sub('.*\[','[',x)
			bpInfo=re.sub('[\[\]{}]','',bpInfo)
			bpInfo=re.sub(' :.*','',bpInfo)
			bpInfo=re.sub("\n","",bpInfo)
			poolSpeciesBranches[bpInfo] = 1;
	f.close()
	print len(poolSpeciesBranches)
	ct = 0;
	branchList = list()
	for line in g:
		y = re.search(':',line)
		if y:
			qtInfo=re.sub('\[.*','',line)
			qtInfo=re.sub(' ','@',qtInfo)
			bpInfo=re.sub('.*\[','[',line)
			bpInfo=re.sub('[\[\]{}]','',bpInfo)
			bpInfo=re.sub(' :.*','',bpInfo)
			pp     = re.sub('.*:','',line)
			if bpInfo in poolSpeciesBranches:
				branchList.append(str(1)+" "+str(pp) +" "+ str(ct/3)+" "+str(ct%3));
				ct += 1
			else:
				branchList.append(str(0)+" "+str(pp) + " "+str(ct/3)+" "+str(ct%3));
				ct += 1;

	return branchList;
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
	outpath = options.out
	sp = options.sp
	if ( not options.inpt  or not options.out or not options.sp):
		sys.exit("Please enter pool of brnaches file, species info file, and output folder location")

	poolSpeciesBranches = dict()

	branches = dict()
	f = open(sp, 'r')
	for x in f:
		y = re.search('^\{',x)
		if y:
			bpInfo=re.sub('.*\[','[',x)
			bpInfo=re.sub('[\[\]{}]','',bpInfo)
			bpInfo=re.sub(' :.*','',bpInfo)
			poolSpeciesBranches[bpInfo] = 1;
	f.close()
	print len(poolSpeciesBranches)
	g = open(inpt,'r')
	ct = 0;
	branchList = list()
	for line in g:
		y = re.search(':',line)
		if y:
			qtInfo=re.sub('\[.*','',line)
			qtInfo=re.sub(' ','@',qtInfo)
			bpInfo=re.sub('.*\[','[',line)
			bpInfo=re.sub('[\[\]{}]','',bpInfo)
			bpInfo=re.sub(' :.*','',bpInfo)
			pp     = re.sub('.*:','',line)
			if bpInfo in poolSpeciesBranches:
				branchList.append(str(1)+" "+str(pp) +" "+ str(ct/3)+" "+str(ct%3));
				ct += 1
			else:
				branchList.append(str(0)+" "+str(pp) + " "+str(ct/3)+" "+str(ct%3));
				ct += 1;
	g.close()
	f = open(outpath+"/ppOfBranches.txt",'w')
	for item in branchList:
		print >> f,item.replace("\n","")
	f.close()
