#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import itertools
import random
import re


def parseQuartetList(filename):
	f = open(filename,'r')
	lq = list()
	for line in f:
		spl = line.split("|")
		b1 = list()
		b2 = list()
		bt1 = spl[0].split("}{")
		
		bt2 = spl[1].split("}{")
		for c in bt1:
			c = c.replace("}","")
			c = c.replace("{","")
			c = c.replace("\n","")
			c = c.replace(" ","")
			b1.append(c.split(","))
		for c in bt2:
			c = c.replace("{","")
			c = c.replace("}","")
			c = c.replace("\n","")
			c = c.replace(" ","")
			b2.append(c.split(","))
		b = (b1,b2)
		lq.append(b)
	return lq
def expandListName(filename):
        src_fpath = os.path.expanduser(os.path.expandvars(filename))
	if not os.path.exists(src_fpath):
		sys.stderr.write('Not found: "%s"' % src_fpath)	
	return src_fpath	
def generateBranchName(c1,c2,sister,remain):
	c = ",".join(sorted(c1 + c2))
	sis = ",".join(sorted(sister+remain))
	if c<sis:
		branchName = c+"|"+sis
	else:
		branchName = sis+"|"+c
	return branchName	

if __name__ == '__main__':
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option("-i", "--file", dest="filename", type="string",
		      help="read list of quartet around internal branches from FILENAME")
	parser.add_option("-q","--quartet",dest="qt",type="string",
			help="read quartet frequencies from FILENAME")
	parser.add_option("-o","--output",dest="out",type="string",
			help="the PATH to write the generated files")
	(options,args) = parser.parse_args()
	if ( not options.filename  or not options.qt):
		sys.exit("Please enter quartet frequency file, and list of quartet around internal branches filename")
	filename = options.filename
	output = options.out
	qt = options.qt
	qt = expandListName(qt)
	filename = expandListName(filename)
	lq = parseQuartetList(filename)
	frqDist = list()
	frqAlter = list()
	for l in lq:
		c1 = sorted(l[0][0])
		c2 = sorted(l[0][1])
		sister = sorted(l[1][0])
		remain = sorted(l[1][1])
		frqTmp = list()
		frqAlterTmp = list()	
		
		for t1 in c1:
			for t2 in c2:
				for tp1 in sister:
					for tp2 in remain:
						sl = sorted([t1,t2,tp1,tp2])
						mainKey = "/".join(sl)
						q1 = {t1,t2}
						q2 = {tp1,tp2}

						if sl[0] in q1:
							k = q1 - {sl[0]}
							f = frq[key][k]
							frqTmp.append(f)
							frqAlterTmp.append(frq[key][tp1])
							frqAlterTmp.append(frq[key][tp2])
						else:
							k = q2 - {sl[0]}
							f = frq[key][k]
							frqTmp.append(f)
							frqAlterTmp.append(frq[key][t1])
							frqAlterTmp.append(frq[key][t2])
		bName = generateBranchName(c1,c2,sister,remain)
		frqDist.append(frqTmp)
		frqAlterTmp.append(frqAlterTmp)
	
