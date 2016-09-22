#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import itertools
import random
import re

def readTable(tmpPath):
        out_path  = os.path.expanduser(os.path.expandvars(tmpPath))
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
	f.close()
        return frq



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
		print b1,b2
		b = (b1,b2)
		lq.append(b)
	f.close()
	return lq
def expandListName(filename):
        src_fpath = os.path.expanduser(os.path.expandvars(filename))
	if not os.path.exists(src_fpath):
		sys.stderr.write('Not found: "%s"' % src_fpath)	
	return src_fpath	
def generateBranchName(c1,c2,sister,remain):
	
	if c1[0]<c2[0]:
		if c1[0]<min(sister[0],remain[0]):
			mainKey = ",".join(sorted(c1+c2)) + "|" + ",".join(sorted(sister+remain))
			if sister[0]<remain[0]:
				tp1 = ",".join(sorted(c1+sister)) + "|" + ",".join(sorted(remain+c2))
				tp2 = ",".join(sorted(c1+remain)) + "|" + ",".join(sorted(sister+c2))
			else:
				tp1 = ",".join(sorted(c1+remain)) + "|" + ",".join(sorted(sister+c2))
				tp2 = ",".join(sorted(c1+sister)) + "|" + ",".join(sorted(remain+c2))
		else:
			mainKey = ",".join(sorted(sister+remain)) + "|" + ",".join(sorted(c1+c2))
			if sister[0]<remain[0]:
				tp1 = ",".join(sorted(sister+c1)) + "|" + ",".join(sorted(remain+c2))
				tp2 = ",".join(sorted(sister+c2)) + "|" + ",".join(sorted(remain+c1))
			else:
				tp1 = ",".join(sorted(remain+c1)) + "|" + ",".join(sorted(sister+c2))
                                tp2 = ",".join(sorted(remain+c2)) + "|" + ",".join(sorted(sister+c1))
	else:
                if c2[0]<min(sister[0],remain[0]):
                        mainKey = ",".join(sorted(c1+c2)) + "|" + ",".join(sorted(sister+remain))
                        if sister[0]<remain[0]:
                                tp1 = ",".join(sorted(c2+sister)) + "|" + ",".join(sorted(remain+c1))
                                tp2 = ",".join(sorted(c2+remain)) + "|" + ",".join(sorted(sister+c1))
                        else:
                                tp1 = ",".join(sorted(c2+remain)) + "|" + ",".join(sorted(sister+c1))
                                tp2 = ",".join(sorted(c2+sister)) + "|" + ",".join(sorted(remain+c1))
                else:
                        mainKey = ",".join(sorted(sister+remain)) + "|" + ",".join(sorted(c1+c2))
                        if sister[0]<remain[0]:
                                tp1 = ",".join(sorted(sister+c2)) + "|" + ",".join(sorted(remain+c1))
                                tp2 = ",".join(sorted(sister+c1)) + "|" + ",".join(sorted(remain+c2))
                        else:
                                tp1 = ",".join(sorted(remain+c2)) + "|" + ",".join(sorted(sister+c1))
                                tp2 = ",".join(sorted(remain+c1)) + "|" + ",".join(sorted(sister+c2))

	return (mainKey,tp1,tp2)	
def readQuartetTable(src_fpath):
        frq = readTable(src_fpath)
        for k in frq:
                sz = sum(frq[k].values())
                for k2 in frq[k]:
                        frq[k][k2] /= sz
        return frq

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
	frq = readQuartetTable(qt)
	frqDist = dict()
	frqAlter = dict()
	for l in lq:
		c1 = sorted(l[0][0])
		c2 = sorted(l[0][1])
		sister = sorted(l[1][0])
		remain = sorted(l[1][1])
		frqTmp = list()
		frqAlterTmp1 = list()	
		frqAlterTmp2 = list()
		(bName,b1Name,b2Name) = generateBranchName(c1,c2,sister,remain)
		print bName,b1Name,b2Name
		a1 = sorted(b1Name.split("|"))
		
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
							k = list(k)[0]
							f = frq[mainKey][k]
							frqTmp.append(f)
							
							a2 = set(a1[0].split(","))-set(c1)-set(c2)
							if tp1 in a2:
								frqAlterTmp1.append(frq[mainKey][tp1])
								frqAlterTmp2.append(frq[mainKey][tp2])
							else:
								frqAlterTmp1.append(frq[mainKey][tp2])
								frqAlterTmp2.append(frq[mainKey][tp1])
						else:
							k = q2 - {sl[0]}
							k = list(k)[0]
							f = frq[mainKey][k]
							frqTmp.append(f)
							a2 = set(a1[0].split(","))-set(sister)-set(remain)
                                                        if t1 in a2:
                                                                frqAlterTmp1.append(frq[mainKey][t1])
                                                                frqAlterTmp2.append(frq[mainKey][t2])
                                                        else:
                                                                frqAlterTmp1.append(frq[mainKey][t2])
                                                                frqAlterTmp2.append(frq[mainKey][t1])
		frqDist[bName] = frqTmp
		frqAlter[bName] = dict()
		frqAlter[bName][b1Name] = frqAlterTmp1
		frqAlter[bName][b2Name] = frqAlterTmp2
	keyL = sorted(frqDist.keys())
	dictMKeys = dict()
	dictAlTop = dict()
	for i in range(0,len(keyL)):
		dictMKeys[keyL[i]] = "Branch-"+str(i)
	tmp = list()
	for key in frqAlter:
		keyt = sorted(frqAlter[key].keys())
		dictAlTop[keyt[0]] = "top-1"
		dictAlTop[keyt[1]] = "top-2"
		
	f = open(output,'w')
	for bName in frqDist:
		for val in frqDist[bName]:
			print >>f, dictMKeys[bName],"top-0:", str(val)
		for nt in frqAlter[bName]:
			for val in frqAlter[bName][nt]:
				print >> f, dictMKeys[bName], dictAlTop[nt]+":",str(val)
	f.close()
