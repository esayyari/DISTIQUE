#! /usr/bin/env python
import sys
import csv
import numpy as np
import pprint
from optparse import OptionParser
     
files=map(str, sys.argv[1:])

table=[]
i=0
for f in files:
    if f is None: continue
    infile=open(f,'r')
    table.append(list())
    for row in csv.reader(infile,delimiter=' '):
	if row=="None": continue
        table[i].append(row)
    infile.close()
    i +=1

tmp=[]
tablef=[]
for z in range(0,len(table)):
    tmp.append(list())
    for r in range(1,len(table[z])):
        for c in range(1,len(table[z][1])):
            table[z][r][c]=float(table[z][r][c])
            tmp[z].append(table[z][r][c])
    tmp[z]= sorted(tmp[z],reverse=True)
for r in range(1,len(table[0])):
    for c in range(1,len(table[0][1])):
        table[0][r][c] /= tmp[0][0]
for z in range(1,len(table)):
    for r in range(1,len(table[z])):
        for c in range(1,len(table[z][1])):
            table[z][r][c] /= tmp[z][0]
            table[0][r][c] += table[z][r][c]
print table[0][0][0]
for r in range(1,len(table[z])):
    print table[0][r][0],
    for c in range(1,len(table[z][1])):
        print '%0.6f' % table[0][r][c],
    print

