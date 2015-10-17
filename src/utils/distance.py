import sys
import numpy as np
import os
from optparse import OptionParser
from  prodDistance import prodDistance
from printDistanceTable import printDistanceTable
from minDistance import minDistance
import pprint

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--file", dest="filename", type="string",
	      help="read data from FILENAME")
parser.add_option("-m","--method",dest="method", default="min",choices=["min", "prod", "minavg", "minmed"],
	      help="define the method to compute distance: min "
		   "prod, minavg,minmed [default: prod]")
parser.add_option("-p","--percentile",dest="p",default=1, type="int",
	      help="for the methods minavg, or minmed, "
		   "find  average or median as the true probability in provided percentile")   
parser.add_option("-c","--seudo_count",dest="ps_count",type="int",default=8,
	      help="pseudo count that will be replaced with the zero frequencies. " 
		   "It will be the power of 10 i.e. 1e-ps_count  [default: 8]")
(options,args) = parser.parse_args()
if (options.method == "minavg" or options.method=="minmed") and ~options.p :
	print("warning, the percentile set to its default value (1)")

filename = options.filename
method = options.method
percentile = options.p
pst=-options.ps_count
ps_count = -pst
f = open(filename, 'r')
frq=dict()
for line in f:
    k=line.split()
    v = dict()
    d = k[0].split('/')
    v[d[1]] = float(k[1])
    v[d[2]] = float(k[2])
    v[d[3]] = float(k[3])
    frq[k[0]] = v
keyDict = sorted(np.unique(("/".join(frq.keys())).split("/")));
mapDict = dict()
def computeDistance(frq, method, percentile):
    if method == 'min':
		return minDistance(frq)
    else:
		return prodDistance(frq)
mapDict = computeDistance(frq,method,percentile)
printDistanceTable(mapDict,keyDict)

	
