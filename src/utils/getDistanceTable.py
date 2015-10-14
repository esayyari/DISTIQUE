#!/usr/bin/env python

import sys
import numpy as np
import os
from optparse import OptionParser
from  prodDistance import prodDistance
from printDistanceTableToFile import printDistanceTableToFile
from minDistance import minDistance
import pprint


def distanceTable(frq,method,outfile):
	percentile = 1
	keyDict = sorted(np.unique(("/".join(frq.keys())).split("/")));
	mapDict = dict()
               	
	def computeDistance(frq, method, percentile):
	    return{
		'min'  : minDistance(frq),
		'prod'  : prodDistance(frq)
	    }[method]
	mapDict = computeDistance(frq,method,percentile)
	printDistanceTableToFile(mapDict,keyDict,outfile)


	
