#!/usr/bin/env python

import dendropy
import itertools
def addQuartetTables(quartTable1,quartTable2):
	for key, val in quartTable1.iteritems():
		for key2, val2 in val.iteritems():
			quartTable1[key][key2] = float(quartTable1[key][key2]) + float(quartTable2[key][key2])
	for key, val in quartTable1.iteritems():
		a = sorted(key.split("/"))
#		print a
		sm = sum(quartTable1[key].values())
#		print key, quartTable1[key][a[1]]/sm, quartTable1[key][a[2]]/sm, quartTable1[key][a[3]]/sm
	return quartTable1
