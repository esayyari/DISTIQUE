#!/usr/bin/env python

import itertools

def expandDict(frq):
	expandedDict = dict()
	s = {0,1,2,3}
	for key in frq:
		l = sorted(key.split("/"))
		t = dict()
		for i in range(1,4):
			t[l[0]+l[i]] = frq[key][l[i]]
			t[l[i]+l[0]] = frq[key][l[i]]
			v = list(s - {0,i})
			t[l[v[0]]+l[v[1]]] = frq[key][l[i]]			
			t[l[v[1]]+l[v[0]]] = frq[key][l[i]]
		expandedDict[key] = t
	return expandedDict
