#!/usr/bin/env python

a = open('./test/distancet.d','r')
b = open('./test/distancetr.d','r')
n = 0
for i, j in zip(a,b):
	if i==j:
		n += 1
print n
	
