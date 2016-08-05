#!/usr/bin/env python
import sys
import os
if ("--help" in sys.argv) or ("-?" in sys.argv) or len(sys.argv) < 3:
    sys.stderr.write("usage: %s [<distance-matrix-1>] [<distance-matrix-2>]\n"%sys.argv[0])
    sys.exit(1)
 
src_fpath1 = os.path.expanduser(os.path.expandvars(sys.argv[1]))
src_fpath2 = os.path.expanduser(os.path.expandvars(sys.argv[2]))

if not os.path.exists(src_fpath1):
    sys.stderr.write('Not found: "%s"' % src_fpath1)
elif not os.path.exists(src_fpath2):
    sys.stderr.write('Not found: "%s"' % src_fpath2)

a = open(src_fpath1,'r')
b = open(src_fpath2,'r')
n = 0
numline = 0
for i, j in zip(a,b):
	numline += 1
	if i==j:
		n += 1
print "number of lines is: " + str(numline)
print "number of equal distance lines is: " + str(n)
	
