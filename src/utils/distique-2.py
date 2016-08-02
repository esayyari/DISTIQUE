#!/usr/bin/env python
import dendropy
import sys
import os
from optparse import OptionParser
from  prodDistance import prodDistance
from minDistance import minDistance
import numpy as np
import itertools
import random
import subprocess
import tableManipulationTools as tbs
import printTools as pr
import toolsTreeTaxa as tstt
import tempfile
import timer as tm
WS_LOC_SHELL= os.environ['WS_HOME']+'/DISTIQUE/src/shell'
WS_LOC_FM = os.environ['WS_HOME']+'/fastme-2.1.4/src'
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--file", dest="filename", type="string",
              help="read quartet table from FILENAME")
parser.add_option("-g","--gene",dest="gt",type="string",
		help="read genetrees from FILENAME")
parser.add_option("-o","--output",dest="out",type="string",
		help="the PATH to write the generated files")
parser.add_option("-t","--threshold",dest="thr",type=float,
		help="the minimum frequency that consensus will use. Default is 0.5",default=0.5)
parser.add_option("-m","--method",dest="method",type=str,
		help="The method to compute the distance of taxa. The default is prod.",default="prod")
parser.add_option("-v","--verbose",dest="verbose",
		help="Verbose",default=1)
parser.add_option("-r","--readFromFile",dest="readFromFile",
		help="Set it to 1 if you already computed the quartet table. If you want the code to compute the quartet tables set it to 0. Default is 0",default=0)
parser.add_option("-a","--averagemethod",dest="av",
		help="The average method to find the average quartet table. Options are geometric mean (gmean), and root mean square (everything else). Default is mean.", default="mean")
parser.add_option("-l",dest="met",
		help = "The method to summerize quartet results around each node, freq, or log. Default is freq", default="freq")
parser.add_option("-u",dest="sumProg",
        help = "The summerize method program to find species tree from distance matrix. The options are ninja, fastme, phydstar. Default is fastme ",default="fastme") 
parser.add_option("-z",dest="sumProgOption",
                  help = "The distance method to build the tree. If sumProg is set to fastme the options are TaxAdd_(B)alME (-s), TaxAdd_(B2)alME (-n), TaxAdd_(O)LSME (-s), TaxAdd_(O2)LSME (-n), B(I)ONJ (default), (N)J. If you path D, the fastme will run without any option (default of fastme). If the  sumProg is set to phydstar, the options are BioNJ, MVR, and NJ. The default is TaxAdd_(B)alME.",default="B")
(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
outpath = options.out
thr = options.thr
thr=options.thr
av = options.av
sumProg = options.sumProg
sumProgOption = options.sumProgOption
if sumProg == "phydstar" and sumProgOption == "":
    sumProgOption = "BioNJ"
elif sumProg == "phydstar":
    sumProgOption = options.sumProgOption

met = options.met
verbose=options.verbose
if options.readFromFile == 1:
	readFromFile = True
elif options.filename:
	readFromFile = True
else:
	readFromFile = False
method = options.method

if ( not options.gt  or not options.out):
	sys.exit("Please enter genetrees file, and output folder location")
if readFromFile:
	frq = tbs.readTable(filename)

src_fpath = os.path.expanduser(os.path.expandvars(gt))

trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')
print "time to compute consensus is: "
tm.tic()
con_tree = trees.consensus(min_freq=thr)   
tm.toc()

ftmpt=tempfile.mkstemp(suffix='.nwk', prefix="consensusTree", dir=outpath, text=None)
con_tree.write(path=ftmpt[1],schema="newick",suppress_rooting=True)

os.close(ftmpt[0])

tstt.labelNodes(con_tree)

(to_resolve,maxPolyOrder) = tstt.findPolytomies(con_tree)
taxa = list()
for e in con_tree.leaf_nodes():
	taxa.append(e.taxon.label)
n = len(con_tree.leaf_nodes())
if verbose:
    print "The summary program is: "+sumProg
    print "The option for this summary program is: "+sumProgOption
    print "Number of taxa is: " + str(n)
    print "the number of polytomies is: "+str(len(to_resolve))
    print "the maximum order of polytomies is: "+str(maxPolyOrder)
if verbose:
	print "computing the total quartet table"
if readFromFile:
	frq = tbs.readTable(filename)
else:
	tm.tic()
	frq = tbs.findQuartetTable(trees,taxa,0,outpath,verbose)

	print "time to find quartet lists: "
	tm.toc()
tm.tic()
for e in con_tree.postorder_node_iter():
	if e in to_resolve:
		val = to_resolve[e]
		(taxa_list,taxa_inv) =  tstt.getTaxaList(to_resolve[e])
		if verbose:
			print "computing the partial quartet table"	
		quartTable = tbs.findTrueAverageTable(frq,taxa_list,av,met)
		if verbose:
			print "computing distance table using the method: "+str(method)
		ftmp3=tempfile.mkstemp(suffix='.d', prefix="distancet.d", dir=outpath, text=None)
		tbs.distanceTable(quartTable,method,ftmp3[1],met)
		ftmp4=tempfile.mkstemp(suffix='.nwk',prefix="distance.d_fastme_tree.nwk",dir=outpath,text=None)
		os.close(ftmp3[0])
		ftmp4=tempfile.mkstemp(suffix='.nwk',prefix="distance.d_fastme_tree.nwk",dir=outpath,text=None)
		FNULL = open(os.devnull,'w')
		
		#subprocess.call([WS_LOC_FM+"/fastme", "-i",ftmp3[1],"-w","none","-o",ftmp4[1],"-I","/dev/null"],stdout=FNULL,stderr=subprocess.STDOUT)
		#os.close(ftmp4[0])

        	tstt.buildTreeFromDistanceMatrix(ftmp3[1],ftmp4[1],sumProg,sumProgOption)
#	        os.close(ftmp3[0])
        	os.close(ftmp4[0])
	        if verbose:
            		print "starting to resolve polytomy"	
	        res= tstt.resolvePolytomy(ftmp4[1],e,verbose)	
print "resolving polytomies takes about: "
tm.toc()
ftmp=tempfile.mkstemp(suffix='.nwk', prefix="distance.d_distique_tree.nwk", dir=outpath, text=None)
con_tree.write(path=ftmp[1],schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
if verbose:
	print "writing the resulting tree as: "+outpath+"/distance.d_distique_tree.nwk"

os.close(ftmp[0])
