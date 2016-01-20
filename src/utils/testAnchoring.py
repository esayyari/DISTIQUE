#!/usr/bin/env python
import dendropy
import sys
import os
import numpy as np
import itertools
import subprocess
import printTools as pr
import tableManipulationTools as tbs
import anchoredTableTools as atbs
import toolsTreeTaxa as tstt
import timer as tm
from optparse import OptionParser
import tableManipulationToolsAnchoring as tbsa
import tempfile
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
parser.add_option("-v","--verbose",dest="verbose",
		help="Verbose",default=1)
parser.add_option("-a","--achs",dest="a",type="string",
		help="Anchors that the table will be based on")
parser.add_option("-s","--sp",dest="sp",
		help="species tree")
parser.add_option("-n","--numStep",dest="num",
		help="The number of anchors, default is 2",default=2)
parser.add_option("-e","--method",dest="am",
		help="The averaging method for finding average quartet table",default="mean")
(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
outpath = options.out
thr = options.thr
sp = options.sp
num = options.num
if (num != "all"):
	num = int(num)
am = options.am
if (options.a):
	ac = sorted(options.a.split(','))
	ac = [(ac[i],ac[i+1])for i in range(0,len(ac)/2)]
	randomSample=False
else:
	randomSample=True
if options.filename:
	readFromFile = True
else:
	readFromFile = False
print readFromFile
verbose=options.verbose
if ( not options.gt  or not options.out):
	sys.exit("Please enter genetrees file, and output folder location")

src_fpath = os.path.expanduser(os.path.expandvars(gt))

tm.tic()
print "reading trees"
trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')
num_taxa = len(trees[0].leaf_nodes())
if (num == "all"):
	num = num_taxa*(num_taxa-1)/2
tm.toc()

tm.tic()
print "majority consensus"
con_tree = trees.consensus(min_freq=thr)   
tstt.labelNodes(con_tree)

ftmp=tempfile.mkstemp(suffix='.nwk', prefix="consensusTree", dir=outpath, text=None)
con_tree.write(path=ftmp[1],schema="newick",suppress_rooting=True) 

taxa = list()
for e in con_tree.leaf_nodes():
	taxa.append(e.taxon.label)
if randomSample:
	ac = tstt.random_combination(itertools.combinations(taxa,2),num)	
n = len(con_tree.leaf_nodes())
if verbose:
	print "Number of taxa is: " + str(n)

tm.toc()

if verbose:
	print "Number of taxa is: " + str(n)
for anch in ac:
	anch = sorted(list(anch))
	print anch
	con_tree_tmp = con_tree.clone(2)
	#(par1,par2,par1_is_Poly,par2_is_Poly,par1_child,par2_child)	
	(to_resolve,maxPolyOrder,con_map) = atbs.findPolytomies(con_tree_tmp,taxa,anch)
		
	for e in to_resolve:
		v=to_resolve[e]
	if verbose:
		print "computing the distance table, anchoring seperately"
		tm.tic()
		[D,frq]=atbs.findAnchoredDistanceTable(anch,trees,taxa,outpath)
		tm.toc()
	for e in con_tree_tmp.postorder_node_iter():
		if e in to_resolve:
			val = to_resolve[e]
			(taxa_list,taxa_inv) =  tstt.getTaxaList(val)
			if verbose:
				print "computing the partial quartet table"
			
			quartTable = tbsa.findTrueAverageTableAnchoring(frq,anch,taxa_list,am)

			if verbose:
				print "computing distance table using the method: "+str(am)
			D=atbs.anchoredDistanceFromFrq(quartTable,anch)
			keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
			fileDistance = "distancet-"+str(anch[0])+"-"+str(anch[1])+".d"
			ftmp3=tempfile.mkstemp(suffix='', prefix=fileDistance, dir=outpath, text=None)
			pr.printDistanceTableToFile(D,keyDict,ftmp3[1])
			ftmp4=tempfile.mkstemp(suffix='', prefix=fileDistance+"_fastme_tree.nwk",dir=outpath,text=None)
			subprocess.call([WS_LOC_FM+"/fastme", "-i",ftmp3[1],"-w","none","-o",ftmp4[1]])
			if verbose:
				print "starting to resolve polytomy"	
			res=atbs.resolvePolytomy(ftmp4[1],e,verbose)	
			if verbose:
				print res
	(num_add,ach_a)=atbs.addAnchores(con_tree_tmp,con_map)
	tstt.prune_tree_trivial_nodes(con_tree_tmp)	
	print "writing the resulting tree as: "+outpath+"/distance-"+str(anch[0])+"-"+str(anch[1])+".d_fastme_tree.nwk"
	ftmp=tempfile.mkstemp(suffix='.nwk', prefix="distance-"+str(anch[0])+"-"+str(anch[1])+".d_fastme_tree.nwk", dir=outpath, text=None)
	con_tree_tmp.write(path=ftmp[1],schema="newick",suppress_rooting=True)
	 
#	res2 = tstt.compareAnchoredRes(outpath+"/distance-"+str(anch[0])+"-"+str(anch[1])+".d_fastme_tree.nwk",taxa,anch,sp,outpath,anch)
#	ach_al = [a.label for a in ach_a]
#	res=tstt.compareAnchoredRes(outpath+"/distance-"+str(anch[0])+"-"+str(anch[1])+".d_fastme_tree.nwk",taxa,ach_al,sp,outpath,anch)
#	if p1_post_child == p1_pre_child and p2_post_child == p2_pre_child:
#		print True
#		print "parent of anchores have the same children"
#	else:
#		print False
#		print "Parent of anchores do not have the same children"
#	print num_add
#	print len(ach_al)
#	print len(con_tree_tmp.leaf_nodes())

#	print res
#	print res2

