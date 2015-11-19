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
parser.add_option("-n","--numStep",dest="num",type="int",
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
tm.toc()

tm.tic()
print "majority consensus"
con_tree = trees.consensus(min_freq=thr)   
tstt.labelNodes(con_tree)

con_tree.write(path=outpath+"/consensusTree.nwk",schema="newick") 

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
			fileDistance = outpath+"/distancet-"+str(anch[0])+"-"+str(anch[1])+".d"
			pr.printDistanceTableToFile(D,keyDict,fileDistance)
			subprocess.call([WS_LOC_FM+"/fastme", "-i",fileDistance,"-w","none","-o",fileDistance+"_fastme_tree.nwk"])
			if verbose:
				print "starting to resolve polytomy"	
			res=atbs.resolvePolytomy(fileDistance+"_fastme_tree.nwk",e,verbose)	
			if verbose:
				print res
	(a_nodes,seed_lab,pa1,pa2,par1_is_Poly,par2_is_Poly,par1_child,par2_child) = con_map	
	p1_pre_child = {xy for xy in par1_child} 
	p2_pre_child = {xy for xy in par2_child}
	#(p1child,p2child)=atbs.addAnchores(con_tree_tmp,con_map)
	#print p1child
	#print p2child
	#p1_post_child = {xy.label for xy in p1child}
	#p2_post_child = {xy.label for xy in p2child}
	#tstt.prune_tree_trivial_nodes(con_tree_tmp)	
	print "writing the resulting tree as: "+outpath+"/distance-"+str(anch[0])+"-"+str(anch[1])+".d_fastme_tree.nwk"
	con_tree_tmp.write(path=outpath+"/distance-"+str(anch[0])+"-"+str(anch[1])+".d_fastme_tree.nwk",schema="newick")
	 
	#res2 = tstt.compareRes(outpath+"/distance-"+str(anch[0])+"-"+str(anch[1])+".d_fastme_tree.nwk",taxa,anch,sp,outpath)
	res=tstt.compareAnchoredRes(outpath+"/distance-"+str(anch[0])+"-"+str(anch[1])+".d_fastme_tree.nwk",taxa,anch,sp,outpath)
#	if p1_post_child == p1_pre_child and p2_post_child == p2_pre_child:
#		print True
#		print "parent of anchores have the same children"
#	else:
#		print False
#		print "Parent of anchores do not have the same children"

	print res
	#print res2

