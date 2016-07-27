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
import subprocess
import printTools as pr
import anchoredTableTools as atbs
import toolsTreeTaxa as tstt
import tableManipulationToolsAnchoring as tbsa
import testDependencies as tD
import warnings
import dendropy
import gc
from compiler.ast import Node
import prodDistance as pd



WS_LOC_SHELL= os.environ['WS_HOME']+'/DISTIQUE/src/shell'
WS_LOC_FM = os.environ['WS_HOME']+'/fastme-2.1.4/src'


usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t","--strategy",dest="strat",type="int",
		help="The version of DISTIQUE to be run 1 (all-paris, prod), 2 (all-paris, max), 3 (Distance-sum), and 4 (Tree-sum), default is DISTIQUE Distance-SUM (3)",default=3)
parser.add_option("-f", "--file", dest="filename", type="string",
	        help="read quartet table from FILENAME")
parser.add_option("-g","--gene",dest="gt",type="string",
		help="read genetrees from FILENAME")
parser.add_option("-o","--output",dest="out",type="string",
		help="the PATH to write the generated files")
parser.add_option("-y","--threshold",dest="thr",type=float,
		help="the minimum frequency that consensus will use. Default is 0.5",default=0.5)
parser.add_option("-v","--verbose",dest="verbose",
		help="Verbose",default=1)
parser.add_option("-u",dest="sumProg",
		help = "The summerize method program to find species tree from distance matrix. The options are ninja, fastme, phydstar. Default is fastme ",default="fastme") 
parser.add_option("-z",dest="sumProgOption",
		help = "The distance method to build the tree. If sumProg is set to fastme the options are TaxAdd_(B)alME (-s), TaxAdd_(B2)alME (-n), TaxAdd_(O)LSME (-s), TaxAdd_(O2)LSME (-n), B(I)ONJ (default), (N)J. The default in this case is B(I)ONJ. if the  sumProg is set to phydstar, the options are BioNJ, MVR, and NJ. The default is TaxAdd_(B)alME.",default="B")
parser.add_option("-s","--sp",dest="sp",
		help="species tree")
parser.add_option("-n","--numStep",dest="num",
		help="The number of anchors, default is 2",default=2)
parser.add_option("-d",dest="debug",
		help = "The flag for indicating that this run is for debugging!", default = False)
parser.add_option("-r",dest="outlier",
		help = "The strategy for outlier removal. The options are pairwise1, pairwise2, consensus10, or consensus3. Default is None", default = "consensus3")
parser.add_option("-x",dest="summary",
		help = "The summary method that will be used to summarize inferred species trees. Default is mrl",default = "mrl")
parser.add_option("-a","--averagemethod",dest="av",
		help="The average method to find the average quartet table. Default is mean.", default="mean")

parser.add_option("-e","--method",dest="am",
		help="The averaging method for finding average quartet table",default="mean")
parser.add_option("-p",dest="met",
		help="The method to summerize quartet results around each node, freq, or log, Default is freq", default="freq")

parser.add_option("-l",dest="fillmethod",
		help="The method to fill empty cells in distance tables, const, rand, or normConst. Default is const", default="const")
parser.add_option("-m","--distmethod",dest="method",type=str,
		help="The method to compute the distance of taxa. The default is prod.",default="prod")
		
(options,args) = parser.parse_args()
strat = options.strat

if (strat == 1):

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
	if options.filename:
		readFromFile = True
	else:
		readFromFile = False
	method = "prod"

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
			subprocess.call([WS_LOC_FM+"/fastme", "-i",ftmp3[1],"-w","none","-o",ftmp4[1],"-I","/dev/null"],stdout=FNULL,stderr=subprocess.STDOUT)
			os.close(ftmp4[0])

			if verbose:
				print "starting to resolve polytomy"	
			res= tstt.resolvePolytomy(ftmp4[1],e,verbose)	
	print "resolving polytomies takes about: "
	tm.toc()
	outfile = outpath+"/distique_all_pairs_prod.nwk"
        if verbose:
            print "writing the resulting tree as: "+outfile
        con_tree.write(path=outfile,schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
if (strat == 2):
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
	if options.filename:
		readFromFile = True
	else:
		readFromFile = False
	method = "min"

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
			subprocess.call([WS_LOC_FM+"/fastme", "-i",ftmp3[1],"-w","none","-o",ftmp4[1],"-I","/dev/null"],stdout=FNULL,stderr=subprocess.STDOUT)
			os.close(ftmp4[0])

			if verbose:
				print "starting to resolve polytomy"	
			res= tstt.resolvePolytomy(ftmp4[1],e,verbose)	
	print "resolving polytomies takes about: "
	tm.toc()
	
	outfile = outpath+"/distique_all_pairs_max.nwk"
        if verbose:
            print "writing the resulting tree as: "+outfile
        con_tree.write(path=outfile,schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
if (strat == 3):
	tm.tic()
	filename = options.filename
	gt = options.gt
	outpath = options.out
	fillmethod = options.fillmethod
	thr = float(options.thr)
	sp = options.sp
	num = options.num
	met = options.met
	sumProg = options.sumProg
	sumProgOption = options.sumProgOption
	if sumProg == "phydstar" and sumProgOption == "":
	    sumProgOption = "BioNJ"
	elif sumProg == "phydstar":
	    sumProgOption = options.sumProgOption

	debugFlag = (options.debug == "1")
	if (num != "all"):
	    num = int(num)
	method = options.am
	am = options.am
	randomSample=True
	if options.filename:
	    readFromFile = True
	    frqT=tbsa.readTable(filename)
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
	tm.toc()
	ftmp=tempfile.mkstemp(suffix='.nwk', prefix="consensusTree", dir=outpath, text=None)
	print con_tree

	taxa = list()
	for e in con_tree.leaf_nodes():
	    taxa.append(e.taxon.label)
	notEnoughSample = True
	skippedPoly = set()
	(to_resolve, _) = tstt.findPolytomies(con_tree)

	acSmall = dict()
	(to_resolvettt,_)= tstt.findPolytomies_with_names(con_tree)
	numAnchors = 0
	if randomSample:
	    ac = tstt.pickAnchors(taxa,to_resolve,num,debugFlag)
	smallAnchs =set()
	for e in to_resolvettt:
	    if len(to_resolvettt[e].keys())<6:
		val = to_resolvettt[e]
		(taxa_list,taxa_inv) =  tstt.getTaxaList(val)
		acSmall[e] = tstt.chooseAnchoresAll(taxa_list,1,debugFlag)
		for acList in acSmall[e]:
		    for anch in acList:
			smallAnchs.add(anch)

			if len(ac) == 0:
			    continue
			ac.append(anch)

	n = len(con_tree.leaf_nodes())
	ac = set(ac)
	if len(ac) == 0:
	    for e in acSmall:
		numAnchors += len(acSmall[e][0])*1
	else:        
	    numAnchors = len(ac)
	i = 0
	distance_tables = dict()            
	if verbose:
	    print "Number of taxa is: " + str(n)
	    print "Number of polytomies is: " + str(len(to_resolve))
	    for k in to_resolve:
		print "The size of polytomy around node "+k.label+" is: "+str(len(to_resolve[k].keys()))
		if len(to_resolve[k].keys())<6:
		    skippedPoly.add(k)
	    print "Filling method is: " + fillmethod
	    print "Distance method is: "+ method
	    print "Averaging distances around polytomies using method: "+ met
	    print "The number of anchors are: "+str(numAnchors)
	    print "The number of anchors for small polytomies is: "+str(len(smallAnchs))
	count_distance_table = dict()
	computedAnchors = dict()
	k=0
	anchPoly = list()
	anchsList = list(ac)
	distance_tables = []
	listPoly = to_resolve.keys()
	for z in range(len(listPoly)):
	    distance_tables.append([])
	if len(ac) !=0:
	    count = 1
	    for z in range(len(listPoly)):
		anchPoly = list()
		e = listPoly[z]
		if e in skippedPoly:
		    continue
		val = to_resolve[e]
		(taxa_list, taxa_inv) = tstt.getTaxaList(val)
		for anch in ac:
		    N1 = taxa_inv[anch[0]]
		    N2 = taxa_inv[anch[1]]
		    if N1 == N2:
			continue
		    N = {N1,N2}
		    if verbose:
			print "The size of polytomy is: "+str(len(taxa_list))
		    tm.tic()
		    anch = sorted(anch)
		    
		    if  readFromFile:
			frq = atbs.findAnchoredDistanceTableFromFile(anch,frqT,taxa,outpath)
		    else:
			frq = atbs.findAnchoredDistanceTableOverallp(e,N,anch,taxa_list,taxa_inv, trees,taxa, outpath,debugFlag)
			
		    print "time elapsing time for counting number of quartets for anch "+anch[0]+" "+anch[1]
		    tm.toc()

		    tm.tic()
		    if readFromFile:
			quartTable = tbsa.findTrueAverageTableAnchoringAddDistancesOverallFromFile(frq,anch,taxa_list,N,method,met)
		    else:
			quartTable = tbsa.findTrueAverageTableAnchoringAddDistancesOverall(frq,anch,taxa_list,N,method,met)
		    [Dtmp,Ctmp] = atbs.anchoredDistanceFromFrqAddDistances(quartTable,anch,taxa_list)
		    atbs.fillEmptyElementsDistanceTable(Dtmp,Ctmp,fillmethod)
		    del quartTable
		    del frq
		    
		    if len(distance_tables[z]) ==0:
			if fillmethod == "normConst":
			    Dmax = max(np.abs(Dtmp.values()))*1.
			    if Dmax == 0:
				Dmax = 1.
			    for kttDtmp in Dtmp:
				Dtmp[kttDtmp] = Dtmp[kttDtmp]/Dmax
			    distance_tables[z] = Dtmp
			    count_distance_table[z] = Ctmp
			else:
			    distance_tables[z] = Dtmp
			    count_distance_table[z] = Ctmp
		    else:
			atbs.addDistanceAnchores(distance_tables[z],Dtmp,count_distance_table[z],Ctmp,fillmethod)
		       
		    print "Computing distance table using anchors "+anch[0]+" and "+anch[1]+" has been finished!"
		    print "The anchor "+str(count)+" out of "+str(len(ac))+" anchors has been finished!"
		    count += 1
		    if (count%50)==0:
			gc.collect()
		    tm.toc()
	countT = 1
	skippedPolyList = list(skippedPoly)
	for e in skippedPolyList:
	    anchPoly = list()
	    z = listPoly.index(e)
	    val = to_resolve[e]
	    (taxa_list,taxa_inv) =  tstt.getTaxaList(val)
	    for achList in acSmall[e.label]:
		if verbose:
		    tm.tic()
		for anch in achList:
		    
		    anch = sorted(list(anch))
		    anch_temp = "/".join(anch)
		    if verbose:
			print anch
		    N1 = taxa_inv[anch[0]]
		    N2 = taxa_inv[anch[1]]
		    if N1 == N2:
			continue
		    N = {N1,N2}
		    if verbose:
			print "The size of polytomy is: "+str(len(taxa_list))      
		    if  readFromFile:
			frq = atbs.findAnchoredDistanceTableFromFile(anch,frqT,taxa,outpath)
		    else:
			frq = atbs.findAnchoredDistanceTableOverallp(e,N,anch,taxa_list,taxa_inv, trees,taxa, outpath,debugFlag)    
		    
		    if readFromFile:
			quartTable = tbsa.findTrueAverageTableAnchoringAddDistancesOverallFromFile(frq,anch,taxa_list,N,method,met)
		    else:
			quartTable = tbsa.findTrueAverageTableAnchoringAddDistancesOverall(frq,anch,taxa_list,N,method,met)            
		    [Dtmp,Ctmp] = atbs.anchoredDistanceFromFrqAddDistances(quartTable,anch,taxa_list)
		    atbs.fillEmptyElementsDistanceTable(Dtmp,Ctmp,fillmethod)
		    if len(distance_tables[z]) == 0 :
			if fillmethod == "normConst":
			    Dmax = max(np.abs(Dtmp.values()))*1.
			    if Dmax == 0:
				Dmax = 1.
			    for kttDtmp in Dtmp:
				Dtmp[kttDtmp] = Dtmp[kttDtmp]/Dmax
			    distance_tables[z] = Dtmp
			    count_distance_table[z] = Ctmp
			else:
			    distance_tables[z] = Dtmp
			    count_distance_table[z] = Ctmp
		    else:
			atbs.addDistanceAnchores(distance_tables[z],Dtmp,count_distance_table[z],Ctmp,fillmethod)
		    if verbose:
			print "The anchor "+str(countT)+" out of "+str(len(smallAnchs))+" anchors has been finished!"
		    countT += 1
		    del quartTable
		    del frq
		print "Computing distance table using anchors "+anch[0]+" and "+anch[1]+" has been finished!"
		
		
		if (countT%50 == 0):
		    gc.collect()
		print "time elapsing  for counting number of quartets for this anchors is: "
		tm.toc()

	if verbose:
	    print "Number of anchors is: "+str(len(ac))
	normalizedD = distance_tables
	normalizedC = count_distance_table
	for z  in range(len(listPoly)):
	    flag=atbs.normalizeDistanceTable(normalizedD[z],normalizedC[z])   
	    
	fileAnch = "listAnchors"
	ftmp3=tempfile.mkstemp(suffix='.txt', prefix=fileAnch, dir=outpath, text=None)
	f = open(ftmp3[1], 'w')
	for anch in ac:
	    anch = sorted(list(anch))
	    anch_temp = "|".join(anch)
	    print >> f, anch_temp
	for achList in acSmall:
	    for anch in achList:
		anch = sorted(list(anch))
		anch_temp = "|".join(anch)
		print >> f, anch_temp
	f.close()
	os.close(ftmp3[0])

	for e in con_tree.postorder_node_iter():
		if e in to_resolve:
		    z=listPoly.index(e)
		    keyDict = sorted(list(np.unique((" ".join(normalizedD[z].keys())).split(" "))))
		    fileDistance = "distancet-"+str(e.label)+".d"
		    ftmp3=tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=outpath, text=None)
		    pr.printDistanceTableToFile(normalizedD[z],keyDict,ftmp3[1])
		    print "writing distance table to "+str(ftmp3[1])
		    os.close(ftmp3[0])
		    ftmp4=tempfile.mkstemp(suffix='.nwk', prefix=fileDistance+"_fastme_tree.nwk",dir=outpath,text=None)
		    tstt.buildTreeFromDistanceMatrix(ftmp3[1],ftmp4[1],sumProg,sumProgOption)
		    os.close(ftmp4[0])
		    if verbose:
			print "starting to resolve polytomy"
		    res=atbs.resolvePolytomy(ftmp4[1],e,verbose)
	tstt.prune_tree_trivial_nodes(con_tree)
	outfile = outpath+"/distique_distance-sum.nwk"
        if verbose:
            print "writing the resulting tree as: "+outfile
        con_tree.write(path=outfile,schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)	
	print "The overall time to infer the species tree is: "
	tm.toc()
if (strat == 4):

	filename = options.filename
	gt = options.gt
	debugFlag = (options.debug == "1")
	outpath = options.out
	thr = float(options.thr)
	sp = options.sp
	num = options.num
	summary = options.summary
	sumProg = options.sumProg
	sumProgOption = options.sumProgOption
	if sumProg == "phydstar" and sumProgOption == "":
	    sumProgOption = "BioNJ"
	elif sumProg == "phydstar":
	    sumProgOption = options.sumProgOption

	met = options.met
	strategy = options.outlier
	print strategy
	if strategy is not None:
	    removeOutliers = True
	else:
	    removeOutliers = False
	if (num != "all"):
	    num = int(num)
	am = options.am
	randomSample=True
	if options.filename:
	    readFromFile = True
	    frqT=tbsa.readTable(filename)
	else:
	    readFromFile = False
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

	os.close(ftmp[0])
	taxa = list()
	for e in con_tree.leaf_nodes():
	    taxa.append(e.taxon.label)


	notEnoughSample = True
	skippedPoly = set()
	(to_resolve, _) = tstt.findPolytomies(con_tree)

	acSmall = dict()
	(to_resolvettt,_)= tstt.findPolytomies_with_names(con_tree)
	numAnchors = 0
	if randomSample:
	    ac = tstt.pickAnchors(taxa,to_resolve,num,debugFlag)
	for e in to_resolvettt:
	    if len(to_resolvettt[e].keys())<6:
		val = to_resolvettt[e]
		(taxa_list,taxa_inv) =  tstt.getTaxaList(val)
		acSmall[e] = tstt.chooseAnchoresAll(taxa_list,1,debugFlag)
		for acList in acSmall[e]:
		    for anch in acList:
			if len(ac) == 0:
			    continue
			ac.append(anch)

	n = len(con_tree.leaf_nodes())
	ac = set(ac)
	numSmallAnchors = 0
	if len(ac) == 0:
	    for e in acSmall:
		numAnchors += len(acSmall[e][0])*1 
	else:
	    for e in acSmall:
		numSmallAnchors += len(acSmall[e][0])
	    numAnchors = len(ac)

	distance_tables = dict()            
	if verbose:
	    print "Number of taxa is: " + str(n)
	    print "Number of polytomies is: " + str(len(to_resolve))
	    for k in to_resolve:
		print "The size of polytomy around node "+k.label+" is: "+str(len(to_resolve[k].keys()))
		if len(to_resolve[k].keys())<6:
		    skippedPoly.add(k)

	    print "Averaging distances around polytomies using method: "+ met
	    print "The number of anchores are: "+str(numAnchors)
	count_distance_table = dict()
	TreeList = dict()
	computedAnchors = dict()
	count = 1
	if ac is not None:
	    count = 1
	    for anch in ac:
		if verbose:
		    tm.tic()
		anch = sorted(list(anch))
		print anch
		for e in to_resolve:
		    if e in skippedPoly:
			continue
		    val = to_resolve[e]
		    (taxa_list,taxa_inv) =  tstt.getTaxaList(val)
		    N1 = taxa_inv[anch[0]]
		    N2 = taxa_inv[anch[1]]
		    if N1 == N2:
			continue
		    N = {N1,N2}
		    if verbose:
			print "The size of polytomy is: "+str(len(taxa_list))
		    tm.tic()
		    
		    if  readFromFile:
			frq = atbs.findAnchoredDistanceTableFromFile(anch,frqT,taxa,outpath)
		    else:
			frq = atbs.findAnchoredDistanceTableOverallp(e,N,anch,taxa_list,taxa_inv, trees,taxa, outpath,debugFlag)
		    if e.label not in TreeList:
			TreeList[e.label] = dendropy.TreeList()
		    if verbose:
			print "computing the partial quartet table"
		    if readFromFile:
			quartTable= tbsa.findTrueAverageTableAnchoringOnDifferentSidesOverallFromFile(frq,anch,taxa_list,N1,N2,am, met)
		    else:
			quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesOverall(frq,anch,taxa_list,N1,N2,am,met)
		    if verbose:
			print "computing distance table using the method: "+str(am)
		    D=atbs.anchoredDistanceFromFrq(quartTable,anch)

		    keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
		    fileDistance = "distancet-"+str(anch[0])+"-"+str(anch[1])+".d"
		    ftmp3=tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=outpath, text=None)
		    pr.printDistanceTableToFile(D,keyDict,ftmp3[1])
		    os.close(ftmp3[0])
		    ftmp4=tempfile.mkstemp(suffix='.nwk', prefix=fileDistance+"_fastme_tree.nwk",dir=outpath,text=None)
		    tstt.buildTreeFromDistanceMatrix(ftmp3[1],ftmp4[1],sumProg,sumProgOption)
		    os.close(ftmp4[0])
		    tree_tmp = dendropy.Tree.get(path=ftmp4[1],schema='newick',rooting="force-unrooted")
		    TreeList[e.label].append(tree_tmp)
		if verbose:             
		    tm.toc()
		    print "The anchor "+str(count)+" out of "+str(len(ac))+" anchors has been finished!"
		count += 1

	if verbose:
	    print "Start finding resolution for polytomies with degree smaller than 6"
	count = 1
	for e in skippedPoly:
	    i = 0
	    val = to_resolve[e]
	    (taxa_list,taxa_inv) =  tstt.getTaxaList(val)
	    for achList in acSmall[e.label]:
		quartTable = dict()
		D = dict()
		Frq = dict()
		if verbose:
		    tm.tic()
		for anch in achList:
		    anch = sorted(list(anch))
		    anch_temp = "/".join(anch)
		    if verbose:
			print anch
			print "time elapsing  for counting number of quartets for this anchors is: "
			tm.toc()
		    N1 = taxa_inv[anch[0]]
		    N2 = taxa_inv[anch[1]]
		    if N1 == N2:
			continue
		    N = {N1,N2}
		    if verbose:
			print "The size of polytomy is: "+str(len(taxa_list))
		    tm.tic()
		    
		    if  readFromFile:
			frq = atbs.findAnchoredDistanceTableFromFile(anch,frqT,taxa,outpath)
		    else:
			frq = atbs.findAnchoredDistanceTableOverallp(e,N,anch,taxa_list,taxa_inv, trees,taxa, outpath,debugFlag)
		    if e.label not in TreeList:
			TreeList[e.label] = dendropy.TreeList()
		    if verbose:
			print "computing the partial quartet table"
		    if readFromFile:
			tbsa.findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomiesOverallFromFile(frq,quartTable,anch,taxa_list,N1,N2,am,met)
		    else:
			tbsa.findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomiesOverall(frq,quartTable,anch,taxa_list,N1,N2,am,met)
		    if verbose:
			print "The anchor "+str(count)+" out of "+str(numSmallAnchors)+" anchors has been finished!"
		    count += 1
		if verbose:
		    print "computing distance table using the method: "+str(am)
		Frq=atbs.anchoredDistanceFromFrqSmallPolytomies(quartTable,am,met)
		D=pd.prodDistance(Frq,met)
		keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
		fileDistance = "distancet-anchList-"+str(i)+".d"
		ftmp3=tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=outpath, text=None)
		pr.printDistanceTableToFile(D,keyDict,ftmp3[1])
		os.close(ftmp3[0])
		ftmp4=tempfile.mkstemp(suffix='.nwk', prefix=fileDistance+"_fastme_tree.nwk",dir=outpath,text=None)
		tstt.buildTreeFromDistanceMatrix(ftmp3[1],ftmp4[1],sumProg,sumProgOption)
		os.close(ftmp4[0])
		tree_tmp = dendropy.Tree.get(path=ftmp4[1],schema='newick')
		if e.label in TreeList:
		    TreeList[e.label].append(tree_tmp)
		else:
		    TreeList[e.label] = dendropy.TreeList()
		    TreeList[e.label].append(tree_tmp)
		
		if verbose:  
		    tm.toc()
	if removeOutliers and verbose:
	    print "removeing outliers"

	for e in TreeList:
	    if removeOutliers:
		tstt.remove_outliers(TreeList[e],strategy,outpath,e,summary)
	for e in con_tree.postorder_node_iter():
	    if e in to_resolve:
		ftmp4=tstt.findMRL(TreeList[e.label],e.label,outpath,summary)
		if verbose:
		    print "starting to resolve polytomy"       
		res=atbs.resolvePolytomy(ftmp4,e,verbose) 
	tstt.prune_tree_trivial_nodes(con_tree)


	outfile = outpath+"/distique_tree-sum.nwk"
	if verbose:    
	    print "writing the resulting tree as: "+outfile
	con_tree.write(path=outfile,schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)

