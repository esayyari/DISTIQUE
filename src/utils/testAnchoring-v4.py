#!/usr/bin/env python
import dendropy
import sys
import os
import numpy as np
import itertools
import subprocess
import printTools as pr
import anchoredTableTools as atbs
import toolsTreeTaxa as tstt
import timer as tm
from optparse import OptionParser
import tableManipulationToolsAnchoring as tbsa
import tempfile
from compiler.ast import Node
import prodDistance as pd
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
parser.add_option("-m",dest="met",
        help="The method to summerize quartet results around each node, freq, or log. Default is freq", default="freq")
parser.add_option("-d",dest="debug",
        help = "The debug flag",default = False)
parser.add_option("-r",dest="outlier",
        help = "The strategy for outlier removal. The options are pairwise1, pairwise2, consensus10, or consensus3. Default is None", default = "consensus3")
(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
debugFlag = (options.debug == "1")
outpath = options.out
thr = float(options.thr)
sp = options.sp
num = options.num
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
if (options.a):
    ac = options.a.split(',')
    ac = [(ac[i],ac[i+1])for i in range(0,len(ac)/2)]
    randomSample=False
else:
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
if len(ac) == 0:
    for e in acSmall:
        numAnchors += len(acSmall[e][0])*1 
else:        
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
count = 0
if ac is not None:
    count = 0
    for anch in ac:
        if verbose:
            tm.tic()
        anch = sorted(list(anch))
        anch_temp = "/".join(anch)
        print anch
        print "time elapsing time for counting number of quartets for this anchors is: "
        tm.tic()
        if not readFromFile:
            [_, frq] = atbs.findAnchoredDistanceTable(anch, trees, taxa, outpath,debugFlag)
            fileFrq = "frequency-"+anch[0]+"-"+anch[1]
            ftmp1 = tempfile.mkstemp(suffix='.frq', prefix=fileFrq, dir=outpath, text=None)
            atbs.writeFrqAnchoredOnFile(frq,ftmp1[1])
            computedAnchors[anch_temp] = ftmp1[1]
            os.close(ftmp1[0])
        else:
            frq = atbs.findAnchoredDistanceTableFromFile(anch,frqT,taxa,outpath)
        tm.toc()
        for e in to_resolve:
            if e in skippedPoly:
                continue
            val = to_resolve[e]
            (taxa_list,taxa_inv) =  tstt.getTaxaList(val)
            if e.label not in TreeList:
                TreeList[e.label] = dendropy.TreeList()
            if verbose:
                print "computing the partial quartet table"
            for nd in taxa_list:
                if anch[0] in taxa_list[nd]:
                    N1 = nd
                if anch[1] in taxa_list[nd]:
                    N2 = nd
            if N1 == N2:
                if verbose:
                        print "The anchors are not in two different clades around the polytomy "+e.label
                continue
            quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSides(frq,anch,taxa_list,N1,N2,am,met)
            if verbose:
                print "computing distance table using the method: "+str(am)
            D=atbs.anchoredDistanceFromFrq(quartTable,anch)
            keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
            fileDistance = "distancet-"+str(anch[0])+"-"+str(anch[1])+".d"
            ftmp3=tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=outpath, text=None)
            pr.printDistanceTableToFile(D,keyDict,ftmp3[1])
            os.close(ftmp3[0])
            ftmp4=tempfile.mkstemp(suffix='.nwk', prefix=fileDistance+"_fastme_tree.nwk",dir=outpath,text=None)
            FNULL = open(os.devnull, 'w')
            subprocess.call([WS_LOC_FM+"/fastme", "-i",ftmp3[1],"-w","none","-o",ftmp4[1],"-I","/dev/null"],stdout=FNULL,stderr=subprocess.STDOUT)
            os.close(ftmp4[0])
            tree_tmp = dendropy.Tree.get(path=ftmp4[1],schema='newick',rooting="force-unrooted")
            TreeList[e.label].append(tree_tmp)
        if verbose:             
            tm.toc()
            print "The anchor "+str(count)+" out of "+str(len(ac))+" anchors has been finished!"
        count += 1

if verbose:
    print "Start finding resolution for polytomies with degree smaller than 6"
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
                tm.tic()
            if not readFromFile:
                if anch_temp in computedAnchors:
                    print "using precomputed quartet table"
                    fname = computedAnchors[anch_temp]
                    frq = atbs.readFrqAnchoredOnFile(fname)
                else:
                    [_, frq] = atbs.findAnchoredDistanceTable(anch, trees, taxa,outpath,debugFlag)
            else:
                frq = atbs.findAnchoredDistanceTableFromFile(anch,frqT,taxa,outpath)
            tbsa.findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomies(frq,quartTable,anch,taxa_list,am,met)
        if verbose:
            print "computing distance table using the method: "+str(am)
        Frq=atbs.anchoredDistanceFromFrqSmallPolytomies(quartTable,am,met)
        print Frq
        D=pd.prodDistance(Frq,met)
        keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
        fileDistance = "distancet-anchList-"+str(i)+".d"
        ftmp3=tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=outpath, text=None)
        pr.printDistanceTableToFile(D,keyDict,ftmp3[1])
        os.close(ftmp3[0])
        ftmp4=tempfile.mkstemp(suffix='.nwk', prefix=fileDistance+"_fastme_tree.nwk",dir=outpath,text=None)
        FNULL = open(os.devnull, 'w')
        subprocess.call([WS_LOC_FM+"/fastme", "-i",ftmp3[1],"-w","none","-o",ftmp4[1],"-I","/dev/null"],stdout=FNULL,stderr=subprocess.STDOUT)
        os.close(ftmp4[0])
        tree_tmp = dendropy.Tree.get(path=ftmp4[1],schema='newick')
        if e.label in TreeList:
            TreeList[e.label].append(tree_tmp)
        else:
            TreeList[e.label] = dendropy.TreeList()
            TreeList[e.label].append(tree_tmp)
        if ac is None:
            if verbose:
                print "The anchor "+str(count)+" out of "+str(len(ac))+" anchors has been finished!"
            count += 1
        if verbose:
            
            tm.toc()
if removeOutliers and verbose:
    print "removeing outliers"

for e in TreeList:
    if removeOutliers:
        tstt.remove_outliers(TreeList[e],strategy,outpath,e)
for e in con_tree.postorder_node_iter():
    if e in to_resolve:
        ftmp4=tstt.findMRL(TreeList[e.label],e.label,outpath)
        if verbose:
            print "starting to resolve polytomy"    
            
        res=atbs.resolvePolytomy(ftmp4,e,verbose)  
tstt.prune_tree_trivial_nodes(con_tree)

ftmp=tempfile.mkstemp(suffix='.nwk', prefix="distique_anchoring_tree.nwk", dir=outpath, text=None)
if verbose:    
    print "writing the resulting tree as: "+ftmp[1]
con_tree.write(path=ftmp[1],schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
os.close(ftmp[0])

