#!/usr/bin/env python
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
import testDependencies as tD
import warnings
import dendropy
import gc
# WS_HOME="/Users/Erfan/Documents/Research"

WS_LOC_SHELL= os.environ['WS_HOME']+'/DISTIQUE/src/shell'
WS_LOC_FM = os.environ['WS_HOME']+'/fastme-2.1.4/src'
# WS_LOC_SHELL = "/Users/Erfan/Documents/Research//DISTIQUE/src/shell"
# WS_LOC_FM = "/Users/Erfan/Documents/Research//DISTIQUE/fastme-2.1.4/src"

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
        help="The method to summerize quartet results around each node, freq, or log, Default is freq", default="freq")
parser.add_option("-l",dest="fillmethod",
        help="The method to fill empty cells in distance tables, const, rand, or normConst. Default is const", default="const")
parser.add_option("-d",dest="debug",
        help = "The flag for indicating that this run is for debugging!", default = False)
parser.add_option("-u",dest="sumProg",
        help = "The summerize method program to find species tree from distance matrix. The options are ninja, fastme, phydstar. Default is fastme ",default="fastme") 
parser.add_option("-z",dest="sumProgOption",
                  help = "The distance method to build the tree. If sumProg is set to fastme the options are TaxAdd_(B)alME (-s), TaxAdd_(B2)alME (-n), TaxAdd_(O)LSME (-s), TaxAdd_(O2)LSME (-n), B(I)ONJ (default), (N)J. The default in this case is B(I)ONJ. if the  sumProg is set to phydstar, the options are BioNJ, MVR, and NJ. The default is BioNJ.",default="B")
tm.tic()
(options,args) = parser.parse_args()
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
if (options.a):
    ac = sorted(options.a.split(','))
    ac = [(ac[i],ac[i+1])for i in range(0,len(ac)/2)]
    randomSample=False
else:
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
    count = 0
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
countT = 0
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
                print "time elapsing  for counting number of quartets for this anchors is: "
                tm.tic()
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

            del quartTable
            del frq
        print "Computing distance table using anchors "+anch[0]+" and "+anch[1]+" has been finished!"
        if verbose:
            print "The anchor "+str(countT)+" out of "+str(len(smallAnchs))+" anchors has been finished!"
        countT += 1
        
        if (countT%50 == 0):
            gc.collect()
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
#             if e not in skippedPoly:
#                 continue
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
print "writing the resulting tree as: "+outpath+"/distance.d_distique_anchoring_tree.nwk"
ftmp=tempfile.mkstemp(suffix='.nwk', prefix="distance.d_distique_anchoring_tree.nwk", dir=outpath, text=None)
con_tree.write(path=ftmp[1],schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
os.close(ftmp[0])
print "The overall time to infer the species tree is: "
tm.toc()
