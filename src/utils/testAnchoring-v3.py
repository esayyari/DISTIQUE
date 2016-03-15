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
(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
outpath = options.out
fillmethod = options.fillmethod
thr = float(options.thr)
sp = options.sp
num = options.num
met = options.met
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
while(notEnoughSample):
    if randomSample:
        ac = tstt.random_combination(itertools.combinations(taxa,2),num)
        ac = tstt.pickAnchors(taxa,to_resolve,num,debugFlag)
    n = len(con_tree.leaf_nodes())
    

    distance_tables = dict()            
    if verbose:
        print "hello!"
        print "Number of taxa is: " + str(n)
        print "Number of polytomies is: " + str(len(to_resolve))
        for k in to_resolve:
            print "The size of polytomy around node "+k.label+" is: "+str(len(to_resolve[k].keys()))
        print "Filling method is: " + fillmethod
        print "Distance method is: "+ method
        print "Averaging distances around polytomies using method: "+ met
        print "The number of anchores are: "+str(len(ac))
    count_distance_table = dict()
    for anch in ac:
        tm.tic()

        anch = sorted(list(anch))
        anch_temp = "/".join(anch)
        print anch
        print "time elapsing  for counting number of quartets for this anchors is: "
        tm.tic()
        if not readFromFile:
            [_, frq] = atbs.findAnchoredDistanceTable(anch, trees, taxa, outpath)
        else:
            frq = atbs.findAnchoredDistanceTableFromFile(anch,frqT,taxa,outpath)
        tm.toc()

        for e in to_resolve:
            val = to_resolve[e]
            
            
            (taxa_list, taxa_inv) = tstt.getTaxaList(val)
            for nd in taxa_list:
                if anch[0] in taxa_list[nd]:
                    N1 = nd
                if anch[1] in taxa_list[nd]:
                    N2 = nd
            if N1 == N2 or len(taxa_list.keys())<6:
                if len(taxa_list.keys())<6:
                    if verbose:
                        print "The size of polytomy around: "+e.label+" was less than 6, resolving it with another method!"
                        skippedPoly.add(e)
                else:
                    if verbose:
                        print "The anchors are not in two different clades around the polytomy "+e.label
                continue
            N = {N1,N2}
            quartTable = tbsa.findTrueAverageTableAnchoringAddDistances(frq,anch,taxa_list,N,method,met)
            if debugFlag:
                outpath1 = outpath + "/quartTable_Main-"+anch[0]+"-"+anch[1]+"-"+e.label+".qt"
                tD.printQuartTable(frq,outpath1)
                outpath2 = outpath + "/quartTable_Main_Processed-"+anch[0]+"-"+anch[1]+"-"+e.label+".qt"
                tD.printQuartetTableAveraged(quartTable,outpath2)
                outpath2 = outpath + "/listTaxa-"+e.label+".lt"
                tD.printTaxaList(taxa_list,outpath2)
            [Dtmp,Ctmp] = atbs.anchoredDistanceFromFrqAddDistances(quartTable,anch,taxa_list)
            if debugFlag:
                outpath2 = outpath +"/DistanceTable-"+anch[0]+"-"+anch[1]+"-"+e.label+".dt"
                outpath3 = outpath +"/CountTable-"+anch[0]+"-"+anch[1]+"-"+e.label+".ct"
                tD.printDistanceTable(Dtmp,Ctmp,outpath2,outpath3)
            atbs.fillEmptyElementsDistanceTable(Dtmp,Ctmp,fillmethod)
            if debugFlag:
                outpath3 = outpath +"/FilledCountTable-"+anch[0]+"-"+anch[1]+"-"+e.label+".fct"
                outpath2 = outpath+"/FilledDistanceTable-"+anch[0]+"-"+anch[1]+"-"+e.label+".fdt"
                tD.printDistanceTable(Dtmp,Ctmp,outpath2,outpath3)
            if e.label in distance_tables:
                atbs.addDistanceAnchores(distance_tables[e.label],Dtmp,count_distance_table[e.label],Ctmp,fillmethod)
            else:
                if fillmethod == "normConst":
                    Dmax = max(np.abs(Dtmp.values()))*1.
                    if Dmax == 0:
                        Dmax = 1.
                    for kttDtmp in Dtmp:
                        Dtmp[kttDtmp] = Dtmp[kttDtmp]/Dmax
                    distance_tables[e.label] = Dtmp
                    count_distance_table[e.label] = Ctmp
                else:
                    distance_tables[e.label] = Dtmp
                    count_distance_table[e.label] = Ctmp
        print "Computing distance table using anchors "+anch[0]+" and "+anch[1]+" has been finished!"
        tm.toc()
    if verbose:
        print "Number of anchors is: "+str(len(ac))
    normalizedD = distance_tables
    normalizedC = count_distance_table
    for e in to_resolve:
        if e in skippedPoly:
            continue
        flag=atbs.normalizeDistanceTable(normalizedD[e.label],normalizedC[e.label])
        if flag:
             warnings.warn("One of the distances has been estimated poorly! repeating sampling anchors")
#              break
#     if flag:
#         continue
    notEnoughSample = False
fileAnch = "listAnchors"
ftmp3=tempfile.mkstemp(suffix='.txt', prefix=fileAnch, dir=outpath, text=None)
f = open(ftmp3[1], 'w')
for anch in ac:
    anch = sorted(list(anch))
    anch_temp = "|".join(anch)
    print >> f, anch_temp
f.close()
os.close(ftmp3[0])

for e in con_tree.postorder_node_iter():
        if e in to_resolve:
            if e in skippedPoly:
                continue
            keyDict = sorted(list(np.unique((" ".join(normalizedD[e.label].keys())).split(" "))))
            fileDistance = "distancet-"+str(e.label)+".d"
            ftmp3=tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=outpath, text=None)
            pr.printDistanceTableToFile(normalizedD[e.label],keyDict,ftmp3[1])
            print "writing distance table to "+str(ftmp3[1])
            os.close(ftmp3[0])
            ftmp4=tempfile.mkstemp(suffix='.nwk', prefix=fileDistance+"_fastme_tree.nwk",dir=outpath,text=None)
            FNULL = open(os.devnull, 'w')
            subprocess.call([WS_LOC_FM+"/fastme", "-i",ftmp3[1],"-w","none","-o",ftmp4[1],"-I","/dev/null"],stdout=FNULL,stderr=subprocess.STDOUT)
            os.close(ftmp4[0])
            if verbose:
                print "starting to resolve polytomy"
            res=atbs.resolvePolytomy(ftmp4[1],e,verbose)
tstt.prune_tree_trivial_nodes(con_tree)
print "writing the resulting tree as: "+outpath+"/distance.d_distique_anchoring_tree.nwk"
ftmp=tempfile.mkstemp(suffix='.nwk', prefix="distance.d_distique_anchoring_tree.nwk", dir=outpath, text=None)
con_tree.write(path=ftmp[1],schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
os.close(ftmp[0])

