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
        help="The method to summerize quartet results around each node, freq, or log", default="log")
parser.add_option("-d",dest="debug",
        help = "The debug flag",default = False)
(options,args) = parser.parse_args()
filename = options.filename
gt = options.gt
debugFlag = (options.debug == "1")
outpath = options.out
thr = options.thr
sp = options.sp
num = options.num
met = options.met
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

os.close(ftmp[0])
taxa = list()
for e in con_tree.leaf_nodes():
    taxa.append(e.taxon.label)
n = len(con_tree.leaf_nodes())
if verbose:
    print "Number of taxa is: " + str(n)
(to_resolve_t,maxPolyOrder)=tstt.findPolytomies(con_tree)
[mapTaxaToPolyNodes,mapPolyNodesToTaxa]= atbs.mapTaxaAroundPoly(to_resolve_t,debugFlag)
if verbose:
    print "Maximum Polytomy size is: "+str(maxPolyOrder)
    

if randomSample:
    ac = tstt.random_combination(itertools.combinations(taxa,2),num)
    ac = tstt.pickAnchors(taxa,to_resolve_t,num,debugFlag)    

exit 
tm.toc()
(to_resolve,_)= tstt.findPolytomies_with_names(con_tree)
if verbose:
    print "Number of taxa is: " + str(n)
TreeList = dict()
for anch in ac:
    tm.tic()
    anch = sorted(list(anch))
    print anch
    con_tree_tmp = con_tree.clone(2)
        
    if verbose:
        print "computing the distance table, anchoring seperately"
        tm.tic()
        [D,frq]=atbs.findAnchoredDistanceTable(anch,trees,taxa,outpath)
        tm.toc()
    skippedPoly = set()
    for e in to_resolve:
        val = to_resolve[e]
        (taxa_list,taxa_inv) =  tstt.getTaxaList(val)
        if e not in TreeList:
            TreeList[e] = dendropy.TreeList()
        if verbose:
            print "computing the partial quartet table"
        
        quartTable = tbsa.findTrueAverageTableAnchoring(frq,anch,taxa_list,am,met)

        if verbose:
            print "computing distance table using the method: "+str(am)
        for node in taxa_list:
            print node
            print taxa_list[node]
            if anch[0] in taxa_list[node]:
                N1 = node
            if anch[1] in taxa_list[node]:
                N2 = node
        if N1 == N2 or len(taxa_list.keys())-2<4:
            if len(taxa_list.keys())-2<4:
                skippedPoly.add(e)
            continue
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
        tree_tmp = dendropy.Tree.get(path=ftmp4[1],schema='newick')
        TreeList[e].append(tree_tmp)              
    tm.toc()

for e in con_tree.postorder_node_iter():
    if e in to_resolve_t:
        if e.label in skippedPoly:
            continue
        ftmp3=tempfile.mkstemp(suffix='.nwk', prefix="distancet-allTreesAroundPoly-"+e.label+".nwk", dir=outpath, text=None)
        TreeList[e.label].write(path=ftmp3[1],schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
        print TreeList[e.label]
        os.close(ftmp3[0])
        ftmp4=tempfile.mkstemp(suffix='.nwk', prefix="distancet-allTreesAroundPoly_MRL_tree"+e.label+".nwk",dir=outpath,text=None)
        FNULL = open(os.devnull, 'w')
        subprocess.call([WS_LOC_SHELL+"/MRL_AroundPoly.sh", "-i",ftmp3[1],"-o",ftmp4[1]],stdout=FNULL,stderr=subprocess.STDOUT)
        os.close(ftmp4[0])
        if verbose:
            print "starting to resolve polytomy"    
            
        res=atbs.resolvePolytomy(ftmp4[1],e,verbose)  
tstt.prune_tree_trivial_nodes(con_tree)

ftmp=tempfile.mkstemp(suffix='.nwk', prefix="distique_anchoring_tree.nwk", dir=outpath, text=None)
if verbose:    
    print "writing the resulting tree as: "+ftmp[1]
con_tree.write(path=ftmp[1],schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
os.close(ftmp[0])

