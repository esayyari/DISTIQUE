#!/usr/bin/env python
import dendropy
import sys
import os
import numpy as np
import itertools
import subprocess
import printTools as pr
import anchoredTableTools as ATBS
import toolsTreeTaxa as tstt
import timer as tm
from optparse import OptionParser
import tableManipulationToolsAnchoring as tbsa
import tempfile
WS_LOC_SHELL = os.environ['WS_HOME']+'/DISTIQUE/src/shell'
WS_LOC_FM = os.environ['WS_HOME']+'/fastme-2.1.4/src'

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--file", dest="filename", type="string",
                  help="read quartet table from FILENAME")
parser.add_option("-g", "--gene", dest="gt", type="string",
                  help="read genetrees from FILENAME")
parser.add_option("-o", "--output", dest="out", type="string",
                  help="the PATH to write the generated files")
parser.add_option("-t", "--threshold", dest="thr", type=float,
                  help="the minimum frequency that consensus will use. Default is 0.5", default=0.5)
parser.add_option("-v", "--verbose", dest="verbose",
                  help="Verbose", default=1)
parser.add_option("-a", "--achs", dest="a", type="string",
                  help="Anchors that the table will be based on")
parser.add_option("-s", "--sp", dest="sp",
                  help="species tree")
parser.add_option("-n", "--numStep", dest="num",
                  help="The number of anchors, default is 2", default=2)
parser.add_option("-e", "--method", dest="am",
                  help="The averaging method for finding average quartet table", default="mean")
parser.add_option("-m", dest="met",
                  help="The method to summerize quartet results around each node, freq, or log", default="log")

(options, args) = parser.parse_args()
filename = options.filename
gt = options.gt
outpath = options.out
thr = options.thr
sp = options.sp
num = options.num
met = options.met
if num != "all":
    num = int(num)
am = options.am
if options.a:
    ac = sorted(options.a.split(','))
    ac = [(ac[i], ac[i+1])for i in range(0, len(ac)/2)]
    randomSample = False
else:
    randomSample = True
if options.filename:
    readFromFile = True
else:
    readFromFile = False
print readFromFile
verbose = options.verbose
if not options.gt or not options.out:
    sys.exit("Please enter genetrees file, and output folder location")

src_fpath = os.path.expanduser(os.path.expandvars(gt))

tm.tic()
print "reading trees"
trees = dendropy.TreeList.get_from_path(src_fpath, 'newick')
num_taxa = len(trees[0].leaf_nodes())
if num == "all":
    num = num_taxa*(num_taxa-1)/2
tm.toc()

tm.tic()
print "majority consensus"
con_tree = trees.consensus(min_freq=thr)
tstt.labelNodes(con_tree)

ftmp = tempfile.mkstemp(suffix='.nwk', prefix="consensusTree", dir=outpath, text=None)
con_tree.write(path=ftmp[1], schema="newick", suppress_rooting=True)

os.close(ftmp[0])
taxa = list()
for e in con_tree.leaf_nodes():
    taxa.append(e.taxon.label)
if randomSample:
    ac = tstt.random_combination(itertools.combinations(taxa, 2), num)
n = len(con_tree.leaf_nodes())
if verbose:
    print "Number of taxa is: " + str(n)

tm.toc()

if verbose:
    print "Number of taxa is: " + str(n)
tree_to_resolve = dict()
to_resolve = dict()
con_map = dict()
to_resolve_name = dict()
map_toresolve_to_toresolvename = dict()
for anch in ac:
    key = "/".join(sorted(list(anch)))
    tree_to_resolve[key] = con_tree.clone(2)
    (to_resolve[key], maxPolyOrder, con_map[key]) = ATBS.findPolytomies(tree_to_resolve[key], taxa, anch)
    to_resolve_name[key] = {[x.label] for x in  to_resolve[key].keys()}
    v = dict()
    for x in to_resolve[key].keys():
        v[x.label] = x
    map_toresolve_to_toresolvename[key] = v

(to_resolve_total, maxPolyOrder) = tstt.findPolytomies_with_names(con_tree)

for e_name in to_resolve_total:

    for anch in ac:
        anch = sorted(list(anch))
        anch_temp = "/".join(anch)
        if e_name not in to_resolve_name[anch_temp]:
            print "this is not a polytomy with respect to anchors: " + str(anch)
            continue
        tm.tic()
        print anch

        v = to_resolve[anch_temp]
        if verbose:
            print "computing the distance table, anchoring seperately"
        tm.tic()
        [D, frq] = ATBS.findAnchoredDistanceTable(anch, trees, taxa, outpath)
        tm.toc()
        
        for e in con_tree_tmp.postorder_node_iter():
            if e in to_resolve:
            val = to_resolve[e]
            (taxa_list, taxa_inv) = tstt.getTaxaList(val)
            if verbose:
                print "computing the partial quartet table"

            quartTable = tbsa.findTrueAverageTableAnchoring(frq, anch, taxa_list, am, met)

            if verbose:
                print "computing distance table using the method: "+str(am)
            D = ATBS.anchoredDistanceFromFrq(quartTable, anch)
            keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
            fileDistance = "distancet-"+str(anch[0])+"-"+str(anch[1])+".d"
            ftmp3 = tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=outpath, text=None)
            pr.printDistanceTableToFile(D, keyDict, ftmp3[1])
            os.close(ftmp3[0])
            ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix=fileDistance+"_fastme_tree.nwk", dir=outpath, text=None)
            FNULL = open(os.devnull, 'w')
            subprocess.call([WS_LOC_FM+"/fastme", "-i", ftmp3[1], "-w", "none", "-o", ftmp4[1], "-I", "/dev/null"], stdout=FNULL, stderr=subprocess.STDOUT)
            os.close(ftmp4[0])
            if verbose:
                print "starting to resolve polytomy"
            print ftmp4[1]
            res = ATBS.resolvePolytomy(ftmp4[1], e, verbose)
    tm.toc()
    (num_add, ach_a) = ATBS.addAnchores(con_tree_tmp, con_map)
    tstt.prune_tree_trivial_nodes(con_tree_tmp)
    print "writing the resulting tree as: " + outpath + "/distance-" + str(anch[0]) + "-" + str(anch[1]) + ".d_distique_anchoring_tree.nwk"
    ftmp = tempfile.mkstemp(suffix='.nwk', prefix="distance-" + str(anch[0]) + "-"+str(anch[1]) + ".d_distique_anchoring_tree.nwk", dir=outpath, text=None)
    con_tree_tmp.write(path=ftmp[1], schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)
    os.close(ftmp[0])

