#!/usr/bin/env python

import dendropy
import sys
import os
from optparse import OptionParser
import numpy as np
import itertools
import random
import subprocess
import tempfile


def expandname(filename):
    src_fpath = os.path.expanduser(os.path.expandvars(filename))
    if not os.path.exists(src_fpath):
        sys.stderr.write('Not found: "%s"' % src_fpath)
    return src_fpath


if __name__ == "__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-s", "--species", dest="sp", type="string",
                      help="read species tree from FILENAME")
    parser.add_option("-o", "--output", dest="out", type="string",
                      help="the PATH to write the generated files")
    (options, args) = parser.parse_args()
    sp = options.sp
    out = options.out

    sp = expandname(sp)

    tree = dendropy.Tree.get(path=sp, schema="newick", store_tree_weights=True)
    tree.encode_bipartitions()
    taxon_namespace = tree.taxon_namespace
    f = open(out, 'w')
    tmpList = list()
    tmp = dict()
    for edge in tree.postorder_internal_edge_iter(exclude_seed_edge=True):
        if edge.length is not None:
            to_print = edge.bipartition.leafset_as_newick_string(taxon_namespace)
            print to_print
            to_print = to_print.replace(";", "")
            to_print = to_print.replace(" ", "")
            to_print = to_print.replace("),(", "|")
            to_print = to_print.replace("(", "")
            to_print = to_print.replace(")", "")
            g = to_print.split("|")
            g[0] = ",".join(sorted(g[0].split(",")))
            g[1] = ",".join(sorted(g[1].split(",")))
            to_print = "|".join(sorted([g[1], g[0]]))
            tmpList.append(to_print)
            tmp[to_print] = edge.length
    tmpList = sorted(tmpList)
    l = [0.2, 0.5, 1, 2]
    for lt in l:
        for i in range(0, len(tmpList)):
            a = 1 - 2. / 3 * np.exp(-tmp[tmpList[i]] * lt)
            print >> f, str(lt) + "X,true,true,R1,Branch-" + str(i) + ",top-0," + str(a)
            b = 1. / 3 * np.exp(-tmp[tmpList[i]] * lt)
            print >> f, str(lt) + "X,true,true,R1,Branch-" + str(i) + ",top-1," + str(b)
            print >> f, str(lt) + "X,true,true,R1,Branch-" + str(i) + ",top-2," + str(b)
