#! /usr/bin/env python

import dendropy
import os
import sys
import re


if ("--help" in sys.argv) or ("-?" in sys.argv) or len(sys.argv) < 3:
    sys.stderr.write("usage: %s [<tree-file-path>] [<out-file-path>]\n"%sys.argv[0])
    sys.exit(1)
 
src_fpath = os.path.expanduser(os.path.expandvars(sys.argv[1]))
if not os.path.exists(src_fpath):
    sys.stderr.write('Not found: "%s"' % src_fpath)
src = open(src_fpath,"r")        
 
dest_fpath = os.path.expanduser(os.path.expandvars(sys.argv[2]))
dest = open(dest_fpath, "w")

print "Will write to file %s" %os.path.abspath(dest_fpath)


trees = dendropy.TreeList()
for tree_file in [src_fpath]:
    trees.read_from_path(
            tree_file,
            'newick')
con_tree = trees.consensus(min_freq=0.5, trees_splits_encoded=False)
for e in con_tree.postorder_internal_node_iter():
    e.label = int(round(float(e.label) * 6)) if e.label is not None else None
con = con_tree.as_string('newick')
dest.write(con)

dest.close()
