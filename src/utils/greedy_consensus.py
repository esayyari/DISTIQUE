#!/usr/bin/env python

import dendropy
import sys
import os
import copy

if __name__ == '__main__':

    if ("--help" in sys.argv) or ("-?" in sys.argv) or ("-h" in sys.argv) or len(sys.argv) < 4:
        sys.stderr.write("usage: %s [threshold] [<out-file-path>] [<trees-file-paths>] \n"%sys.argv[0])
        sys.exit(1)
    
    threshold = float(sys.argv[1])
    resultsFile = sys.argv[2]

    trees = dendropy.TreeList()
   
    for fpath in sys.argv[3:]:
        trees.read_from_path(fpath, 'newick')

    con_tree = trees.consensus(min_freq=threshold)     
        
    for e in con_tree.postorder_node_iter():
        if hasattr(e,"support"):
            e.support=e.support*100
            e.label=str(e.support)
    to_resolve=list()
    taxa = set(con_tree.leaf_nodes())
    for e in con_tree.postorder_node_iter():
	if  len(e.child_nodes())>3 or (len(e.child_nodes())==3 and e.parent_node>0) :
		print len(e.child_nodes())
		to_resolve.append(list())
		tmp_set = set()
		for x in e.child_node_iter():
			to_resolve[len(to_resolve)-1].append(x.leaf_nodes())
			tmp_set = tmp_set | set(x.leaf_nodes())
		to_resolve[len(to_resolve)-1].append(list(taxa-tmp_set))
    for a in to_resolve:
	for b in a:
		for l in b:
			print l.taxon.__str__(),
		print
		print
	print
	print
				
    con_tree.write(path=resultsFile,schema='newick',suppress_rooting=True)    
