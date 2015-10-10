#/usr/bin/env python

import dendr pth in sys.argv[3:]:
        trees.read_from_path(fpath, 'newick')

    con_tree = trees.consensus(min_freq=threshold)

    for e in con_tree.postorder_node_iter():
        if hasattr(e,"support"):
            e.support=e.support*100
            e.label=str(e.support)
    to_resolve=dict()
    tmp = list()
    taxa = set(con_tree.leaf_nodes())
    p = 0
    for e in con_tree.postorder_node_iter():
        n = len(e.child_nodes())
        tmp_set = set()
        if n>2:
                e.label = "poly"+str(j)
                j += 1y
import sys
import os
import copy
import random
from collections import defaultdict
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
    to_resolve=dict()
    tmp = list()
    taxa = set(con_tree.leaf_nodes())
    p = 0
    for e in con_tree.postorder_node_iter():
	n = len(e.child_nodes())
	tmp_set = set()
	if n>2:
		e.label = "poly"+str(j)
		j += 1
		v = dict()
		i = 0
		for c in e.child_node_iter():
			c.label=str("child"+str(i))
			if len(tmp)>0:
				t = tmp.pop()
				v["child"+str(i)] = t
				tmp_set = tmp_set | t
				i += 1
		tmp.append(tmp_set);
		if len(taxa-tmp_set)>0:
			e.parent_node.label = str("child"+str(n))
			t = taxa-tmp_set
			v["child"+str(n)] = t
		to_resolve[e] = v
	else:
		for i in range(0,n):
			if len(tmp)>0:
				tmp_set = tmp_set | tmp.pop()
		if e.is_leaf():
			tmp_set.add(e)
		tmp.append(tmp_set)
    for key,vals in to_resolve.iteritems():
	tmp_list = list()
	for key2, vals2 in vals.iteritems():
		a = random.sample(set(vals2),1)
		print a[0].taxon.__str__(),
		tmp_list.append(a[0])
	print
				
    for key, vals in to_resolve.iteritems():
	lch = key.child_nodes()
	if key.parent_node:
		lch.append(key.parent_node)
	
	for i in lch:
		print 
		print i.label+":",
		for l in vals[i.label]:
			print l.taxon.__str__(),
		print
	print
	for key2, vals2 in vals.iteritems():
		for l in vals2:
			print l.taxon.__str__(),
		print
	print
				
    con_tree.write(path=resultsFile,schema='newick',suppress_rooting=True)    
