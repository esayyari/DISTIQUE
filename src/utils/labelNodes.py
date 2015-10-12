#!/usr/bin/env python
import dendropy
def labelNodes(tree):
	i = 0
	name = list()
	for e in tree.postorder_node_iter():
		if e.taxon != None:
			e.label = e.taxon.label
			name.append(e.label)
		else:
			e.label = "node"+str(i)
			i += 1
	return 
