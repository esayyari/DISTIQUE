#!/usr/bin/env python
import dendropy

def findPolytomies(con_tree):
	to_resolve = dict()
	maxPolyOrder = 0
	tmp = list()
	node_dict = dict()
	taxa = set(con_tree.leaf_nodes())
	for e in con_tree.postorder_node_iter():
		tmp_set = set()
		n = len(e.child_nodes())
		if n>3 or ( n==3 and e.parent_node != None ):
			if e.parent_node != None:
				sz = n + 1
			else:
				sz = n
			maxPolyOrder = max(maxPolyOrder,sz)
			v = dict()
			for c in e.child_nodes():
				if len(tmp)>0:
					t = tmp.pop()
					
					v[c.label] = node_dict[c.label]
					tmp_set = tmp_set | t
			tmp.append(tmp_set);
			node_dict[e.label] = tmp_set
			if len(taxa-tmp_set)>0:
				t = taxa-tmp_set
				v[e.parent_node.label] = t
			to_resolve[e] = v
		else:
			for i in range(0,n):
				if len(tmp)>0:
					tmp_set = tmp_set | set(tmp.pop())
			if e.is_leaf():
				tmp_set.add(e)
			node_dict[e.label] = tmp_set
			tmp.append(tmp_set)


	return (to_resolve,maxPolyOrder)
