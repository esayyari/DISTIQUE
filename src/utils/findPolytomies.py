#!/usr/bin/env python
import dendropy

def findPolytomies(con_tree):
	to_resolve = dict()
	tmp = list()
	taxa = set(con_tree.leaf_nodes())
	for e in con_tree.postorder_node_iter():
		tmp_set = set()
		n = len(e.child_nodes())
		if n>3 or ( n==3 and e.parent_node != None ):
			v = dict()
			for c in e.child_node_iter():
#				if len(tmp)>0:
#					t = tmp.pop()
#					v[t.label] = t
				tmp_set = tmp_set | set(c.leaf_nodes())
				v[c.label] = set(c.leaf_nodes())
		#	tmp.append(tmp_set);
			if len(taxa-tmp_set)>0:
				t = taxa-tmp_set
				v[e.parent_node.label] = t
			to_resolve[e] = v
		#else:
		#	for i in range(0,n):
		#		if len(tmp)>0:
		#			tmp_set = tmp_set | set(tmp.pop())
		#	if e.is_leaf():
		#		tmp_set.add(e)
		#	tmp.append(tmp_set)


	return to_resolve
