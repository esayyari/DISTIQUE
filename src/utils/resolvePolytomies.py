#!/usr/bin/env python
import sys
import dendropy
import os
from buildTree import buildTree
import subprocess
def resolvePolytomy(pathToTree,node,otr,verbose):
	src_fpath = os.path.expanduser(os.path.expandvars(pathToTree))
	if not os.path.exists(src_fpath):
   		 sys.stderr.write('Not found: "%s"' % src_fpath)
	tlist = dendropy.TreeList()
	tlist = dendropy.TreeList.get(path=src_fpath,schema="newick")
	sp_tree = tlist[0]
	adjacent_list = set()
	dict_children=dict()
	tn = 0
	for t in node.adjacent_nodes():		
		dict_children[t.label] = t
		adjacent_list.add(t.label)
	if node.parent_node is not None:
		label = node.parent_node.label
		filter = lambda taxon: True if taxon.label==label else False
		nd = sp_tree.find_node_with_taxon(filter)
		if nd is not None:
			sp_tree.reroot_at_edge(nd.edge, update_bipartitions=False)
	stack = list()
	for e in sp_tree.postorder_node_iter():
		n = len(e.child_nodes())
		if len(node.child_nodes())<=2:
			continue
		if n > 0:
			tmp_next = str()
			t = node.insert_new_child(n+1)
			for i in range(0,n):
				if len(stack)>0:
					tmp = stack.pop()
					if dict_children[tmp]  in node.child_nodes():
						node.remove_child(dict_children[tmp])
					elif len(node.child_nodes())==1:
						continue
					t.add_child(dict_children[tmp])
					tmp_next += tmp
			t.label = tmp_next
			dict_children[tmp_next] = t
			stack.append(tmp_next)
		elif e.is_leaf():
			e_t = e.taxon.__str__()
			e_t = e_t.split("'")
			stack.append(e_t[1])	
	if verbose:
	#	inferedTree = buildTree(adjacent_list,otr,node)
	#	inferedTree.write(path="./tmp.nwk",schema="newick")
	#	tns = dendropy.TaxonNamespace()
	#	tree1 = dendropy.Tree.get_from_path(src_fpath,"newick",taxon_namespace=tns,rooting="force-unrooted")
	#	tree2 = dendropy.Tree.get_from_path("tmp.nwk","newick",taxon_namespace=tns,rooting="force-unrooted")
	#	res = dendropy.calculate.treecompare.false_positives_and_negatives(tree1,tree2)
		print "done"	
		return None 
	else:
		print "done"
		return None 
		

	
