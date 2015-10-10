import dendropy
import os

def resolvePolytomy(pathToTree,node,otr):
	src_fpath = os.path.expanduser(os.path.expandvars(pathToTree))
	if not os.path.exists(src_fpath):
   		 sys.stderr.write('Not found: "%s"' % src_fpath)
	tlist = dendropy.TreeList()
	tlist = dendropy.TreeList.get(path=src_fpath,schema="newick")
	sp_tree = tlist[0]
	dict_children=dict()
	tn = 0
	for t in node.adjacent_nodes():		
		dict_children[t.label] = t
	if node.parent_node>0:
		filter = lambda taxon: True if taxon.label==node.parent_node.label else False
		nsp_tmp = sp_tree.find_node_with_taxon(filter)
		sp_tree.reroot_at_edge(nsp_tmp.edge, update_bipartitions=False)
	stack = list()
#	print nsp_tmp.taxon.label
	for e in sp_tree.postorder_node_iter():
		n = len(e.child_nodes())
		if len(node.child_nodes())<=2:
			return
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
					print tmp
			t.label = tmp_next
			dict_children[tmp_next] = t
			stack.append(tmp_next)
		elif e.is_leaf():
			e_t = e.taxon.__str__()
			e_t = e_t.split("'")
			stack.append(e_t[1])
#	node.add_child(t)	
	print "done"	
	return 
		

	
