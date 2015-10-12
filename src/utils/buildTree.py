#!/usr/bin/env python
import dendropy
def buildTree(setNodeLabels,tree,center):
	inferedTree = tree.clone(2)
	taxa = dendropy.TaxonNamespace()
	for node in inferedTree.postorder_node_iter():
		if node.label == center.label:
			center = node
			break
	inferedTree.reroot_at_node(center,update_bipartitions=False, suppress_unifurcations=False)
	listNodes = list()	
	for node in inferedTree.preorder_node_iter():
		if node.label in setNodeLabels:
			listNodes.append(node)
	for node in listNodes:
			node.clear_child_nodes()
			if node.taxon==None:
				vt = node.label
				tmp = dendropy.Taxon(label=vt)
				inferedTree.taxon_namespace.add_taxon(tmp)
				taxa.add_taxon(tmp)
				node.taxon = tmp
			else:
				tmp = node.taxon
				taxa.add_taxon(tmp)
	
	inferedTree.retain_taxa(taxa,update_bipartitions=True)
	inferedTree.deroot()		
	#for node in inferedTree:
	#	node.label=None 
	return inferedTree
