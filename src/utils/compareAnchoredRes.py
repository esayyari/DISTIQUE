import dendropy
def compareAnchoredRes(tree,taxa,achs,sp,outpath):
	taxa = set(taxa)-{achs[0],achs[1]}
	tns = dendropy.TaxonNamespace()
	
	
#	taxa = tree.taxon_namespace
	tree1 = dendropy.Tree.get_from_path(sp,"newick",taxon_namespace=tns,rooting="force-unrooted")
	inferedTree = tree1.clone(2)
	inferedTree.retain_taxa_with_labels(taxa, update_bipartitions=True)
	inferedTree.deroot()		
	inferedTree.write(path=outpath+"/sp.nwk",schema="newick")
	tns = dendropy.TaxonNamespace()
	tree1 = dendropy.Tree.get_from_path(outpath+'/sp.nwk',"newick",taxon_namespace=tns,rooting="force-unrooted")
	tree2 = dendropy.Tree.get_from_path(tree,"newick",taxon_namespace=tns,rooting="force-unrooted")
	res = dendropy.calculate.treecompare.false_positives_and_negatives(tree1,tree2)

	
	return res
