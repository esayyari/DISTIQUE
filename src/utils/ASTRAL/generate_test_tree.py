#!/usr/bin/python

import dendropy
from dendropy.simulate import treesim
import sys

ntaxa=int(sys.argv[1])
num_reps=int(sys.argv[2])
outpath=sys.argv[3]
t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=ntaxa)
t.write(path=outpath+"/test_tree.species_tree.trees",schema="newick",suppress_rooting=True,suppress_edge_lengths=True)

gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=t.taxon_namespace,
        num_contained=1)
gene_trees = dendropy.TreeList()
for rep in range(num_reps):
	gene_tree = treesim.contained_coalescent_tree(containing_tree=t,gene_to_containing_taxon_map=gene_to_species_map)
	gene_trees.append(gene_tree)
gene_trees.write(path=outpath+"/test_tree.gene_trees.trees",schema="newick",suppress_rooting=True,suppress_edge_lengths=True)
