#!/usr/bin/env python
import dendropy
import sys
import os
import numpy as np
import itertools
from collections import defaultdict
import random
import subprocess
import tempfile
from dendropy.calculate import treecompare
def buildTree(setNodeLabels,tree,center):
    inferedTree = tree.clone(2)
    taxa = dendropy.TaxonNamespace()
    for node in inferedTree.postorder_node_iter():
        if node.taxon.label == center.taxon.label:
            center = node
            break
    inferedTree.reroot_at_node(center,update_bipartitions=False, suppress_unifurcations=False)
    listNodes = list()
    for node in inferedTree.preorder_node_iter():
        if node.taxon.label in setNodeLabels:
            listNodes.append(node)
    for node in listNodes:
            node.clear_child_nodes()
            if node.taxon==None:
                vt = node.taxon.label
                tmp = dendropy.Taxon(label=vt)
                inferedTree.taxon_namespace.add_taxon(tmp)
                taxa.add_taxon(tmp)
                node.taxon = tmp
            else:
                tmp = node.taxon
                taxa.add_taxon(tmp)

    inferedTree.retain_taxa(taxa,update_bipartitions=True)
    inferedTree.deroot()
    return inferedTree
def compareRes(tree,taxa,anch,sp,outpath):
    tns = dendropy.TaxonNamespace()
    tns = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get_from_path(sp,"newick",taxon_namespace=tns,rooting="force-unrooted")
    tree2 = dendropy.Tree.get_from_path(tree,"newick",taxon_namespace=tns,rooting="force-unrooted")
    res = treecompare.false_positives_and_negatives(tree1,tree2)


    return res

def compareAnchoredRes(tree,taxa,achs,sp,outpath,trueAnch):
    taxa = set(taxa)-set(achs)
    tns = dendropy.TaxonNamespace()
    anch = trueAnch
    tree1 = dendropy.Tree.get_from_path(sp,"newick",taxon_namespace=tns,rooting="force-unrooted")
    inferedTree = tree1.clone(2)
    inferedTree.retain_taxa_with_labels(taxa, update_bipartitions=True)
    inferedTree.deroot()
    ftmp1=tempfile.mkstemp(suffix='.nwk', prefix="sp.nwk-"+str(anch[0])+"-"+str(anch[1]), dir=outpath, text=None)
    inferedTree.write(path=ftmp1[1],schema="newick",suppress_rooting=True)
    tree2 = dendropy.Tree.get_from_path(tree,"newick",taxon_namespace=tns,rooting="force-unrooted")
    inferedTree = tree2.clone(2)
    inferedTree.retain_taxa_with_labels(taxa, update_bipartitions=True)
    inferedTree.deroot()
    ftmp2=tempfile.mkstemp(suffix='.nwk', prefix=tree+".retained", dir=outpath, text=None)
    inferedTree.write(path=ftmp2[1],schema="newick",suppress_rooting=True)
    tns = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get_from_path(ftmp1[1], taxon_namespace=tns,rooting="force-unrooted")
    tree2 = dendropy.Tree.get_from_path(ftmp2[1],"newick",taxon_namespace=tns,rooting="force-unrooted")
    res = treecompare.false_positives_and_negatives(tree1,tree2)
    return res

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
                    if c.taxon is not None:
                        v[c.taxon.label] = node_dict[c.taxon.label]
                    else:
                        v[c.label] = node_dict[c.label]
                    tmp_set = tmp_set | t
            tmp.append(tmp_set);
            if e.taxon is not None:
                node_dict[e.taxon.label] = tmp_set
            else:
                node_dict[e.label] = tmp_set
            if len(taxa-tmp_set)>0:
                t = taxa-tmp_set
                if e.taxon is not None:
                    v[e.parent_node.taxon.label] = t
                else:
                    v[e.parent_node.label] = t
            to_resolve[e] = v
        else:
            for i in range(0,n):
                if len(tmp)>0:
                    tmp_set = tmp_set | set(tmp.pop())
            if e.is_leaf():
                tmp_set.add(e)
            if e.taxon is not None:
                node_dict[e.taxon.label] = tmp_set
            else:
                node_dict[e.label] = tmp_set
            tmp.append(tmp_set)


    return (to_resolve,maxPolyOrder)
def getTaxaList(taxaDict):
        taxa = dict()
        inv_taxa=dict()
        taxa_list=dict()
        n = len(taxaDict)
        for key1,vals1 in taxaDict.iteritems():
                v = list()
                for t in vals1:
                        if t.taxon is not None:
                            k = t.taxon.label
                        else:
                            k = t.label
                        v.append(k)
                taxa_list[key1] = v
    	for k, v in taxa_list.iteritems():
                for v2 in v:
                        inv_taxa[v2] = k
		return (taxa_list,inv_taxa)
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

def resolvePolytomy(pathToTree,node,otr,verbose):
    src_fpath = os.path.expanduser(os.path.expandvars(pathToTree))
    if not os.path.exists(src_fpath):
            sys.stderr.write('Not found: "%s"' % src_fpath)
    tlist = dendropy.TreeList()
    tlist = dendropy.TreeList.get(path=src_fpath,schema="newick")
    sp_tree = tlist[0]
    adjacent_list = set()
    dict_children=dict()
    for t in node.adjacent_nodes():
	if t.taxon is not None:
	        dict_children[t.taxon.label] = t
		adjacent_list.add(t.taxon.label)
	else:
		dict_children[t.label] = t
		adjacent_list.add(t.label)
    if node.parent_node is not None:
	if node.parent_node.taxon is not None:
	        label = node.parent_node.taxon.label
	else:
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
            children = set(node.child_nodes())
            for i in range(0,n):
                if len(stack)>0:
                    tmp = stack.pop()
                    if dict_children[tmp]  in children:
                        node.remove_child(dict_children[tmp])
                    elif len(node.child_nodes())==1:
                        continue
                    t.add_child(dict_children[tmp])
                    tmp_next += tmp
	    if t.taxon is not None:
	            t.taxon.label = tmp_next
	    else:
		    t.label = tmp_next
            dict_children[tmp_next] = t
            stack.append(tmp_next)
        elif e.is_leaf():
            e_t = e.taxon.__str__()
            e_t = e_t.split("'")
            stack.append(e_t[1])
    if verbose:

        return None
    else:
        return None


def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(xrange(n), r))
    return tuple(pool[i] for i in indices)

def prune_tree_trivial_nodes(tree):
    tree.update_bipartitions()
#	for e in tree.postorder_node_iter():

#		ch = e.child_nodes()
#		p = e.parent_node
#		if len(ch)==1 and (p is not None):
#			p.add_child(ch[0])
#			e.remove_child(ch[0])
#			p.remove_child(e)
    return

def findPolytomies_with_names(con_tree):
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
            maxPolyOrder = max(maxPolyOrder, sz)
            v = dict()
            for c in e.child_nodes():
                if len(tmp)>0:
                    t = tmp.pop()
		    if c.taxon is not None:
	                    v[c.taxon.label] = node_dict[c.taxon.label]
		    else:
			    v[c.label] = node_dict[c.label]
                    tmp_set = tmp_set | t
            tmp.append(tmp_set);
	    if e.taxon is not None:
	            node_dict[e.taxon.label] = tmp_set
	    else:
		    node_dict[e.label] = tmp_set
            if len(taxa-tmp_set)>0:
                t = taxa-tmp_set
		if e.parent_node.taxon is not None:
	                v[e.parent_node.taxon.label] = t
		else:
			v[e.parent_node.label] = t
	    
            to_resolve[e.label] = v
        else:
            for i in range(0,n):
                if len(tmp)>0:
                    tmp_set = tmp_set | set(tmp.pop())
            if e.is_leaf():
                tmp_set.add(e)
	    if e.taxon is not None:
	            node_dict[e.taxon.label] = tmp_set
	    else:
       		    node_dict[e.label] = tmp_set
            tmp.append(tmp_set)


    return to_resolve, maxPolyOrder
