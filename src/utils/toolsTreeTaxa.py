#!/usr/bin/env python
import dendropy
import sys
import os
import numpy as np
import random
import tempfile
from scipy.spatial import distance
from dendropy.calculate import treecompare
import itertools
from scipy.stats import mstats
import subprocess
WS_LOC_SHELL= os.environ['WS_HOME']+'/DISTIQUE/src/shell'
WS_LOC_FM = os.environ['WS_HOME']+'/fastme-2.1.4/src'
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
            for _ in range(0,n):
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
        inv_taxa=dict()
        taxa_list=dict()
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
        filt = lambda taxon: True if taxon.label==label else False
        nd = sp_tree.find_node_with_taxon(filt)
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
            for _ in range(0,n):
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
    random.seed()
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = (random.sample(xrange(n), r))
    return [pool[i] for i in indices]





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
            for _ in range(0,n):
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
def remove_outliers(treeList,strategy,outpath,e):
    print "the strategy is: "+strategy
    if len(treeList)<10:
        print "number of trees is "+str(len(treeList))+". This is not enough for outlier removal!"
        return treeList
    if strategy == "consensus10" or strategy == "consensus3":
        ftmp=findMRL(treeList,e,outpath)
        ref_tree = dendropy.Tree.get(path=ftmp,schema="newick")
        treeList.append(ref_tree)
        d = list()
        
        for tree in treeList:
            tree.encode_bipartitions()
            ref_tree.encode_bipartitions()
            res = treecompare.false_positives_and_negatives(ref_tree,tree)
            d.append(res[1])
        if strategy == "consensus3":
            mean = np.mean(d)
#             mean = mstats.mode(d)
#             mean = mean[0]
            print "the mean distance to consensus tree was: "+str(mean)
            st = np.std(d)
            print "the std of distances to consensus tree was: "+str(st)
            for i in range(len(d)-1,0,-1):
                k = (d[i] -mean)/st
                if d[i] > mean + 2.*st:
                    print "deleting "+str(i)+"th tree!"
                    print "d[i] to delete: "+str(d[i])
                    del treeList[i]
        else:
            sortIdx = np.argsort(d,0)
            print len(sortIdx)
            print sortIdx
            m = int(len(sortIdx)/5.)
            print "deleting "+str(m)+" of the trees"
            idx = sorted([x for x in sortIdx[len(sortIdx)-m:len(sortIdx)]],reverse=True)
            print idx
            print d
            for i in idx:
                print "deleting the tree "+str(i)+"the. The distance to consensus tree was: "+str(d[i])
                del treeList[i]
    elif strategy == "pairwise1" or strategy == "pairwise2" or strategy == "pairwise3":
        D = np.ndarray(shape=(len(treeList),len(treeList)),dtype=float)
        for i in range(0,len(treeList)):
            D[i][i] = 0.
            for j in range(i+1,len(treeList)):
                tree1 = treeList[i]
                tree2 = treeList[j]
                tree1.encode_bipartitions()
                tree2.encode_bipartitions()
                res1 = treecompare.false_positives_and_negatives(tree1,tree2)
                D[i][j] = res1[1]
                D[j][i] = res1[0]
        if strategy == "pairwise1":
            d = np.mean(D,1)
            
            C = np.cov(D)
            v =[distance.mahalanobis(D[:,i],d,C) for i in range(0,len(treeList))]
            print v
            sortIdx = np.argsort(v,0)
            m = int(len(sortIdx)*0.15)
            idx = sorted([x for x in sortIdx[len(sortIdx)-m:len(sortIdx)]],reverse=True)
            for i in idx:
                print "deleting the tree "+str(i)+"the. The distance to consensus tree was: "+str(v[i])
                del treeList[i]
        elif strategy == "pairwise3":
            d = np.mean(D,0)
            
            sortIdx = np.argsort(d,0)
            print len(sortIdx)
            print sortIdx
            m = int(len(sortIdx)/5.)
            print "deleting "+str(m)+" of the trees"
            idx = sorted([x for x in sortIdx[len(sortIdx)-m:len(sortIdx)]],reverse=True)
            print idx
            print d
            for i in idx:
                print "deleting the tree "+str(i)+"the. The distance to consensus tree was: "+str(d[i])
                del treeList[i]
        else:
            d = np.mean(D,0)
            print d
            mean = np.mean(d)
            st = np.std(d)
            idx = list()
            for k in range(len(d)-1,0,-1):
                if d[k] > mean+1.5*st:
                    print "deleting the tree "+str(k)+"the. The distance to consensus tree was: "+str(d[k])
                    del treeList[k]

    return treeList
def findPolytomiesNames(con_tree):
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
            to_resolve[e.label] = v
        else:
            for _ in range(0,n):
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
def pickAnchors(taxa,to_resolve,num,debugFlag):
    c = num
    anchs = list()
    taxa = set(taxa)
    for _ in range(0,c):
        for e,v in to_resolve.iteritems():
            S = set(v.keys())
            N = set(v.keys())
            while len(N)>0:
                if len(N)>1:
                    nodes=random.sample(N,2)
                else:
                    Ntmp = list(N)
                    S = S - N
                    ntmp = random.sample(S,1)
                    nodes=[Ntmp[0],ntmp[0]]
                if debugFlag:
                    "the nodes we picked around polytomy "+e.label+" are: "+nodes[0]+" " +nodes[1]
                a1 = random.sample(v[nodes[0]],1)
                a2 = random.sample(v[nodes[1]],1)

                if debugFlag:
                    print "the anchors for these choice of nodes are: "+a1[0].taxon.label+", and "+a2[0].taxon.label
                anchs.append((a1[0].taxon.label,a2[0].taxon.label))
                N = N - set(nodes)
                
    return anchs
def chooseAnchoresAll(list_taxa,num,debugFlag):
    ac = list()
    taxa = set(list_taxa.keys())
    for _ in range(0,num):
        actmp = list()
        N=itertools.combinations(taxa,2)
        for nodeAnchs in N:
            if len(list_taxa[nodeAnchs[0]]) == 1:
                a1 = list(list_taxa[nodeAnchs[0]])
            else:
                a1 = random.sample(list_taxa[nodeAnchs[0]],1)
            if len(list_taxa[nodeAnchs[1]]) == 1:
                a2 = list(list_taxa[nodeAnchs[1]])
            else:
                a2 = random.sample(list_taxa[nodeAnchs[1]],1)
            actmp.append((a1[0],a2[0]))
        ac.append(actmp)
    if debugFlag:
        for actmp in ac:
            print "printing another list!"
            for a in actmp:
                print "printing taxa in list"
                print a 
    return ac
def findMRL(treeList,e,outpath):
    ftmp3=tempfile.mkstemp(suffix='.nwk', prefix="distancet-allTreesAroundPoly-"+e+".nwk", dir=outpath, text=None)
    treeList.write(path=ftmp3[1],schema="newick",suppress_rooting=True,suppress_internal_node_labels=True)
    os.close(ftmp3[0])
    ftmp4=tempfile.mkstemp(suffix='.nwk', prefix="distancet-allTreesAroundPoly_MRL_tree"+e+".nwk",dir=outpath,text=None)
    FNULL = open(os.devnull, 'w')
    subprocess.call([WS_LOC_SHELL+"/MRL_AroundPoly.sh", "-i",ftmp3[1],"-o",ftmp4[1]],stdout=FNULL,stderr=subprocess.STDOUT)
    os.close(ftmp4[0])
    return ftmp4[1]