import dendropy
import itertools
import sys
import os
import tableManipulationTools as tbs
import printTools as pr
import toolsTreeTaxa as tstt
import numpy as np
def anchoredDistance(**kwargs):
    readFromTable=False
    for k,v in kwargs.iteritems():
        if k == "qfile":
            qfile=v
            readFromTable = True
        elif k == "achs":
            achs = sorted(v)
        elif k == "gt":
            trees = v
        elif k == "outfile":
            outfile = v
        elif k == "wrkPath":
            out = v
        elif k == "taxa":
            taxa = v
    if(readFromTable):
        frq=tbs.readQuartetTable(qfile)
        D = anchoredDistanceFromFrq(frq,achs)
    else:
        [D,frq] = findAnchoredDistanceTable(achs,trees,taxa,out)
    keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))));
    pr.printDistanceTableToFile(D,keyDict,outfile)
    return frq
def anchoredDistanceFromFrq(frq,achs):
    D = dict()
    for k,v in frq.iteritems():
        kt = sorted(k.split("/"))
        if ((achs[0] in kt ) and( achs[1] in kt)):
            s = sorted(list(set(kt)-{achs[0],achs[1]}))
            if len(s) == 1:
                c0 = kt.count(achs[0])
                c1 = kt.count(achs[1])
                if c0 == 2:
                    s[1] = achs[0]
                elif c1 == 2:
                    s[1] = achs[1]
            elif len(s) == 0:
                s[0] = achs[0]
                s[1] = achs[1]
            key1 = s[0]+" "+s[1]
            key2 = s[1]+" "+s[0]
            D[key1] = -np.log(frq[k])
            D[key2] = D[key1]
    return D
def buildEmptyQuartets(anch,taxa,n):
    taxa = set(taxa)-set(anch)
    Q = dict()
    for pair in chooseTaxa(taxa):
        p1 = sorted(list(pair))
        p2 = sorted(anch)
        p = "/".join(sorted(p1+p2))
        Q[p] = [0.5,n+1.5]
    return Q

def chooseTaxa(taxa):
    return set(itertools.combinations(taxa, 2))
def genKey(p1,p2):
    if p1[0]<p2[0]:
            p = p1[0]+' '+p1[1]+' | '+p2[0]+' '+p2[1]
    else:
        p = p2[0]+' '+p2[1]+' | '+p1[0]+' '+p1[1]
    return p
def findAnchoredQuartets(anch,trees,taxa,out):
    anch = sorted(anch)
    n = len(trees)
    Q = buildEmptyQuartets(anch,taxa,n)
    Q2 = list()
    for tree in trees:
        rerooted=reroot(tree,anch)
        node = rerooted[0]
        root = rerooted[1]
        while(node.parent_node.label!=root.label):
            node_pre = node
            node = node.parent_node
            chs = node.child_nodes()
            chs_n = len(chs)
            if len(chs)>2:
                for i in range(0,chs_n):
                    ch = chs[i]
                    if (ch == node_pre):
                        continue
                    else:
                        listTaxa = list(ch.leaf_nodes())
                        Q= addQuartets(ch, listTaxa,Q,Q2,anch)
                    for j in range(i+1,chs_n):
                        if (chs[i] == chs[j]) or (chs[j]==node_pre):
                            continue
                        else:
                            listTaxatmp = [listTaxa,list(chs[j].leaf_nodes())]
                            Q=removeFromQuartetLentreesh(Q,listTaxatmp,anch)
            else:
                for ch in chs:
                    if (ch==node_pre):
                        continue
                    else:
                        listTaxa = list(ch.leaf_nodes())
                        Q= addQuartets(ch, listTaxa,Q,Q2,anch)
    return Q
def reroot(tree,anch):
    tstt.labelNodes(tree)
    filter = lambda taxon: True if taxon.label==anch[0] else False
    root = tree.find_node_with_taxon(filter)
    filter = lambda taxon: True if taxon.label==anch[1] else False
    node = tree.find_node_with_taxon(filter)
    tree.reroot_at_node(root, update_bipartitions=True, suppress_unifurcations=False)
    return [node,root]
def findAllChildrenPairs(listTaxa):
    listTaxaLabels = list()
    for t in listTaxa:
        listTaxaLabels.append(t.label)
    pairs = chooseTaxa(listTaxaLabels)
    return pairs
def addQuartets( ch, listTaxa,Q,Q2,anch):
    pairs = findAllChildrenPairs(listTaxa)
    for p in pairs:
        key=genKey(list(p),anch)
        t = sorted(list(p)+anch)
        Q2.append(key)
        Q["/".join(t)][0] += 1
    return Q
def removeFromQuartetLentreesh(Q,listTaxa,anch):
    x = list()
    for et in listTaxa:
        T = list()
        for t in et:
            T.append(t.label)
        x.append(T)
    keySet=[p for p in itertools.product(*x)]
    for p in keySet:
        key = sorted(list(p) + anch)
        Q["/".join(key)][1] -= 1
    return Q
def findAnchoredDistanceTable(achs,trees,taxa,out):
    frq=findAnchoredQuartets(achs,trees,taxa,out)
    D = dict()
    for k,v in frq.iteritems():
        kt = sorted(k.split("/"))
        if ((achs[0] in kt ) and( achs[1] in kt)):
            s = sorted(list(set(kt)-{achs[0],achs[1]}))
            key1 = s[0]+" "+s[1]
            key2 = s[1]+" "+s[0]
            D[key1]=-np.log(frq[k][0]/float(frq[k][1]))
            D[key2]=D[key1]
        return [D,frq]

def resolvePolytomy(pathToTree,node,verbose):
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
def findPolytomies(con_tree,taxa,anch):
    to_resolve = dict()
    maxPolyOrder = 0
    tmp = list()
    node_dict = dict()
    anch = sorted(anch)
    filter = lambda taxon: True if taxon.label==anch[0] else False
    nd1 = con_tree.find_node_with_taxon(filter)
    filter = lambda taxon: True if taxon.label==anch[1] else False
    nd2 = con_tree.find_node_with_taxon(filter)
    par1 = nd1.parent_node
    par2 = nd2.parent_node
    isSibiling = False
    if par1 == par2:
        isSibiling = True
        par1_is_Poly = True
        par2_is_Poly = True

    n1 = len(par1.adjacent_nodes())
    n2 = len(par2.adjacent_nodes())
    if ((n1>4 and not isSibiling )or (n1==4 and par1.parent_node is None and not isSibiling)):
        par1_is_Poly = True
    else:
        par1_is_Poly = False
    if ((n2>4 and not isSibiling )or (n2==4 and par2.parent_node is None and not isSibiling)):
        par2_is_Poly = True
    else:
        par2_is_Poly = False

    if not isSibiling:
        par1_child = {n.label for n in par1.leaf_iter()}
        par2_child = {n.label for n in par2.leaf_iter()}
    else:
        par1_child = {n.label for n in nd1.parent_node.leaf_iter()}
        par2_child = {n.label for n in nd2.parent_node.leaf_iter()}
    con_tree.prune_taxa([nd1.taxon,nd2.taxon], update_bipartitions=False, suppress_unifurcations=False)
    anch_nodes = {nd1,nd2}
    anch = set(anch)
    taxa = set(con_tree.leaf_nodes())
    for e in con_tree.postorder_node_iter():
        tmp_set = set()
        n = len(e.child_nodes())
        if n>3 or (n==3 and e.parent_node is not None and len(e.parent_node.child_nodes())>1):
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
            tmp.append(tmp_set)
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


    return (to_resolve,maxPolyOrder,(list(anch_nodes),con_tree.seed_node.label,par1,par2,par1_is_Poly,par2_is_Poly,par1_child,par2_child))

def addAnchores(con_tree,con_map):
    (anch,con_tree.seed_node.label,par1,par2,par1_is_Poly,par2_is_Poly,par1_child,par2_child)=con_map
    p1 = False
    p2 = False
    isSibiling = (par1_child == par2_child)
    nu = 0
    ach_a = set(anch)
    if anch[0].label>anch[1].label:
        t = anch[0]
        anch[0] = anch[1]
        anch[1] = t
    if len(par1.adjacent_nodes())>0 and (not par1_is_Poly) and (not isSibiling):
        par1.insert_child(0,anch[0])
        nu += 1
        ach_a -=  {anch[0]}
        p1 = True
    if (len(par2.adjacent_nodes())>0) and (not par2_is_Poly) and (not isSibiling):
        par2.insert_child(0,anch[1])
        p2 = True
        ach_a -=  {anch[1]}
        nu += 1
        if p1 or par1_is_Poly:
            print "adding anchores without traversing"
            return (nu,ach_a)
    if isSibiling and not par1_is_Poly:
        t=par1.insert_new_child(0)
        t.add_child(anch[0])
        t.add_child(anch[1])
        nu += 2
        ach_a = {}
        return (nu,ach_a)
    if (par1_is_Poly and par2_is_Poly) or (par2_is_Poly and p1):
        return (nu,ach_a)
    return (nu,ach_a)
def findPolytomies_with_name(con_tree,taxa,anch):
    to_resolve = dict()
    maxPolyOrder = 0
    tmp = list()
    node_dict = dict()
    anch = sorted(anch)
    filter = lambda taxon: True if taxon.label==anch[0] else False
    nd1 = con_tree.find_node_with_taxon(filter)
    filter = lambda taxon: True if taxon.label==anch[1] else False
    nd2 = con_tree.find_node_with_taxon(filter)
    par1 = nd1.parent_node
    par2 = nd2.parent_node
    isSibiling = False
    if par1 == par2:
        isSibiling = True
        par1_is_Poly = True
        par2_is_Poly = True

    n1 = len(par1.adjacent_nodes())
    n2 = len(par2.adjacent_nodes())
    if ((n1>4 and not isSibiling )or (n1==4 and par1.parent_node is None and not isSibiling)):
        par1_is_Poly = True
    else:
        par1_is_Poly = False
    if ((n2>4 and not isSibiling )or (n2==4 and par2.parent_node is None and not isSibiling)):
        par2_is_Poly = True
    else:
        par2_is_Poly = False

    if not isSibiling:
        par1_child = {n.label for n in par1.leaf_iter()}
        par2_child = {n.label for n in par2.leaf_iter()}
    else:
        par1_child = {n.label for n in nd1.parent_node.leaf_iter()}
        par2_child = {n.label for n in nd2.parent_node.leaf_iter()}
    con_tree.prune_taxa([nd1.taxon,nd2.taxon], update_bipartitions=False, suppress_unifurcations=False)
    anch_nodes = {nd1,nd2}
    anch = set(anch)
    taxa = set(con_tree.leaf_nodes())
    for e in con_tree.postorder_node_iter():
        tmp_set = set()
        n = len(e.child_nodes())
        if n>3 or (n==3 and e.parent_node is not None and len(e.parent_node.child_nodes())>1):
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
            tmp.append(tmp_set)
            node_dict[e.label] = tmp_set
            if len(taxa-tmp_set)>0:
                t = taxa-tmp_set
                v[e.parent_node.label] = t
            to_resolve[e.label] = v
        else:
            for i in range(0,n):
                if len(tmp)>0:
                    tmp_set = tmp_set | set(tmp.pop())
            if e.is_leaf():
                tmp_set.add(e)
            node_dict[e.label] = tmp_set
            tmp.append(tmp_set)


    return (to_resolve,maxPolyOrder,(list(anch_nodes),con_tree.seed_node.label,par1,par2,par1_is_Poly,par2_is_Poly,par1_child,par2_child))

def findIfPolytomyAnch(con_tree_tmp,anch,e):
    flag = 0
    return flag
def addDistanceAnchores(D1,D2,C1,C2):
    return
def anchoredDistanceFromFrqAddDistances(quartTable, anch,fillmethod):
    return [D,countElem]
