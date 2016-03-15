import dendropy
import itertools
import sys
import os
import tableManipulationTools as tbs
import printTools as pr
import toolsTreeTaxa as tstt
import numpy as np
import random
import warnings
import timer as tm
from scipy import stats
from numpy import mean, sqrt, square
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
    for k in frq:
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
#         print "time for re-rooting is: "
#         tm.tic()
        rerooted=reroot(tree,anch)
#         tm.toc()
        node = rerooted[0]
        root = rerooted[1]
#         print "time for counting is:"
#         tm.tic()
        while(node.parent_node is not root):
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
#         tm.toc()
    return Q
def reroot(tree,anch):
    #tstt.labelNodes(tree)
    filt = lambda node: True if node.taxon.label  in anch else False
    nodes = [x for x in tree.leaf_node_iter(filt)]
    if tree.seed_node is not None and tree.seed_node.taxon is not None:
        if tree.seed_node.taxon.label in anch:
            if nodes[0] is not tree.seed_node:
                return [nodes[0],tree.seed_node]
            else:
                raise Exception("Oops! there is not leaf with the same label as one of the anchors!")
    node = nodes[1]
    root = nodes[0]
#     print "This is the main rerooting!"
#     tm.tic()

    tree.reroot_at_node(root, update_bipartitions=False, suppress_unifurcations=False)
#     tm.toc()
    return [node,root]
def findAllChildrenPairs(listTaxa):
    listTaxaLabels = list()
    for t in listTaxa:
        listTaxaLabels.append(t.taxon.label)
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
            T.append(t.taxon.label)
        x.append(T)
    keySet=[p for p in itertools.product(*x)]
    for p in keySet:
        key = sorted(list(p) + anch)
        Q["/".join(key)][1] -= 1
    return Q
def findAnchoredDistanceTable(achs,trees,taxa,out):
    frq=findAnchoredQuartets(achs,trees,taxa,out)
    D = dict()
    for k in frq:
        kt = sorted(k.split("/"))
        if ((achs[0] in kt ) and( achs[1] in kt)):
            s = sorted(list(set(kt)-{achs[0],achs[1]}))
            key1 = s[0]+" "+s[1]
            key2 = s[1]+" "+s[0]
            D[key1]=-np.log(frq[k][0]/float(frq[k][1]))
            D[key2]=D[key1]
        return [D,frq]
def findAnchoredDistanceTableFromFile(anch,frqT,taxa,outpath):
    frq=findAnchoredQuartetsFromFile(anch,frqT,taxa,outpath)
    return frq
def findAnchoredQuartetsFromFile(anch,frqT,taxa,out):
    anch = sorted(anch)
    anchS = set(anch)
    a=frqT.keys()
    b=frqT[a[0]].keys()
    n = frqT[a[0]][b[0]][1]
    Q = buildEmptyQuartets(anch,taxa,n)
    for key in Q:
        l = key.split("/")
        anchComp = sorted(list(set(l) - anchS))
        if l[0] in anchS:
            tq = anch[1]
        else:
            tq = anchComp[1]
        Q[key][0] = frqT[key][tq][0]
        Q[key][1] = frqT[key][tq][1]
    return Q
def resolvePolytomy(pathToTree,node,verbose):
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
            filt = lambda node: True if node.taxon.label==label else False
            if sp_tree.seed_node.taxon is not None and sp_tree.seed_node.taxon.label == label:
                nd = sp_tree.seed_node
            else:
                ndl= [x for x in sp_tree.leaf_node_iter(filt)]
                nd = ndl[0]
            if nd is not None:
                sp_tree.reroot_at_edge(nd.edge, update_bipartitions=False)

        stack = list()
        for e in sp_tree.postorder_node_iter():
                n = len(e.child_nodes())
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
def findPolytomies(con_tree,taxa,anch):
    to_resolve = dict()
    maxPolyOrder = 0
    tmp = list()
    node_dict = dict()
    anch = sorted(anch)
    filt = lambda taxon: True if taxon.label in anch else False
    nodes = [x for x in con_tree.leaf_node_iter(filt)]
    if ( con_tree.seed_node.taxon is not None and con_tree.seed_node.taxon.label in anch ):
        nd1 = con_tree.seed_node
        nd2 = nodes[0]
    else:
        nd1 = nodes[0]
        nd2 = nodes[1]
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
        par1_child = {n.taxon.label for n in par1.leaf_iter()}
        par2_child = {n.taxon.label for n in par2.leaf_iter()}
    else:
        par1_child = {n.taxon.label for n in nd1.parent_node.leaf_iter()}
        par2_child = {n.taxon.label for n in nd2.parent_node.leaf_iter()}
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
                    if c.taxon is not None:
                        v[c.taxon.label] = node_dict[c.taxon.label]
                    else:
                        v[c.label] = node_dict[c.label]
                    tmp_set = tmp_set | t
            tmp.append(tmp_set)
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
            to_resolve[e] = v
        else:
            for _ in range(0,n):
                if len(tmp)>0:
                    tmp_set = tmp_set | set(tmp.pop())
            if e.is_leaf():
                tmp_set.add(e)
                node_dict[e.taxon.label] = tmp_set
            else:
                node_dict[e.label] = tmp_set
            tmp.append(tmp_set)

    if con_tree.seed_node.taxon is not None:
        seedLab = con_tree.seed_node.taxon.label
    else:
        seedLab = con_tree.seed_node.label
    return (to_resolve,maxPolyOrder,(list(anch_nodes),seedLab,par1,par2,par1_is_Poly,par2_is_Poly,par1_child,par2_child))

def addAnchores(con_tree,con_map):
    (anch,seedLabel,par1,par2,par1_is_Poly,par2_is_Poly,par1_child,par2_child)=con_map
    if con_tree.seed_node.taxon is not None:
        con_tree.seed_node.taxon.label = seedLabel
    else:
        con_tree.seed_node.label = seedLabel
    p1 = False
    isSibiling = (par1_child == par2_child)
    nu = 0
    ach_a = set(anch)
    if anch[0].taxon.label>anch[1].taxon.label:
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
    filt = lambda taxon: True if taxon.label==anch[0] else False
    nd1 = con_tree.find_node_with_taxon(filt)
    filt = lambda taxon: True if taxon.label==anch[1] else False
    nd2 = con_tree.find_node_with_taxon(filt)
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
        par1_child = {n.taxon.label for n in par1.leaf_iter()}
        par2_child = {n.taxon.label for n in par2.leaf_iter()}
    else:
        par1_child = {n.taxon.label for n in nd1.parent_node.leaf_iter()}
        par2_child = {n.taxon.label for n in nd2.parent_node.leaf_iter()}
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

                    v[c.taxon.label] = node_dict[c.taxon.label]
                    tmp_set = tmp_set | t
            tmp.append(tmp_set)
            node_dict[e.taxon.label] = tmp_set
            if len(taxa-tmp_set)>0:
                t = taxa-tmp_set
                v[e.parent_node.taxon.label] = t
            to_resolve[e.taxon.label] = v
        else:
            for _ in range(0,n):
                if len(tmp)>0:
                    tmp_set = tmp_set | set(tmp.pop())
            if e.is_leaf():
                tmp_set.add(e)
            node_dict[e.taxon.label] = tmp_set
            tmp.append(tmp_set)


    return (to_resolve,maxPolyOrder,(list(anch_nodes),con_tree.seed_node.taxon.label,par1,par2,par1_is_Poly,par2_is_Poly,par1_child,par2_child))


def addDistanceAnchores(D,Dtmp,C,Ctmp,fillmethod):
    if fillmethod == "normConst":
        Dmax = max(np.abs(Dtmp.values()))*1.
        if Dmax == 0:
            Dmax = 1.
        for kttDtmp in Dtmp:
            Dtmp[kttDtmp] /= Dmax
    for key in D:
        C[key] += Ctmp[key]
        D[key] += Dtmp[key]
    return
def anchoredDistanceFromFrqAddDistances(frq,achs,taxa_list):
    D = dict()
    C = dict()
    anch_set = {achs[0],achs[1]}
    for k in frq:
        kt = sorted(k.split("/"))
        s = sorted(list(set(kt)-anch_set))
        if len(s) == 1:
            c0 = kt.count(achs[0])
            c1 = kt.count(achs[1])
            if c0 == 2:
                s.append(achs[0])
            elif c1 == 2:
                s.append(achs[1])
            key1 = s[0] + " " + s[1]
            key2 = s[1] + " " + s[0]
            D[key1] = -np.inf
            D[key2] = D[key1]
            C[key1] = 0
            C[key2] = C[key1]
        elif len(s) == 0 and len(kt) > 0:
            s = list()
            s.append(achs[0])
            s.append(achs[1])
            key1 = s[0] + " "+s[1]
            key2 = s[1] + " "+s[0]
            D[key1] = -np.inf
            D[key2] = D[key1]
            C[key1] = 0
            C[key2] = C[key1]
        else:
            key1 = s[0]+" "+s[1]
            key2 = s[1]+" "+s[0]
            if frq[k]<0:
                D[key1] = -np.inf
                D[key2] = D[key1]
                C[key1] = 0
                C[key2] = C[key1]
            else:
                D[key1] = -np.log(frq[k])
                D[key2] = D[key1]
                C[key1] = 1
                C[key2] = C[key1]
    return [D,C]
def fillEmptyElementsDistanceTable(D,C,fillmethod):
    random.seed()
    for key in D:
        if D[key] == -np.inf:
            if fillmethod == "rand":
                l = -random.uniform(-1,0)
                if l>1./3.:
                    D[key] = -np.log(3./2.*(1-l))
                else:
                    D[key] = -np.log(3.*l)
                C[key] = 1
            else:
                D[key] = 0
                C[key] = 0
    return
def normalizeDistanceTable(D,C):
    flag = any(C[key] == 0 or D[key]<0 for key in D)
    if not flag:
        for key in D:
            D[key] += 10.
            D[key] /= C[key]
    return flag

def check_four_point(D,listTaxa):
    Dv = list()
    sT = {0,1,2,3}
    for i in range(0,len(listTaxa)):
        l = sorted(listTaxa[i])
        d = dict()
        for j in range(1,4):
            s = sorted(list(sT-{0,j}))
            key2 = l[s[0]] + " " + l[s[1]]
            key1 = l[0] + " " + l[j]
            d[j] = D[key2] + D[key1]
        Dv.append(max(d))
    return Dv
def findDistanceTable(Q,method,met):
    D = dict()
    for q, v2 in Q.iteritems():
        if met == "log":
            if method == "gmean":
                vtt = -np.log(np.exp(-stats.gmean(v2)))
            elif method == "mean":
                vtt = -np.log(np.exp(-mean(v2)))
            else:
                vtt = -np.log(np.exp(-sqrt(mean(square(v2)))))
        elif met == "freq":
            if method == "gmean":
                vtt = -np.log(stats.gmean(v2))
            elif method == "mean":
                vtt = -np.log(mean(v2))
            else:
                vtt = -np.log(sqrt(mean(square(v2))))
        D[q] = vtt        
    return D
def sumofLogsDistanceTable(D):
    D_Final = dict()
    for k in D:
        l = k.split("|")
        k1 = l[0].split(" ")
        k2 = l[1].split(" ")
        
        if l[0] in D_Final:
            D_Final[k1[0]+" "+k1[1]] += D[k]
            D_Final[k1[1]+" "+k1[0]] += D[k]
        else:
            D_Final[k1[0]+" "+k1[1]] = D[k]
            D_Final[k1[1]+" "+k1[0]] = D[k]
        if l[1] in D_Final:
            D_Final[k2[0]+" "+k2[1]] += D[k]
            D_Final[k2[1]+" "+k2[0]] += D[k]
        else:
            D_Final[k2[0]+" "+k2[1]] = D[k]
            D_Final[k2[1]+" "+k2[0]] = D[k]
    return D_Final
def mapTaxaAroundPoly(to_resolve,debugFlag):
    mapTaxaToPolyNodes = dict()
    mapPolyNodesToTaxa = dict()
    for key in to_resolve:
        v = to_resolve[key]
        for nodes in v:
            for taxon in v[nodes]:
                if taxon.taxon.label in mapTaxaToPolyNodes:
                    mapTaxaToPolyNodes[taxon.taxon.label].add(nodes)
                else:
                    mapTaxaToPolyNodes[taxon.taxon.label] = {nodes}
            mapPolyNodesToTaxa[nodes] = {x.taxon.label for x in v[nodes]}
    if debugFlag:
        print "the mapping from taxa to Poly nodes"
        for key in mapTaxaToPolyNodes:
            print key+" of length: "+str(len(mapTaxaToPolyNodes[key]))+" maps to: ",mapTaxaToPolyNodes[key]
        print "the mapping from Poly nodes to taxa"
        for key in mapPolyNodesToTaxa:
            print key+" of lenght: "+str(len(mapPolyNodesToTaxa[key]))+" maps to: ",mapPolyNodesToTaxa[key]
    return [mapTaxaToPolyNodes,mapPolyNodesToTaxa]
            
