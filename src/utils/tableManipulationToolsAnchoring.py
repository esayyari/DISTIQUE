#! /usr/bin/env python
import itertools
from scipy import stats
from numpy import mean, sqrt, square
import numpy as np
import timer as tm
import random
import os
from numpy.lib.index_tricks import nd_grid
def generateKey(taxa_list):
    chosen = list()
    for  v in taxa_list.values():
        rt = random.sample(v,1)
        chosen.append(rt[0])
    allQuartetComb = itertools.combinations(chosen,4)
    origKeys = ['/'.join(sorted(q)) for q in allQuartetComb ]
    return origKeys


def findTrueAverageTableAnchoring(frq,anch,list_taxa,method,met):
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    TotalKey = dict()
    n = len(lst_taxa)
    numG = max(v[1] for v in frq.values())
    for i in range(0,n):
        for j in range(i+1,n):
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:
                    lab_taxon_i = taxon_i
                    lab_taxon_j = taxon_j
                    lab_taxon_k = anch[0]
                    lab_taxon_z = anch[1]
                    key_orig = "/".join(sorted([lab_taxon_i,lab_taxon_j,lab_taxon_k,lab_taxon_z]))

                    l = sorted([lst_taxa[i],lst_taxa[j],anch[0],anch[1]])
                    key_inv = "/".join(l)
                    if key_orig in frq:
                        v = frq[key_orig]
                    else:
                        v = list()
                        v.append(0.5)
                        v.append(numG)
                    v_inv = float(v[0])/v[1]
                    if key_inv in TotalKey:
                        if (met=="freq"):
                            vt = TotalKey[key_inv]
                            vt.append(v_inv)
                        elif met == "log":
                            vt = TotalKey[key_inv]
                            vt.append(-np.log(1.*v_inv))
                    else:
                        if (met == "freq"):

                            vt = list()
                            vt.append(v_inv)
                        elif met == "log":
                            vt = list()
                            vt.append(-np.log(1.*v_inv))
                    TotalKey[key_inv] = vt
    TotalKeyf = dict()
    for q,v2 in TotalKey.iteritems():
        l = set(q.split("/"))
        l = list(l - set(anch))
        if met == "log":
            if method == "gmean":
                vtt = np.exp(-stats.gmean(v2))
            elif method == "mean":
                vtt = np.exp(-mean(v2))
        else:
            vtt = np.exp(-sqrt(mean(square(v2))))
        if met == "freq":
            if method == "gmean":
                vtt = (stats.gmean(v2))
            elif method == "mean":
                vtt = (mean(v2))
            else:
                vtt = (sqrt(mean(square(v2))))
        TotalKeyf[q] = vtt
    return TotalKeyf

def genKey(p1,p2):
    if p1[0]<p2[0]:
            p = p1[0]+' '+p1[1]+' | '+p2[0]+' '+p2[1]
    else:
        p = p2[0]+' '+p2[1]+' | '+p1[0]+' '+p1[1]
    return p
def findTrueAverageTableAnchoringAddDistances(frq, anch, list_taxa, N,method, met):
    tm.tic()
    [TotalKeyf,_]=initializeQuartetTable( anch, list_taxa)
    anch = sorted(list(anch))
    lst_taxa = list_taxa.keys()
    TotalKey = dict()
    n = len(lst_taxa)
    numG = max(v[1] for v in frq.values())
    skipClades = N
    for i in range(0, n):
        if lst_taxa[i] in skipClades:
            continue
        for j in range(i+1, n):
            if lst_taxa[j] in skipClades:
                continue
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:
                    lab_taxon_i = taxon_i
                    lab_taxon_j = taxon_j
                    p = sorted([lab_taxon_i,lab_taxon_j])
                    key_orig = genKey(p,anch)
                    l = sorted([lst_taxa[i], lst_taxa[j], anch[0], anch[1]])
                    key_inv = "/".join(l)
                    if key_orig in frq:
                        v = frq[key_orig]
                    else:
                        v = list()
                        v.append(0.5)
                        v.append(numG)	
                    v_inv = float(v[0])/v[1]
                    if key_inv in TotalKey:
                        if met == "freq":
                            vt = TotalKey[key_inv]
                            vt.append(v_inv)
                        elif met == "log":
                            vt = TotalKey[key_inv]
                            vt.append(-np.log(1.*v_inv))
                    else:
                        if met == "freq":
                            vt = list()
                            vt.append(v_inv)
                        elif met == "log":
                            vt = list()
                            vt.append(-np.log(1.*v_inv))
                    TotalKey[key_inv] = vt
    for q, v2 in TotalKey.iteritems():
        if met == "log":
            if method == "gmean":
                vtt = np.exp(-stats.gmean(v2))
            elif method == "mean":
                vtt = np.exp(-mean(v2))
                
            else:
                vtt = np.exp(-sqrt(mean(square(v2))))
        elif met == "freq":
            if method == "gmean":
                vtt = (stats.gmean(v2))

            elif method == "mean":
                
                vtt = (mean(v2))
                
            else:
                vtt = (sqrt(mean(square(v2))))
                
        TotalKeyf[q] = vtt
    tm.toc()
    return TotalKeyf
def initializeQuartetTable( anch, list_taxa):
    D = dict()
    C = dict()
    lst_taxa = list(list_taxa.keys())
    n = len(lst_taxa)
    for i in range(0, n):
        for j in range(i+1, n):
                l = sorted([lst_taxa[i], lst_taxa[j], anch[0], anch[1]])
                key_inv = "/".join(l)
                D[key_inv] = -1.
                C[key_inv] = 0
    return [D,C]
def findTrueAverageTableAnchoringAddAnchores(quartTable,frq, anch, list_taxa, method, met):
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    n = len(lst_taxa)
    numG = max(v[1] for v in frq.values())
    for k_inv,v in list_taxa.iteritems():
        if anch[0] in v:
            N1 = k_inv
        if anch[1] in v:
            N2 = k_inv
    skipClades = {N1,N2}
    if N1 == N2:
        return
    for i in range(0, n):
        if lst_taxa[i] in skipClades:
            continue
        for j in range(i+1, n):
            if lst_taxa[j] in skipClades:
                continue
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:
                    lab_taxon_i = taxon_i
                    lab_taxon_j = taxon_j
                    oS = {taxon_i,taxon_j}
                    lab_taxon_k = anch[0]
                    lab_taxon_z = anch[1]
                    key_orig = "/".join(sorted([lab_taxon_i, lab_taxon_j, lab_taxon_k, lab_taxon_z]))
                    l = sorted([lst_taxa[i], lst_taxa[j]])
                    N = sorted([N1,N2])
                    key_inv = N[0]+" "+N[1]+"|"+l[0]+" "+l[1]
                    if anch[0] in oS or anch[1] in oS:
                        raise Exception("Oops! the anchors shouldn't be here!")
                    else:
                        if key_orig in frq:
                            v = frq[key_orig]
                        else:
                            v = list()
                            v.append(0.5)
                            v.append(numG)    
                        v_inv = float(v[0])/v[1]
                        if key_inv in quartTable:
                            if met == "freq":
                                vt = quartTable[key_inv]
                                vt.append(v_inv)
                            elif met == "log":
                                vt = quartTable[key_inv]
                                vt.append(-np.log(1.*v_inv))
                        else:
                            if met == "freq":
                                vt = list()
                                vt.append(v_inv)
                            elif met == "log":
                                vt = list()
                                vt.append(-np.log(1.*v_inv))
                        quartTable[key_inv] = vt
    return 
def findTrueAverageTableAnchoringOnDifferentSides(frq,anch,list_taxa,N1,N2,method, met):
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    TotalKey = dict()
    n = len(lst_taxa)
    N = {N1,N2}
    numG = max(v[1] for v in frq.values())
    for i in range(0,n):
        if lst_taxa[i] in N:
            continue 
        for j in range(i+1,n):
            if lst_taxa[j] in N:
                continue
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:

                    p = sorted([taxon_i,taxon_j])
                    key_orig = genKey(p,anch)

                    l = sorted([lst_taxa[i],lst_taxa[j],anch[0],anch[1]])
                    key_inv = "/".join(l)
                    if key_orig in frq:
                        v = frq[key_orig]
                    else:
                        v = list()
                        v.append(0.5)
                        v.append(numG)
                    v_inv = float(v[0])/v[1]
                    if key_inv in TotalKey:
                        if (met=="freq"):
                            vt = TotalKey[key_inv]
                            vt.append(v_inv)
                        elif met == "log":
                            vt = TotalKey[key_inv]
                            vt.append(-np.log(1.*v_inv))
                    else:
                        if (met == "freq"):
                            vt = list()
                            vt.append(v_inv)
                        elif met == "log":
                            vt = list()
                            vt.append(-np.log(1.*v_inv))
                    TotalKey[key_inv] = vt
    TotalKeyf = dict()
    for q,v2 in TotalKey.iteritems():
        l = set(q.split("/"))
        l = list(l - set(anch))
        if met == "log":
            if method == "gmean":
                vtt = np.exp(-stats.gmean(v2))
            elif method == "mean":
                vtt = np.exp(-mean(v2))
        else:
            vtt = np.exp(-sqrt(mean(square(v2))))
        if met == "freq":
            if method == "gmean":
                vtt = (stats.gmean(v2))
            elif method == "mean":
                vtt = (mean(v2))
            else:
                vtt = (sqrt(mean(square(v2))))
        TotalKeyf[q] = vtt
    return TotalKeyf
def findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomies(frq,TotalKeyf,anch,list_taxa,method,met):
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    TotalKey = dict()
    n = len(lst_taxa)
    numG = max(v[1] for v in frq.values())
    for k_inv,v in list_taxa.iteritems():
        if anch[0] in v:
            N1 = k_inv
        if anch[1] in v:
            N2 = k_inv
    N = {N1,N2}
    NL = sorted(list(N))
    for i in range(0,n):
        if lst_taxa[i] in N:
            continue
        for j in range(i+1,n):
            if lst_taxa[j] in N:
                continue
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:
                    p=sorted([taxon_i,taxon_j])
                    key_orig = genKey(p, anch)

                    l = sorted([lst_taxa[i],lst_taxa[j],N1,N2])
                    key_inv = "/".join(l)
                    if key_orig in frq:
                        v = frq[key_orig]
                        v_inv = float(v[0])/v[1]
                    else:
                        v_inv = 0.5/numG
                    if key_inv in TotalKey:
                        if (met=="freq"):
                            vt = TotalKey[key_inv]
                            vt.append(v_inv)
                        elif met == "log":
                            vt = TotalKey[key_inv]
                            vt.append(-np.log(1.*v_inv))
                    else:
                        if (met == "freq"):
                            vt = list()
                            vt.append(v_inv)
                        elif met == "log":
                            vt = list()
                            vt.append(-np.log(1.*v_inv))
                    TotalKey[key_inv] = vt
    for q,v2 in TotalKey.iteritems():
        l = sorted(list(q.split("/")))
        O = sorted(list(set(l)-N))
        if l[0] == NL[0]:
            key = NL[1]
        elif l[0] == O[0]:
            key = O[1]
        if q in TotalKeyf:
            vtmp=TotalKeyf[q]
            if key in vtmp:
                for x in v2:
                    vtmp[key].append(x)
            else:
                vtmp[key] = v2
        else:
            vtmp = dict()
            vtmp[key] = v2 
            TotalKeyf[q] = vtmp
    return TotalKeyf
def readTable(tmpPath):
    out_path  = os.path.expanduser(os.path.expandvars(tmpPath))
    f = open(out_path, 'r')
    frq = dict()
    for line in f:
        k=line.split()
        v = dict()
        d = k[0].split('/')
        s = float(k[1])+float(k[2])+float(k[3])
        v[d[1]] = [float(k[1]),s]
        v[d[2]] = [float(k[2]),s]
        v[d[3]] = [float(k[3]),s]
        frq[k[0]] = v
    return frq
def findTrueAverageTableAnchoringAddDistancesOverall(frq, anch, list_taxa, N,method, met):
    tm.tic()
    [TotalKeyf,_]=initializeQuartetTable( anch, list_taxa)
    anch = sorted(list(anch))
    lst_taxa = list_taxa.keys()
    TotalKey = dict()

    n = len(lst_taxa)
    skipClades = N
    for i in range(0, n):
        if lst_taxa[i] in skipClades:
            continue
        for j in range(i+1, n):
            if lst_taxa[j] in skipClades:
                continue
            l = sorted([lst_taxa[i], lst_taxa[j], anch[0], anch[1]])
            key_inv = "/".join(l)
            key_orig = genKey(anch,sorted([lst_taxa[i],lst_taxa[j]]))
            v = frq[key_orig]
            if len(v) == 1:
                v.append(1)   
            if key_inv in TotalKey:
                    vt = TotalKey[key_inv]
                    vt[0] += v[0]
                    vt[1] += v[1]
            else:
                    vt = list()
                    vt = v
            TotalKey[key_inv] = vt
    for q, v2 in TotalKey.iteritems():
        vtt = v2[0]/v2[1]
        TotalKeyf[q] = vtt
    tm.toc()
    return TotalKeyf
def findTrueAverageTableAnchoringAddDistancesOverallFromFile(frq, anch, list_taxa, N,method, met):
    tm.tic()
    [TotalKeyf,_]=initializeQuartetTable( anch, list_taxa)
    anch = sorted(list(anch))
    lst_taxa = list_taxa.keys()
    TotalKey = dict()
    n = len(lst_taxa)
    skipClades = N
    for i in range(0, n):
        if lst_taxa[i] in skipClades:
            continue
        for j in range(i+1, n):
            if lst_taxa[j] in skipClades:
                continue
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:
                    lab_taxon_i = taxon_i
                    lab_taxon_j = taxon_j
                    p = sorted([lab_taxon_i,lab_taxon_j])
                    key_orig = genKey(p,anch)
                    l = sorted([lst_taxa[i], lst_taxa[j], anch[0], anch[1]])
                    key_inv = "/".join(l)        
                    v = frq[key_orig]
                    if len(v) == 1:
                        v.append(1)
                    else:
                        v[0] -= 0.5
                        v[1] -= 1.5   
                    if key_inv in TotalKey:
                            vt = TotalKey[key_inv]
                            vt[0] += v[0]
                            vt[1] += v[1]
                    else:
                            vt = list()
                            vt = v
                    TotalKey[key_inv] = vt
    for q, v2 in TotalKey.iteritems():
        vtt = (v2[0]+0.5)/(v2[1]+1.5)
        TotalKeyf[q] = vtt
    tm.toc()
    return TotalKeyf
def findTrueAverageTableAnchoringOnDifferentSidesOverall(frq,anch,list_taxa,N1,N2,method, met):
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    TotalKey = dict()
    n = len(lst_taxa)
    N = {N1,N2}
    for i in range(0,n):
        if lst_taxa[i] in N:
            continue 
        for j in range(i+1,n):
            if lst_taxa[j] in N:
                continue
            p = sorted([lst_taxa[i],lst_taxa[j]])
            key_orig = genKey(p,anch)

            l = sorted([lst_taxa[i],lst_taxa[j],anch[0],anch[1]])
            key_inv = "/".join(l)
            v = frq[key_orig]
            v_inv = float(v[0])/v[1]
            if key_inv in TotalKey:
                if (met=="freq"):
                    vt = TotalKey[key_inv]
                    vt.append(v_inv)
                elif met == "log":
                    vt = TotalKey[key_inv]
                    vt.append(-np.log(1.*v_inv))
            else:
                if (met == "freq"):
                    vt = list()
                    vt.append(v_inv)
                elif met == "log":
                    vt = list()
                    vt.append(-np.log(1.*v_inv))
            TotalKey[key_inv] = vt
    TotalKeyf = dict()
    for q,v2 in TotalKey.iteritems():
        l = set(q.split("/"))
        l = list(l - set(anch))
        if met == "log":
            if method == "gmean":
                vtt = np.exp(-stats.gmean(v2))
            elif method == "mean":
                vtt = np.exp(-mean(v2))
        else:
            vtt = np.exp(-sqrt(mean(square(v2))))
        if met == "freq":
            if method == "gmean":
                vtt = (stats.gmean(v2))
            elif method == "mean":
                vtt = (mean(v2))
            else:
                vtt = (sqrt(mean(square(v2))))
        TotalKeyf[q] = vtt
    return TotalKeyf
def findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomiesOverall(frq,TotalKeyf,anch,list_taxa,N1,N2,method,met):
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    TotalKey = dict()
    n = len(lst_taxa)
    
    N = {N1,N2}
    NL = sorted(list(N))
    for i in range(0,n):
        if lst_taxa[i] in N:
            continue
        for j in range(i+1,n):
            if lst_taxa[j] in N:
                continue
            p=sorted([lst_taxa[i],lst_taxa[j]])
            key_orig = genKey(p, anch)

            l = sorted([lst_taxa[i],lst_taxa[j],N1,N2])
            key_inv = "/".join(l)
            v = frq[key_orig]
            v_inv = float(v[0])/v[1]
            if key_inv in TotalKey:
                if (met=="freq"):
                    vt = TotalKey[key_inv]
                    vt.append(v_inv)
                elif met == "log":
                    vt = TotalKey[key_inv]
                    vt.append(-np.log(1.*v_inv))
            else:
                if (met == "freq"):
                    vt = list()
                    vt.append(v_inv)
                elif met == "log":
                    vt = list()
                    vt.append(-np.log(1.*v_inv))
            TotalKey[key_inv] = vt
    for q,v2 in TotalKey.iteritems():
        l = sorted(list(q.split("/")))
        O = sorted(list(set(l)-N))
        if l[0] == NL[0]:
            key = NL[1]
        elif l[0] == O[0]:
            key = O[1]
        if q in TotalKeyf:
            vtmp=TotalKeyf[q]
            if key in vtmp:
                for x in v2:
                    vtmp[key].append(x)
            else:
                vtmp[key] = v2
        else:
            vtmp = dict()
            vtmp[key] = v2 
            TotalKeyf[q] = vtmp
    return TotalKeyf

    
def findTrueAverageTableAnchoringOnDifferentSidesOverallFromFile(frq,anch,list_taxa,N1,N2,method, met):
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    TotalKey = dict()
    n = len(lst_taxa)
    N = {N1,N2}
    for i in range(0,n):
        if lst_taxa[i] in N:
            continue 
        for j in range(i+1,n):
            if lst_taxa[j] in N:
                continue
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:
                    lab_taxon_i = taxon_i
                    lab_taxon_j = taxon_j
                    p = sorted([lab_taxon_i,lab_taxon_j])
                    key_orig = genKey(p,anch)
                    l = sorted([lst_taxa[i], lst_taxa[j], anch[0], anch[1]])
                    key_inv = "/".join(l)        
                    v = frq[key_orig]
                    if len(v) == 1:
                        v.append(1)
                    else:
                        v[0] -= 0.5
                        v[1] -= 1.5   
                    if key_inv in TotalKey:
                            vt = TotalKey[key_inv]
                            vt[0] += v[0]
                            vt[1] += v[1]
                    else:
                            vt = list()
                            vt = v
                    TotalKey[key_inv] = vt
    TotalKeyf = dict()
            
        
    for q,v2 in TotalKey.iteritems():
        l = set(q.split("/"))
        l = list(l - set(anch))
        
        TotalKeyf[q] = (v2[0]+0.5)/(v2[1]+1.5)
    return TotalKeyf


def findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomiesOverallFromFile(frq,TotalKeyf,anch,list_taxa,N1,N2,method,met):
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    TotalKey = dict()
    n = len(lst_taxa)
    
    N = {N1,N2}
    NL = sorted(list(N))
    for i in range(0,n):
        if lst_taxa[i] in N:
            continue
        for j in range(i+1,n):
            if lst_taxa[j] in N:
                continue
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:
                    lab_taxon_i = taxon_i
                    lab_taxon_j = taxon_j
                    p = sorted([lab_taxon_i,lab_taxon_j])
                    key_orig = genKey(p,anch)
                    l = sorted([lst_taxa[i], lst_taxa[j],N1,N2])
                    key_inv = "/".join(l)        
                    v = frq[key_orig]
                    if len(v) == 1:
                        v.append(1)
                    else:
                        v[0] -= 0.5
                        v[1] -= 1.5   
                    if key_inv in TotalKey:
                            vt = TotalKey[key_inv]
                            vt[0] += v[0]
                            vt[1] += v[1]
                    else:
                            vt = list()
                            vt = v
                    TotalKey[key_inv] = vt
                    

    for q,v2 in TotalKey.iteritems():
        l = sorted(list(q.split("/")))
        O = sorted(list(set(l)-N))
        vtt = (v2[0]+0.5)/(v2[1]+1.5)
        if l[0] == NL[0]:
            key = NL[1]
        elif l[0] == O[0]:
            key = O[1]
        if q in TotalKeyf:
            vtmp=TotalKeyf[q]
            if key in vtmp:
                vtmp[key].append(vtt)
            else:
                vtmp[key] = list()
                vtmp[key].append(vtt)
        else:
            vtmp = dict()
            vtmp[key] = list()
            vtmp[key].append(vtt) 
        TotalKeyf[q] = vtmp
    return TotalKeyf
