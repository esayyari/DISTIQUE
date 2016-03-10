#! /usr/bin/env python
import itertools
from scipy import stats
from numpy import mean, sqrt, square
import numpy as np
import timer as tm
import random
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
                    if (taxon_i in set(anch) or taxon_j in set(anch)):
                        continue
                    else:
                        if key_orig in frq:
                            v = frq[key_orig]
                        else:
                            v[0] = 0.5
                            v[1] = numG
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


def findTrueAverageTableAnchoringAddDistances(frq, anch, list_taxa, method, met):
    tm.tic()
    [TotalKeyf,_]=initializeQuartetTable( anch, list_taxa)
    anch = sorted(list(anch))
    lst_taxa = list(list_taxa.keys())
    TotalKey = dict()
    n = len(lst_taxa)
    numG = max(v[1] for v in frq.values())
    for i in range(0, n):
        for j in range(i+1, n):
            for taxon_i in list_taxa[lst_taxa[i]]:
                for taxon_j in list_taxa[lst_taxa[j]]:
                    lab_taxon_i = taxon_i
                    lab_taxon_j = taxon_j
                    lab_taxon_k = anch[0]
                    lab_taxon_z = anch[1]
                    key_orig = "/".join(sorted([lab_taxon_i, lab_taxon_j, lab_taxon_k, lab_taxon_z]))
                    l = sorted([lst_taxa[i], lst_taxa[j], anch[0], anch[1]])
                    key_inv = "/".join(l)
                    if (key_orig.split("/")).count(anch[0]) > 1 or (key_orig.split("/")).count(anch[1]) > 1:
                        continue
                    else:
                        if key_orig in frq:
                            v = frq[key_orig]
                        else:
                                v[0] = 0.5
                                v[1] = numG	
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

