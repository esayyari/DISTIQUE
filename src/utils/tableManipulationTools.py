#! /usr/bin/env python
import dendropy
import itertools
import os
import printTools as pr
from scipy import stats
from numpy import mean, sqrt, square
from  prodDistance import prodDistance
from minDistance import minDistance
import copy
import numpy as np
import tempfile
import random


def convergencedQuartTable(qTable1,qTable2,eps,verbose):
	kl = KLdivergence(qTable1,qTable2)
	if kl<eps:
		if verbose==1:
			print kl
		return True
	else:
		return False
	

def KLdivergence(qTable1,qTable2):
	kld = dict()
	for key in qTable1:
		v1 = qTable1[key]
		sz1 = sum(v1.values())
		v2 = qTable2[key]
		sz2 = sum(v2.values())
		kld[key] = 0
		for key2 in v1:
			kld[key] += float(v1[key2])/float(sz1)*np.log((float(v1[key2])*sz2)/(sz1*float(v2[key2])))
	totalKL = sum(kld.values())
	return totalKL
def expandDict(frq):
	expandedDict = dict()
	s = {0,1,2,3}
	for key in frq:
		l = sorted(key.split("/"))
		t = dict()
		for i in range(1,4):
			t[l[0]+l[i]] = frq[key][l[i]]
			t[l[i]+l[0]] = frq[key][l[i]]
			v = list(s - {0,i})
			t[l[v[0]]+l[v[1]]] = frq[key][l[i]]			
			t[l[v[1]]+l[v[0]]] = frq[key][l[i]]
		expandedDict[key] = t
	return expandedDict

def partialQuartetTable(quartTable,origKeys,inv_taxa):	
	expandedQTable = expandDict(quartTable)
	pQuartTable = dict()
	inv_keys = dict()
	for key in origKeys:
		tmp = key.split("/")
		inv_tmp = list()
		for k in tmp:
			it = inv_taxa[k]
			inv_tmp.append(it)
			inv_keys[it] = k
		inv_tmp_sorted = sorted(inv_tmp)
		dummyKey = "/".join(inv_tmp_sorted)
		t = dict()
		for i in range(1,4):
			t[inv_tmp_sorted[i]] = expandedQTable[key][inv_keys[inv_tmp_sorted[0]]+inv_keys[inv_tmp_sorted[i]]]	
		pQuartTable[dummyKey] = t
	return pQuartTable
def findTrueAverageTable(frq,list_taxa,method,met):
	n = len(list_taxa)
# 	print "n is: " + str(n)
	lst_taxa = list(list_taxa.keys())
	TotalKey = dict()
	s = {1,2,3}
	for i in range(0,n):
		for j in range(i+1,n):
			for k in range(j+1,n):
				for z in range(k+1,n):
					for taxon_i in list_taxa[lst_taxa[i]]:
						for taxon_j in list_taxa[lst_taxa[j]]:
							for taxon_k in list_taxa[lst_taxa[k]]:
								for taxon_z in list_taxa[lst_taxa[z]]:
									keyt = "/".join(sorted([taxon_i,taxon_j,taxon_k,taxon_z]))
									lab_taxon_i = taxon_i
									lab_taxon_j = taxon_j
									lab_taxon_k = taxon_k
									lab_taxon_z = taxon_z
									tmp_dict = dict()
									tmp_dict[lst_taxa[i]] = lab_taxon_i
									tmp_dict[lst_taxa[j]] = lab_taxon_j
									tmp_dict[lst_taxa[k]] = lab_taxon_k
									tmp_dict[lst_taxa[z]] = lab_taxon_z
									key_orig = "/".join(sorted([lab_taxon_i,lab_taxon_j,lab_taxon_k,lab_taxon_z]))
									l = sorted([lst_taxa[i],lst_taxa[j],lst_taxa[k],lst_taxa[z]])
									key_inv = "/".join(l)
									v = frq[key_orig]
									v_inv = dict()
									sz = 0
									for q in range(1,4):
										q1 = sorted([tmp_dict[l[0]],tmp_dict[l[q]]])
										stmp = list(s-{q})
										q2 = sorted([tmp_dict[l[stmp[0]]],tmp_dict[l[stmp[1]]]])
										if q1[0]<q2[0]:
											v_inv[l[q]] = v[q1[1]]	
											sz +=  v[q1[1]]
										else:
											v_inv[l[q]] = v[q2[1]]
											sz += v[q2[1]]
									if key_inv in TotalKey:
										vt = TotalKey[key_inv] 
										for keyt in vt.keys():
											if met == "freq":
												vt[keyt].append(float(v_inv[keyt])/sz)
											elif met == "log":
												prob = float(v_inv[keyt])/sz
												if prob >=1.:
													prob = 1. - 10.**(-10)
													vt[keyt].append(-np.log(prob))
									else:
										vt = dict()
										for q in v_inv:
											if met == "freq":
	
												vt[q] = list()
												vt[q].append(float(v_inv[q])/sz)
											elif met == "log":
												prob = float(v_inv[q])/sz
												if prob >= 1.:
													prob = 1. - 10.**(-10)
												vt[q] = list()
												vt[q].append(-np.log(prob))
									TotalKey[key_inv] = vt
									
	TotalKeyf = dict()
	for q,v in TotalKey.iteritems():
		vtt = dict()
		for q2,v2 in v.iteritems():
			if met == "log":
				if method == "gmean":
					vtt[q2] = np.exp(-stats.gmean(v2))
				elif method == "mean":
					b = np.exp(-mean(v2))
					vtt[q2] = (b)
				else:
					vtt[q2] = np.exp(-sqrt(mean(square(v2))))
			if met == "freq":
				if method == "gmean":
					vtt[q2] = (stats.gmean(v2))
				elif method == "mean":
					vtt[q2] = (mean(v2))
				else:
					vtt[q2] = (sqrt(mean(square(v2))))
		TotalKeyf[q] = vtt
	return TotalKeyf
										
									

def normalizeDistanceTable(D): 
    Max = 1.
    for key in D:
        Max = np.max([D[key],Max])
    for key in D:
    	D[key] += Max
        D[key] /= Max
    return
def distanceTable(frq,method,outfile,met):
	keyDict = sorted(np.unique(("/".join(frq.keys())).split("/")));
	mapDict = dict()       	
	if method == 'min':
		mapDict =minDistance(frq,met)
	elif method == 'prod':
		mapDict =prodDistance(frq,met)
	pr.printDistanceTableToFile(mapDict,keyDict,outfile)


	
def generateKey(taxa_list):
	chosen = list()
	for v in taxa_list.values():
		rt = random.sample(v,1)
		chosen.append(rt[0])
	allQuartetComb = itertools.combinations(chosen,4)
	origKeys = ['/'.join(sorted(q)) for q in allQuartetComb ]
	return origKeys
		
		
		
def readQuartetTable(gt):
	src_fpath = os.path.expanduser(os.path.expandvars(gt))
	frq = readTable(src_fpath)
	for k in frq:
		sz = sum(frq[k].values())
		for k2 in frq[k]:
			frq[k][k2] /= sz
	return frq
def readTable(tmpPath):
	out_path  = os.path.expanduser(os.path.expandvars(tmpPath))
	f = open(out_path, 'r')
	frq = dict()
	for line in f:
		k=line.split()
		v = dict()
		d = k[0].split('/')
		v[d[1]] = float(k[1])
		v[d[2]] = float(k[2])
		v[d[3]] = float(k[3])
		frq[k[0]] = v
	return frq

		
#!/usr/bin/env python

def findQuartetTable(trees,origKeys,keyType,tmpPath,verbose):
	
	WS_LOC_SH = os.environ['WS_HOME']+'/DISTIQUE/src/shell/'
	out_path  = os.path.expanduser(os.path.expandvars(tmpPath))
	taxa = set()
	if keyType==1:
		for key in origKeys:	
			l = key.split("/")
			taxa = taxa | set(l)
	else:
		taxa = origKeys
	treeList = dendropy.TreeList()
	for tree in trees:
		t = copy.deepcopy(tree)
		t.retain_taxa_with_labels(taxa, update_bipartitions=False, suppress_unifurcations=True)
		treeList.append(t)
	outfile = "partialQuartetTable.nwk"
	ftmp = out_path+"/"+outfile
	treeList.write(path = ftmp,schema="newick",suppress_rooting=True)
	command = WS_LOC_SH+"compute.quartet.freq.sh -o "+out_path+" -w 1 -f "+ftmp

	os.system(command)
	frq = 	readTable(out_path+'/quartet_tmp.q')
	if verbose:
		print "The quartet table has been computed"
	return frq
def findTrueAverageTableAnchoring(frq,anch,list_taxa,method):
	n = len(set(list_taxa)-set(anch))	
	anch = sorted(list(anch))
	lst_taxa = list(list_taxa.keys())
	TotalKey = dict()
	s = {1,2,3}
	for i in range(0,n):
		for j in range(i+1,n):
			for k in range(j+1,n):
				for z in range(k+1,n):
					for taxon_i in list_taxa[lst_taxa[i]]:
						for taxon_j in list_taxa[lst_taxa[j]]:
							for taxon_k in list_taxa[lst_taxa[k]]:
								for taxon_z in list_taxa[lst_taxa[z]]:
									keyt = "/".join(sorted([taxon_i,taxon_j,taxon_k,taxon_z]))
									lab_taxon_i = taxon_i
									lab_taxon_j = taxon_j
									lab_taxon_k = taxon_k
									lab_taxon_z = taxon_z
									tmp_dict = dict()
									tmp_dict[lst_taxa[i]] = lab_taxon_i
									tmp_dict[lst_taxa[j]] = lab_taxon_j
									tmp_dict[lst_taxa[k]] = lab_taxon_k
									tmp_dict[lst_taxa[z]] = lab_taxon_z
									key_orig = "/".join(sorted([lab_taxon_i,lab_taxon_j,lab_taxon_k,lab_taxon_z]))
									l = sorted([lst_taxa[i],lst_taxa[j],lst_taxa[k],lst_taxa[z]])
									key_inv = "/".join(l)
									v = frq[key_orig]
									v_inv = dict()
									for q in range(1,4):
										q1 = sorted([tmp_dict[l[0]],tmp_dict[l[q]]])
										stmp = list(s-{q})
										q2 = sorted([tmp_dict[l[stmp[0]]],tmp_dict[l[stmp[1]]]])
										if q1[0]<q2[0]:
											v_inv[l[q]] = v[q1[1]]	
										else:
											v_inv[l[q]] = v[q2[1]]
									if key_inv in TotalKey:
										vt = TotalKey[key_inv] 
										for keyt in vt.keys():
											vt[keyt].append(v_inv[keyt])
									else:
										vt = dict()
										for q in v_inv:
											vt[q] = list()
											vt[q].append(v_inv[q])
									TotalKey[key_inv] = vt
									
	TotalKeyf = dict()
	for q,v in TotalKey.iteritems():
		vtt = dict()
		for q2,v2 in v.iteritems():
			if method == "gmean":
				vtt[q2] = stats.gmean(v2)
			elif method == "mean":
				vtt[q2] = mean(v2)
			else:
				vtt[q2] = sqrt(mean(square(v2)))
		TotalKeyf[q] = vtt
	return TotalKeyf
def listQuartetsAroundBipartitions(sp):
    leaves = sp.leaf_nodes()
    setLeaves = set(leaves)
    newRoot = random.sample(leaves)
    sp.reroot_at_edge(newRoot.edge,update_bipartitions=True)
    for nd in sp.postorder_node_iter():
        if nd == root:
            break
        if nd.parent == root
            siblingNodes = nd.sibling_nodes()
            for sibl in siblingNodes:
                if size(sibl.leaf_nodes())>
    
    return