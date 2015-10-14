#!/usr/bin/env python
import dendropy
def findTrueAverageTable(frq,list_taxa):
	n = len(list_taxa)
	print list_taxa
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
											vt[keyt] += v_inv[keyt]
									else:
										TotalKey[key_inv] = v_inv
	return TotalKey
										
									

