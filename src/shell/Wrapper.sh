#!/bin/bash
s=fm
e=bionj
c=8
d=prod
o=1
mean='awk "BEGIN{s=0;n=0}{n+=1;s+=\$1}END{print s/n}"'
true_sp=../../../data/11-taxon/10_taxon/10-taxon/higher-ILS/true-speciestrees
gt=../../../data/11-taxon/10_taxon/10-taxon/higher-ILS/true-genetrees


for i in {1..20}; do
./firstTest.sh -s $s -e $e -c  $true_sp/R"$i".true.tre -d $d -o $o $gt/R"$i"/true.nwk
done
#cat $gt/R*/true.nwk_d_"$e".t_res_stat.txt> $gt/res_stat.txt
cat $gt/R*/true.nwk_d_fastme_tree.nwk_res_stat.txt> $gt/res_stat.txt
rm $gt/R*/true.nwk_*
rm $gt/R*/log*
cat $gt/res_stat.txt | awk '{print $3}' | mean
