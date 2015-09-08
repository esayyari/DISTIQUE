#!/bin/bash
source ~/.bashrc
DIR=`pwd`
source $DIR/setup.sh
s=pd
e=nj
c=8
d=prod
o=1
true_sp=../../../data/11-taxon/10_taxon/10-taxon/higher-ILS/true-speciestrees
gt=../../../data/11-taxon/10_taxon/10-taxon/higher-ILS/true-genetrees

for i in {1..20}; do
./summarize.tree.sh -w 1 -s $s -e $e  -d $d -o $o -r $WS_LOC_RES/R"$i" $gt/R"$i"/true.nwk 
./compare.tree.sh -s $true_sp/R"$i".true.tre -g $WS_LOC_RES/R"$i"/distance.d_"$e".t > $WS_LOC_RES/R"$i"/res.txt
#./compare.tree.sh -s $true_sp/R"$i".true.tre -g $WS_LOC_RES/R"$i"/distance.d_fastme_tree.nwk>$WS_LOC_RES/R"$i"/res.txt

done
cat $WS_LOC_RES/R*/res.txt > $WS_LOC_RES/res_stat.txt
