#!/bin/bash

for i in {1..20}; do
./firstTest.sh -m pd -d mvr -c  ../../../data/11-taxon/10_taxon/10-taxon/higher-ILS/true-speciestrees/R"$i".true.tre ../../../data/11-taxon/10_taxon/10-taxon/higher-ILS/true-genetrees/R"$i"/true.nwk
done

cat ../../../data/11-taxon/10_taxon/10-taxon/higher-ILS/true-genetrees/R*/true.nwk_d_mrv.t_res_stat.txt> ../../../data/11-taxon/10_taxon/10-taxon/higher-ILS/true-genetrees/res_stat.txt
