#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

source $DIR/setup.sh

test $# == 4 || { echo USAGE: dirname-genelist dirname-BestMLTrees rep-count output-dir; exit 1; }
genelist=$1
MLlist=$2
rep=$3
out=$4
best=RAxML_bipartitions.final.f200
bs=RAxML_bootstrap.all
for i in {0..2022}; do
ls $genelist/*/bin.$i.txt/*/$bs >> $out/genelist
cat $MLlist/*/bin.$i.txt/*/$best >> $out/Best
done
java -jar $WS_HOME/Astral/astral.4.7.8.jar -k bootstraps_norun -i $out/Best -b $out/genelist -o $out/BS -r $rep
