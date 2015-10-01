#!/bin/bash
DIR=$( dirname ${BASH_SOURCE[0]} )
set -e
wrkroot=~/data/outputs/new/avian
source $DIR/setup.sh
sp=~/data/Avian/avian-model-species.tre
resroot=~/data/results/Newresults/avian
method=prod_fm
dataset=avian
for x in `ls $wrkroot`; do
	for i in {1..20}; do
		if [ ! -d $resroot/$x ]; then
			mkdir -p $resroot/$x
		fi
		if [ -s $wrkroot/$x/R$i/distance.d_fastme_tree.nwk ]; then
		$DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/distance.d_fastme_tree.nwk > $resroot/$x/$dataset-$x-$method-R$i.txt;
		printf "working on $resroot/$x/$dataset-$x-$method-R$i.txt has been finished.\n"
		fi
	done

done
