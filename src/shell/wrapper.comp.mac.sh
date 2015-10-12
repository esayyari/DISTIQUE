#!/bin/bash
DIR=$( dirname ${BASH_SOURCE[0]} )
set -e
wrkroot=~/Desktop/data/mammalian
sp=~/Documents/Research/data/mammalian/mammalian-model-species.tre
resroot=~/Desktop/results/mammalian/njst
method=njst
dataset=mammalian
for x in `find $wrkroot -maxdepth 1 -type d -name "*X*"`; do
	for i in {1..20}; do
		x=$(echo $x | sed -e 's/^.*mammalian\///')
		ls $wrkroot/$x/R$i
		
		if [ ! -d $resroot/$x ]; then
			mkdir -p $resroot/$x
			printf "$resroot/$x has been created"
		else
			printf "$resroot/$x already exists"
		fi
		if [ -s $wrkroot/$x/R$i/njst/distance.d_njst.nwk ]; then
		$DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/distance.d_njst.nwk > $resroot/$x/$dataset-$x-$method-R$i.txt;
		printf "working on $resroot/$x/$dataset-$x-$method-R$i.txt has been finished.\n"
		else
			echo	$wrkroot/$x/R$i/distance.d_fastme_tree.nwk was not found!
		fi
	done

done
