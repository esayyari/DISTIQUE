#!/bin/bash
DIR=$( pwd )
set -e
wrkroot=/oasis/projects/nsf/uot136/esayyari/outputs
source $DIR/setup.sh
dataset=11-taxon
sp=/oasis/projects/nsf/uot136/esayyari/data/$dataset/$dataset-model-species.tre
resroot=./oasis/projects/nsf/uot136/esayyari/results/$dataset
method=prod_fm
m=min
for x in `find $wrkroot/$dataset -maxdepth 1 -type d -name "model*"`; do
	for i in {1..50}; do
		x=$( echo $x | sed -e 's/^.*outputs\///' )
		echo $x
		tmp_x=$(echo $x | sed -e 's/\//-/g' );
		echo $tmp_x	
		sp=/oasis/projects/nsf/uot136/esayyari/data/$x/$i/S_relabeled.s
		
		for y in `ls $wrkroot/$x/$i`; do
			y=$( echo $y | sed -e 's/^.*\/$x\/[0-9]*\///')	
			if [ -s $wrkroot/$x/R$i/distance.d_fastme_tree.nwk ]; then
				printf "mkdir -p $resroot/$x; module load python; module load scipy; $DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/distance.d_fastme_tree.nwk > $resroot/$x/$dataset-prod_fm-fastme-$tmp_x-R$i.txt\n" >>list.comp.txt;
				printf "mkdir -p $resroot/$x; working on $resroot/$x/$dataset-prod_fm-fastme-$tmp_x-$method-R$i.txt has been finished.\n"
			else
				printf "mkdir -p $resroot/$x; $wrkroot/$x/R$i/distance.d_fastme_tree.nwk was not found \n"
			       
			fi
			if [ -s $wrkroot/$x/R$i/njst/distance.d_njst_tree.nwk ]; then
				printf "mkdir -p $resroot/$x; module load python; module load scipy; $DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/njst/distance.d_njst_tree.nwk > $resroot/$x/$dataset-astrid-fastme-$tmp_x-R$i.txt\n" >> list.comp.txt
				printf "mkdir -p $resroot/$x; working on $resroot/$x/$dataset-astrid-fastme-$tmp_x-$method-R$i.txt has been finished.\n"
			else
				printf "mkdir -p $resroot/$x; $wrkroot/$x/R$i/njst/distance.d_njst_tree.nwk was not found \n"
			fi
			if [ -s $wrkroot/$x/R$i/astral/distance.d_astral_tree.nwk ]; then
				printf "mkdir -p $resroot/$x; module load python; module load scipy; $DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/astral/distance.d_astral_tree.nwk > $resroot/$x/$dataset-astral-$tmp_x-R$i.txt\n" >> list.comp.txt
				printf "mkdir -p $resroot/$x; working on $resroot/$x/$dataset-astral-$tmp_x-$method-R$i.txt has been finished.\n"
			else
				printf "mkdir -p $resroot/$x; $wrkroot/$x/R$i/astral/distance.d_astral_tree.nwk was not found \n"
		       fi
			if [ -s $wrkroot/$x/R$i/distique-cons/min/min/distance.d_distique_tree.nwk ]; then
				printf "mkdir -p $resroot/$x; $DIR/compare.tree.sh -s $sp -g $wrkroot/$x/R$i/distique-cons/min/min/distance.d_distique_tree.nwk > $resroot/$x/$dataset-distique-cons-min-$tmp_x-R$i.txt\n">> list.comp.txt
			else
				printf "mkdir -p $resroot/$x; wrkroot/$x/R$i/distique-cons/min/min/distance.d_distique_tree.nwk was not found \n"
			fi
			if [ -s $wrkroot/$x/R$i/distique-cons/prod/prod/distance.d_distique_tree.nwk ]; then
				printf "mkdir -p $resroot/$x; $DIR/compare.tree.sh -s $sp -g $wrkroot/$x/R$i/distique-cons/prod/prod/distance.d_distique_tree.nwk > $resroot/$x/$dataset-distique-cons-prod-$tmp_x-R$i.txt\n">> list.comp.txt
			else
				printf "mkdir -p $resroot/$x; $wrkroot/$x/R$i/distique-cons/prod/prod/distance.d_distique_tree.nwk was not found \n"
			fi
			if [ -s $wrkroot/$x/R$i/distique-2/prod/prod/distance.d_distique_tree.nwk ]; then
				printf "mkdir -p $resroot/$x; $DIR/compare.tree.sh -s $sp -g $wrkroot/$x/R$i/distique-2/prod/prod/distance.d_distique_tree.nwk > $resroot/$x/$dataset-distique-2-prod-$tmp_x-R$i.txt\n">> list.comp.txt
			else
				printf "mkdir -p $resroot/$x; $wrkroot/$x/R$i/distique-2/prod/prod/distance.d_distique_tree.nwk was not found\n"
			fi
			if [ -s $wrkroot/$x/R$i/distique-2/min/min/distance.d_distique_tree.nwk ]; then
				printf "mkdir -p $resroot/$x; $DIR/compare.tree.sh -s $sp -g $wrkroot/$x/R$i/distique-2/min/min/distance.d_distique_tree.nwk > $resroot/$x/$dataset-distique-2-min-$tmp_x-R$i.txt\n">> list.comp.txt
			else
				printf "mkdir -p $resroot/$x; $wrkroot/$x/R$i/distique-2/min/min/distance.d_distique_tree.nwk was not found"
			fi

		#	printf "mkdir -p $resroot/$x; module load python; module load scipy; $DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/distance.d_fastme_tree.nwk > $resroot/$x/$dataset-prod_fm-fastme-$tmp_x-R$i.txt\n" >>list.comp.txt;
		#	printf "mkdir -p $resroot/$x; working on $resroot/$x/$dataset-prod_fm-fastme-$tmp_x-$method-R$i.txt has been finished.\n"
                 #       printf "mkdir -p $resroot/$x; $wrkroot/$x/R$i/distance.d_fastme_tree.nwk was not found \n"
        	done               
	done

done
