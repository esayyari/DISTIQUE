#!/bin/bash
DIR=$( pwd )
set -e
wrkroot=/oasis/projects/nsf/uot136/esayyari/outputs
source $DIR/setup.sh
dataset=mammalian
sp=/oasis/projects/nsf/uot136/esayyari/data/$dataset/$dataset-model-species.tre
resroot=/oasis/projects/nsf/uot136/esayyari/results/$dataset
method=prod_fm
for x in `find $wrkroot/$dataset -maxdepth 1 -type d -name "*X*"`; do
	for i in {1..20}; do
		x=$( echo $x | sed -e 's/^.*outputs\///' )
		tmp_x=$(echo $x | sed -e 's/\//-/g' );
		if [ ! -d $resroot/$x ]; then
			mkdir -p $resroot/$x
		fi
		if [ -s $wrkroot/$x/R$i/distance.d_fastme_tree.nwk ]; then
			printf "module load python; module load scipy; $DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/distance.d_fastme_tree.nwk > $resroot/$x/$dataset-prod_fm-fastme-$tmp_x-R$i.txt\n" >>list.comp.txt;
			printf "working on $resroot/$x/$dataset-prod_fm-fastme-$tmp_x-$method-R$i.txt has been finished.\n"
		else
                        printf "$wrkroot/$x/R$i/distance.d_fastme_tree.nwk was not found \n"
                       
		fi
                if [ -s $wrkroot/$x/R$i/njst/distance.d_njst_tree.nwk ]; then
                	printf "module load python; module load scipy; $DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/njst/distance.d_njst_tree.nwk > $resroot/$x/$dataset-astrid-fastme-$tmp_x-R$i.txt\n" >> list.comp.txt
               		printf "working on $resroot/$x/$dataset-astrid-fastme-$tmp_x-$method-R$i.txt has been finished.\n"
		else
			printf "$wrkroot/$x/R$i/njst/distance.d_njst_tree.nwk was not found \n"
                fi
		if [ -s $wrkroot/$x/R$i/astral/distance.d_astral_tree.nwk ]; then
                        printf "module load python; module load scipy; $DIR/compare.tree.sh  -s $sp -g $wrkroot/$x/R$i/astral/distance.d_astral_tree.nwk > $resroot/$x/$dataset-astral-$tmp_x-R$i.txt\n" >> list.comp.txt
                        printf "working on $resroot/$x/$dataset-astral-$tmp_x-$method-R$i.txt has been finished.\n"
                else
                        printf "$wrkroot/$x/R$i/astral/distance.d_astral_tree.nwk was not found \n"
                fi
	done

done
