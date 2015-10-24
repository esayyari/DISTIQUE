#!/bin/bash

DIR=~/repository/DISTIQUE/src/shell/
res=/oasis/projects/nsf/uot136/esayyari/results
path=/oasis/projects/nsf/uot136/esayyari/data/11-taxon/
input_path=/oasis/projects/nsf/uot136/esayyari/outputs/11-taxon/
dataset=11-taxon
for data_dir in `ls $input_path`; do
	mkdir -p $res/$data_dir
	echo $data_dir
	for i in {1..50}; do
		for y in `ls $input_path/$data_dir/$i/`;do
#			if [ -s $input_path/$data_dir/$i/$y/njst/distance.d_njst_tree.nwk ]; then
#				printf "mkdir -p ./$input_path/$data_dir/$i/$y/njst; $DIR/compare.tree.sh -g $input_path/$data_dir/$i/$y/njst/distance.d_njst_tree.nwk -s $path/$data_dir/$i/S_relabeled_tree.trees > ./$input_path/$data_dir/$i/$y/njst/11-taxon-njst-11-taxon-$data_dir-R$i.txt\n">>tasks.massive
#			else
#				printf "the file $input_path/$data_dir/$i/$y/njst/distance.d_njst_tree.nwk was not found\n"
#			fi
#			if [ -s $input_path/$data_dir/$i/$y/astral/distance.d_astral_tree.nwk ]; then
#				printf "mkdir -p ./$input_path/$data_dir/$i/$y/astral; $DIR/compare.tree.sh -g $input_path/$data_dir/$i/$y/astral/distance.d_astral_tree.nwk -s $path/$data_dir/$i/S_relabeled_tree.trees > ./$input_path/$data_dir/$i/$y/astral/11-taxon-astral-11-taxon-$data_dir-R$i.txt\n">>tasks.massive
#			else
#				printf "the file $input_path/$data_dir/$i/$y/astral/distance.d_astral_tree.nwk was not found\n"
#			fi
#			if [ -s $input_path/$data_dir/$i/$y/distance.d_fastme_tree.nwk ]; then
#			printf "mkdir -p ./$input_path/$data_dir/$i/$y/distique-1; $DIR/compare.tree.sh -g $input_path/$data_dir/$i/$y/distance.d_fastme_tree.nwk -s $path/$data_dir/$i/S_relabeled_tree.trees > ./$input_path/$data_dir/$i/$y/distique-1/11-taxon-distique-1-prod-11-taxon-$data_dir-R$i.txt\n">>tasks.massive
 #                       else
  #                              printf "the file $input_path/$data_dir/$i/$y/distance.d_fastme_tree.nwk was not found\n"
#			fi
			#if [ -s $input_path/$data_dir/$i/$y/distique-2/min/min/distance.d_distique_tree.nwk ]; then
			#printf "mkdir -p ./$input_path/$data_dir/$i/$y/distique-2/min/min; $DIR/compare.tree.sh -g $input_path/$data_dir/$i/$y/distique-2/min/min/distance.d_distique_tree.nwk -s $path/$data_dir/$i/S_relabeled_tree.trees > ./$input_path/$data_dir/$i/$y/distique-2/min/min/11-taxon-distique-2-min-11-taxon-$data_dir-R$i.txt\n">>tasks.massive
                        #else
                         #       printf "the file $input_path/$data_dir/$i/$y/distique-2/min/min/distance.d_distique_tree.nwk was not found\n"
			#fi
			#if [ -s $input_path/$data_dir/$i/$y/distique-2/prod/prod/distance.d_distique_tree.nwk ]; then
			#printf "mkdir -p ./$input_path/$data_dir/$i/$y/distique-2/prod/prod/; $DIR/compare.tree.sh -g $input_path/$data_dir/$i/$y/distique-2/prod/prod/distance.d_distique_tree.nwk -s $path/$data_dir/$i/S_relabeled_tree.trees > ./$input_path/$data_dir/$i/$y/distique-2/prod/prod/11-taxon-distique-2-prod-11-taxon-$data_dir-R$i.txt\n">>tasks.massive
                        #else
                         #       printf "the file $input_path/$data_dir/$i/$y/distique-2/prod/prod/distance.d_distique_tree.nwk was not found\n"
			#fi
			#if [ -s $input_path/$data_dir/$i/$y/distique-cons/prod/prod/distance.d_distique_tree.nwk ]; then
			#	printf "mkdir -p ./$input_path/$data_dir/$i/$y/distique-cons/prod/prod/; $DIR/compare.tree.sh -g $input_path/$data_dir/$i/$y/distique-cons/prod/prod/distance.d_distique_tree.nwk -s $path/$data_dir/$i/S_relabeled_tree.trees > ./$input_path/$data_dir/$i/$y/distique-cons/prod/prod/11-taxon-distique-cons-prod-11-taxon-$data_dir-R$i.txt\n">>tasks.massive
                        #else
                         #       printf "the file $input_path/$data_dir/$i/$y/distique-cons/prod/prod/distance.d_distique_tree.nwk was not found\n"
			#fi
			#if [ -s $input_path/$data_dir/$i/$y/distique-cons/min/min/distance.d_distique_tree.nwk ]; then
			#	printf "mkdir -p ./$input_path/$data_dir/$i/$y/distique-cons/min/min/; $DIR/compare.tree.sh -g $input_path/$data_dir/$i/$y/distique-cons/min/min/distance.d_distique_tree.nwk -s $path/$data_dir/$i/S_relabeled_tree.trees > ./$input_path/$data_dir/$i/$y/distique-cons/min/min/11-taxon-distique-cons-prod-11-taxon-$data_dir-R$i.txt\n">>tasks.massive
                        #else
                         #       printf "the file $input_path/$data_dir/$i/$y/distique-cons/min/min/distance.d_distique_tree.nwk was not found\n"
			#fi
			if [ -s $input_path/$data_dir/$i/$y/cons/distance.d_cons_tree.nwk ]; then
				printf "mkdir -p ./$input_path/$data_dir/$i/$y/cons; $DIR/compare.tree.sh -g $input_path/$data_dir/$i/$y/cons/distance.d_cons_tree.nwk -s $path/$data_dir/$i/S_relabeled_tree.trees > .//$input_path/$data_dir/$i/$y/cons/11-taxon-cons-11-taxon-$data_dir-$y-R$i.txt\n">>list.comp.txt
                        else
                                printf "the file $input_path/$data_dir/$i/$y/distique-cons/min/min/distance.d_distique_tree.nwk was not found\n"
			fi
		done
	done		 

done

