#!/bin/bash

DIR=~/repository/DISTIQUE/src/shell/
res=~/data/results
path=~/data/11-taxon/
input_path=~/data/outputs/unzip
 #11-taxon_model.10.5400000.0.000000037,11,astral_fasttree_genetrees,10,100,9,.5
for data_dir in `ls $input_path`; do
	dataset=$( echo $data_dir | grep -o "\(10-taxon\|11-taxon\|\)" )
	mkdir -p $res/$data_dir
	if [ "$dataset" == "10-taxon" ]; then
#	out.fm.bionj.min.8.110-taxon.zip-dir/10-taxon/higher-ILS/estimated-genetrees/R2/distance.d_fastme_tree.nwk	
#	true-speciestrees/R1.true.tre	
		allfiles=$( find $input_path/$data_dir -name "distance\.d_*.nwk" -o -name "distance\.d_*.t") 
	
		for file in $allfiles; do
			R=$( echo $file | grep -o "R[0-9]*" )
		 	subpath=$(  echo $file | sed -e "s/^.*\.zip-dir\///g" | sed -e "s/distance.*\///g" | sed -e "s/\/R.*//g" | sed -e "s/estimated-ge.*\|true-ge.*//g")
			mod=$( echo $file | grep -o "fm\|pd" )
			mod2=$( echo $file | grep -o "min\|prod")
			gene=$( echo $file | grep -o "lower.*\|higher.*" | sed -e "s/\/distance.d.*//g" | sed -e "s/\//\./g")
			
			sp=~/data/$subpath/true-speciestrees/$R.true.tre
			$DIR/compare.tree.sh -s $sp -g $file > $res/$data_dir/$dataset-$dataset-$mod2"_"$mod-$gene.txt
			printf "$sp\n $file\n"
		done	
	elif [ "$dataset" == "11-taxon" ]; then
		allfiles=$( find $input_path/$data_dir -name "distance\.d_*.nwk" -o -name "distance\.d_*.t") 
		for file in $allfiles; do
			
		 	subpath=$(  echo $file | sed -e "s/^.*\.zip-dir\///g" | sed -e "s/distance.*\///g" | sed -e "s/\/relabeled.*$//g" )
			model=$( echo $subpath | sed -e "s/\/.*//g" )
			mod=$( echo $file | grep -o "fm\|pd" )
			mod2=$( echo $file | grep -o "min\|prod")
			gene=$( echo $file | grep -o "[0-9]*\/relabeled_.*"| sed -e "s/\/distance.d.*//g" | sed -e "s/\//\./g")
			R=$( echo $gene | sed -e 's/\..*//' )
			sp=~/data/$dataset/$subpath/S_relabeled_tree.trees
			$DIR/compare.tree.sh -s $sp -g $file> $res/$data_dir/$dataset"_"$model-$dataset-$mod2"_"$mod-$gene-$R.txt
			printf "$sp\n $file\n"
			echo "$res/$dataset"_"$model-$dataset-$mod2"_"$mod-$gene-$R.txt"
		done	

	fi
done


