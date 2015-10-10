#!/bin/bash
DIR=$( dirname "${BASH_SOURCE[0]}" )
data_dir=~/Documents/Research/results/Astrid/11-taxon/njst
out=~/Documents/Research/results/Astrid/11-taxon/njst/results
#if [ ! -d "$data_dir/unzip" ]; then	
#mkdir -p $data_dir/unzip
#fi
for data_name in `ls $data_dir`; do
#	unzip $data_dir/$data_name -d $data_dir/unzip/
	
       	out_file=$(echo $data_name | sed -e 's/out\.fm\.\(bionj\|nj\|mvr\|unj\)\.\(min\|prod\)\.8\.1//g;s/\(\.zip\)$//g;s/-\(lowILS\|highILS\|medILS\|veryHighILS\)$//g')
	fm_res=distance.d_njst.nwk
	path_to_true_sp_tree=~/Documents/Research/data/11-taxon/$out_file/
	if [ -s $out/$data_name-res.txt ]; then
	rm $out/$data_name-res.txt
	fi
 	if [ ! -d $out/$data_name-dir ]; then
	mkdir -p $out/$data_name-dir
	fi
#	for file in `find $data_dir/unzip/$data_name-dir/ -name $pd_res -o -name $fm_res`; do
	for file in `find $data_dir/$data_name  -name $fm_res`; do
		path=$( dirname ${file} )
		tr_sp_tmp=$( echo $path | sed -e "s/.*\($data_name-dir\/\)//" | sed -e "s/.*$out_file//g" )
		save_name=$( echo $path | tr '/' '.' | sed -e 's/^\.\|\(home\.esayyari\.data\.outputs\.unzip\.out\)\.//g')
		echo $save_name	
		tr_sp_path=$path_to_true_sp_tree/$tr_sp_tmp/
		tr_sp_path2=$( echo $tr_sp_path | sed -e "s/\/relabeled.*//" )
		
		echo $tr_sp_path2	
		tr_spt=$tr_sp_path2
		R=$( echo $tr_sp_path2 | sed -e 's/.*\/\///')
		echo $R
		true_tree=$( find $tr_sp_path2 -name "S_relabeled_tree.trees")
		printf "Working on $true_tree and $file has been finished\n"
				
		$DIR/compare.tree.sh -s $true_tree -g $file >> $out/$data_name-dir/$save_name-$R"."txt
#		$DIR/compare.tree.sh -s $true_tree -g $file >> $out/$data_name-dir/$save_name.txt
	done
done
