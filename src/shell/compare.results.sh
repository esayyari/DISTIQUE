#!/bin/bash
DIR=$( dirname "${BASH_SOURCE[0]}" )
data_dir=~/data/outputs
out=~/data/results
if [ ! -d "$data_dir/unzip" ]; then	
mkdir -p $data_dir/unzip
fi
for data_name in `ls $data_dir`; do
#	unzip $data_dir/$data_name -d $data_dir/unzip/
	
       	out_file=$(echo $data_name | sed -e 's/out\.\(fm\|pd\)\.\(bionj\|nj\|mvr\|unj\)\.\(min\|prod\)\.8\.1//g;s/\(\.zip\)$//g;s/-\(lowILS\|highILS\|medILS\|veryHighILS\)$//g')
	sum_method=$( echo $data_name | grep -o '\(fm\|pd\)')
	if [ "$sum_method" == "pd" ]; then
	sum_opt=$( echo $data_name | grep -o '\(bionj\|nj\|mvr\|unj\)' )
	pd_res="distance.d_$sum_opt.t"
	fi
	fm_res="distance.d_fastme_tree.nwk"
	path_to_true_sp_tree=~/data/$out_file/
	if [ -s $out/$data_name-res.txt ]; then
	rm $out/$data_name-res.txt
	fi
 	if [ ! -d $out/$data_name-dir ]; then
	mkdir -p $out/$data_name-dir
	fi
	for file in `find $data_dir/unzip/$data_name-dir/ -name $pd_res -o -name $fm_res`; do
		path=$( dirname ${file} )
		tr_sp_tmp=$( echo $path | sed -e "s/.*\($data_name-dir\/\)//g" | sed -e "s/.*$out_file//g" )
		save_name=$( echo $path | tr '/' '.' | sed -e 's/^\.\|\(home\.esayyari\.data\.outputs\.unzip\.out\)\.//g')
		echo $save_name	
		tr_sp_path=$path_to_true_sp_tree/$tr_sp_tmp/
		
		if [ "$out_file" ==  "11-taxon" ]; then
		tr_sp_path2=$( echo $tr_sp_path | sed -e "s/relabeled.*//g" )
		
		tr_spt=$tr_sp_path2
		true_tree=$( find $tr_sp_path2 -name "S_relabeled_tree.trees")
		printf "Working on $true_tree and $file has been finished\n"
				
		$DIR/compare.tree.sh -s $file -g $file >> $out/$data_name-dir/$save_name-chk.txt
#		$DIR/compare.tree.sh -s $true_tree -g $file >> $out/$data_name-dir/$save_name.txt
		elif [ "$out_file" == "10-taxon" ]; then
		tr_sp_tmp1=$( echo $tr_sp_path | sed -e 's/\(true-genetrees\|estimated-genetrees\).*$//g' )
		tmp_name=$( echo $tr_sp_path | grep -o '\(true-genetrees\|estimated-genetrees\)')
		R=$(echo $path | grep -o '\(R[0-9]*\)' )
		tr_sp=$tr_sp_tmp1/true-speciestrees/$R.true.tre
#		$DIR/compare.tree.sh -s $tr_sp -g $file >> $out/$data_name-dir/$save_name.txt
		$DIR/compare.tree.sh -s $file -g $file >> $out/$data_name-dir/$save_name-chk.txt
		printf "Working on $tr_sp and $file has been finished\n"
		fi
	done
done
