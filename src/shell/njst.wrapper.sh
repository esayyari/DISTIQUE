#!/bin/bash

DIR=$( dirname ${BASH_SOURCE[0]} && pwd )
data_out=$WS_HOME/results/tests/11-taxon/njst
data_dir=$WS_HOME/data/11-taxon
#data_dir=$WS_HOME/data//11-taxon/model.10.5400000.0.000000037/11/
NJstS=$WS_HOME/ASTRID/src

for files in `find $data_dir -name 'all_fasttree.genetrees'`; do
	out=$( dirname ${files})
	mkdir -p $out/njst
	if [ -s "$out/njst/distance.csv" ]; then
		tmp_x=$(echo $out | sed -e 's/^.*11-taxon\/*//')
		echo $tmp_x
		mkdir -p $data_out/$tmp_x
		mv $out/njst/* $data_out/$tmp_x/
	else
	python $NJstS/ASTRID.py	-i $files -m fastme2 -o $out/njst/distance.d_njst.nwk -c $out/njst/distance.csv
		tmp_x=$(echo $out | sed -e 's/^.*11-taxon\/*//')
		echo $tmp_x
		mkdir -p $data_out/$tmp_x
                mv $out/njst/* $data_out/$tmp_x/
	fi
	printf "working on $files has been finished\n"
#	if [ -s "$data_out/$tmp_x/njst/distance.csv" ]; then
#		printf "CACHE: $data_out/$tmp_x/njst/distance.csv\n"
#	else
#		printf "error: CACHE does not found"
#		exit 1
#	fi
	
done
