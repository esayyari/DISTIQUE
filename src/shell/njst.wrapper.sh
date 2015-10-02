#!/bin/bash

DIR=$( dirname ${BASH_SOURCE[0]} && pwd )

data_dir=$WS_HOME/data/11-taxon
#data_dir=$WS_HOME/data//11-taxon/model.10.5400000.0.000000037/11/
NJstS=$WS_HOME/ASTRID/src

for files in `find $data_dir -name 'all_fasttree.genetrees'`; do
	out=$( dirname ${files})
	mkdir -p $out/njst
	python $NJstS/ASTRID.py	-i $files -m fastme2 -o $out/njst/distance.d_njst.nwk -c CACHE
	
	mv $DIR/CACHE* $out/njst/
	cp $out/njst/CACHE $out/njst/distance.csv
	printf "working on $file has been finished\n"
	if [ -s "$out/njst/distance.csv" ]; then
		printf "CACHE: $out/njst/distance.csv\n"
	else
		printf "error: CACHE does not found"
		exit 1
	fi
	
done
