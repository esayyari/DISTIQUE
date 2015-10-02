#!/bin/bash

DIR=$( dirname ${BASH_SOURCE[0]} && pwd )

data_dir=~/Desktop/data
Astral=$WS_HOME/Astral

for files in `find $data_dir -name 'genetrees.gt'`; do
	out=$( dirname ${files})
	mkdir -p $out/astral
	gtime -po $out/astral/time.info  java -Xmx5000M -jar  $Astral/astral.4.7.8.jar -i $files  -o $out/astral/distance.d_astral.nwk 
	
	printf "working on $file has been finished\n"
	
done
