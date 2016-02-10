#!/bin/bash

#set -x
#set -e
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/setup.sh
source $DIR/setupData.sh

show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-g GENE TREE] [-n NUMBER OF ANCHORES] [ -o OUTPUT DIR ] [ -r REPEATS]
Estimate species tree using DISTIQUE  method. GNEE TREES FILE is the file of gene trees concatenated.
	-g 	GENE TREE 		The gene tree file for which we want to infer species tree
	-n 	NUMBER OF ANCHORES	The number of anchors 
	-o  	OUTPUT DIR		The output directory, to copy the results
	-r 	REPEATS			The number of repeating the experiments
EOF
}
n=10
m=log
if [ $# -lt 1 ]; then
	show_help
	exit 0;
fi  
while getopts "hg:n:o:r:m:" opt; do
	case $opt in 
	h) 
		show_help
		exit 0;
		;;
	g)
		gt=$OPTARG
		;;
	n) 	nt=$OPTARG
		;;
	o)
		o=$OPTARG
		;;
	r) 
		rt=$OPTARG
		;;
	m)
		m=$OPTARG
		;;
	esac
done
out_final=`mktemp -d`
END=$rt
START=1
for (( c=$START; c<=$END; c++ ))
do
out=`mktemp -d`
$WS_LOC_PUTIL/testAnchoring.py -g $gt -o $out -n $nt -m $m
$WS_LOC_SH/mrl.sh $out
cat $out/distance-*.nwk > $out/anchored_trees
cat $out/anchored_mrl_tree.nwk >> $out/anchored_trees
tar czvf  $out/results$nt-$c.tar.gz $out/*
mv $out/results$nt-$c.tar.gz $out_final
done
tar czvf $out_final/results$nt-$m.tar.gz $out_final/*
mv $out_final/results$nt-$m.tar.gz $o
rm -r $out

