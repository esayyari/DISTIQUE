#!/bin/bash

#set -x
#set -e
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/setup.sh
source $DIR/setupData.sh

show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-g GENE TREE] [-n NUMBER OF ANCHORES] [ -o OUTPUT DIR ] 
Estimate species tree using DISTIQUE  method. GNEE TREES FILE is the file of gene trees concatenated.
	-g 	GENE TREE 		The gene tree file for which we want to infer species tree
	-n 	NUMBER OF ANCHORES	The number of anchors 
	-o  	OUTPUT DIR		The output directory, to copy the results
EOF
}
n=10
if [ $# -lt 1 ]; then
	show_help
	exit 0;
fi  
while getopts "hg:n:o:" opt; do
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
	esac
done

out=`mktemp -d`
$WS_LOC_PUTIL/testAnchoring.py -g $gt -o $out -n $nt 
$WS_LOC_SH/mrl.sh $out
cat $out/distance-*.nwk > $out/anchored_trees
cat $out/anchored_mrl_tree.nwk >> $out/anchored_trees
tar czvf  $out/results.tar.gz $out/*
mv $out/results.tar.gz $o
rm -r $out
