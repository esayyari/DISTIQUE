#!/bin/bash

set -x
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/setup.sh
source $DIR/setupData.sh

show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-g GENE TREE] [-n NUMBER OF ANCHORES] [ -o OUTPUT DIR ] [ -r REPEATS] [ -u SUMMET][-z SUMMETOPT][-v VERSION]
Estimate species tree using DISTIQUE  method. GNEE TREES FILE is the file of gene trees concatenated.
	-g 	GENE TREE 		The gene tree file for which we want to infer species tree
	-n 	NUMBER OF ANCHORES	The number of anchors 
	-o  	OUTPUT DIR		The output directory, to copy the results
	-r 	REPEATS			The number of repeating the experiments
EOF
}
n=10
m=log
l=const
if [ $# -lt 1 ]; then
	show_help
	exit 0;
fi  
rt=1
while getopts "hg:n:o:r:m:l:u:z:v:" opt; do
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
	l)
		l=$OPTARG
		;;
	u) 	
		u=$OPTARG
		;;
	z)
		z=$OPTARG
		;;
	v)
		v=$OPTARG
		;;
	esac
done
out_final=`mktemp -d`
rt=1
END=$rt
START=1
for (( c=$START; c<=$END; c++ ))
do
out=`mktemp -d`
if [ "$v" == "5" or "$v" == "3" ]; then
$WS_LOC_PUTIL/testAnchoring-v$v.py -g $gt -o $out -n $nt -u $u -z $z

else
$WS_LOC_PUTIL/testAnchoring-v$v.py -g $gt -o $out -n $nt -u $u -z $z
$WS_LOC_SH/mrl.sh $out
cat $out/distance-*.nwk* > $out/anchored_trees
cat $out/anchored_mrl_tree.nwk >> $out/anchored_trees
fi
tar czvf  $out/results$v-$nt-$z-$u.tar.gz $out/*
mv $out/results$v-$nt-$z-$u.tar.gz $out_final
done
tar czvf $o/results$v-$nt-$z-$u.tar.gz $out_final/*
#mv $out_final/results$nt-$m.tar.gz $o
rm -r $out

