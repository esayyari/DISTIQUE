#!/bin/bash
DIR=$( dirname "${BASH_SOURCE[0]}" )
show_help(){
cat << EOF
USAGE: ${0##*/} [-h] [-s SPECIES TREE] [-g GENE TREE ] 
Compares two trees. SPECIES TREE is the species tree or the reference tree, while GENE TREE is the tree infered from gene trees.
EOF
}


while getopts "hs:g:" opt; do
	case $opt in
	h)
		show_help
		exit 0
		;;
	s)
		s=$OPTARG
		;;
	g)
		g=$OPTARG
		;;
	'?')	
		printf "Unknown input option"
		show_help
		;;
	esac
done

if [ -z "$s" ]; then
	printf "Enter a species tree"
	exit 1
fi

if [ -z "$g" ]; then
	printf "Enter an infered tree"
	exit 1
fi
res_root=$( dirname "${g}")
tmp1=`mktemp`
tmp2=`mktemp`
head -n 1 "$s">$tmp1
head -n 1 "$g">$tmp2
$WS_GLB_SH/compareTrees.missingBranch $tmp1 $tmp2
			
