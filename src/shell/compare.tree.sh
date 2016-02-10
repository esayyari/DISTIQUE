#!/bin/bash
DIR=$( dirname "${BASH_SOURCE[0]}" )
source $DIR/setup.sh
source $DIR/setupData.sh
show_help(){
cat << EOF
USAGE: ${0##*/} [-h] [-s SPECIES TREE] [-g ESTIMATED SPECIES TREE ]
Compares two trees. SPECIES TREE is the species tree or the reference tree, while ESTIMATED SPECIES TREE is the tree infered from gene trees.
EOF
}

while getopts "hs:g:i:" opt; do
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
	i)
		i=$OUTARG
		;;
        '?')
                printf "Unknown input option"
                show_help
                ;;
        esac
done

if [ -z "$s" ]; then
        printf "Enter a species tree"
	show_help
        exit 1
fi

if [ -z "$g" ]; then
        printf "Enter an infered tree"
	show_help
        exit 1
fi
res_root=$( dirname "${g}")
tmp1=`mktemp`
tmp2=`mktemp`
head -n 1 $s | sed -e 's/\[&U\] \|\[&R\] //'>$tmp1
head -n 1 $g | sed -e 's/\[&U\] \|\[&R\] //'>$tmp2
#for x in `head -n 1 $s`; do
# echo -n "$x" > $tmp1 #| sed -e 's/:-[0-9]*\.[0-9]*\|:[0-9]*\.[0-9]*//g;s/;//g'>$tmp1
#done
#for x in `head -n 1 $g`; do
# echo -n "$x" | sed -e 's/:-[0-9]*\.[0-9]*\|:[0-9]*\.[0-9]*//g;s/;//g'>$tmp2
#done
$WS_GLB_SH/compareTrees.missingBranch $tmp1 $tmp2 -simplify

rm $tmp1
rm $tmp2
