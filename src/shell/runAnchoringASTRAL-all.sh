#!/bin/bash
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
set -x

show_help(){
cat << EOF
USAGE: ${0##*/} [-h] [-i DATA PATH] [-o OUTPATH] [-n NUMROUNDS] [-v VERSION]
Runs DISTIQUE Anchoring on given gene tree available at DATA PATH. The results will be written at OUTPATH. NUMROUNDS is number of rounds for sampling around each polytomy, while VERSION indicates which version of the code to use.
EOF
}
n=1000
z=B
while getopts "hg:f:o:n:z:" opt; do
        case $opt in
        h)
                show_help
                exit 0
                ;;
        g)
                g=$OPTARG
                ;;
	f)
		f=$OPTARG
		;;
        o)
                o=$OPTARG
                ;;
	n)
		n=$OPTARG
		;;
	
	z) 	
		z=$OPTARG
		;;
        '?')
                printf "Unknown input option"
                show_help
                ;;
        esac
done

dir=$(dirname $g)
WS_LOC_UTIL=$WS_HOME/DISTIQUE/src/utils/
tmpDIR=`mktemp -d`
echo $tmpDIR

gt1000half=`mktemp  ${tmpDIR}/genetrees"$n".half.XXXXX` || exit 1

head -n $n $g > $gt1000half
sp=$dir/s_tree.trees


res1000half=`mktemp -d $tmpDIR/"$n"genes.half.XXXXX` || exit 1

if [ -s $gt1000half ]; then
/usr/bin/time -p $WS_LOC_UTIL/distique.py -g $gt1000half -o $res1000half -f $f -u fastme -z $z > $res1000half/results.log 2>&1
echo "working on $dir/estimatedgenetre.halfresolved$n has been finished!"
fi

for x in `find $tmpDIR -name "distance.*"`; do
	y=$(dirname $x)
	$WS_HOME/DISTIQUE/src/shell/compare.tree.sh -s $sp -g $x > $y/results-distique-all.score
done

y=$(echo $o | sed -e 's/^.*\///')
tar czf $o/testDISTIQUE-all-$n-$z-fastme.tar.gz $tmpDIR/* 
rm -r $tmpDIR
