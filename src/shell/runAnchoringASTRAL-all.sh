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

m=(1000 200 50)
for n in ${m[@]}; do
gt1000half=`mktemp  ${tmpDIR}/genetrees"$n".half.XXXXX` || exit 1

head -n $n $g > $gt1000half
sp=$dir/s_tree.trees
quart=$(find $dir/distique-all/distique-all/"$n".genes.half* -name "quartet_tmp.q" || exit 1)

res1000half=`mktemp -d $tmpDIR/"$n"genes.half.XXXXX` || exit 1

if [ -s $gt1000half ]; then
/usr/bin/time -p $WS_LOC_UTIL/distique.py -t 1 -m prod -a mean -l freq -g $gt1000half -f $quart -o $res1000half  > $res1000half/results.log 2>&1
echo "/usr/bin/time -p $WS_LOC_UTIL/distique.py -t 1 -m prod -a mean -l freq -g $gt1000half -f $quart -o $res1000half" >> $res1000half/results.log
DATE=`date +%Y-%m-%d:%H:%M:%S`;  info=$(cd ${WS_HOME}/DISTIQUE; git show-ref  --head --tags);
echo $DATE >> $res1000half/results.log
echo version info: >> $res1000half/results.log
echo $info >> $res1000half/results.log
echo "working on $dir/estimatedgenetre.halfresolved$n has been finished!"
fi
done
for x in `find $tmpDIR -name "distique_*"`; do
	y=$(dirname $x)
	$WS_HOME/DISTIQUE/src/shell/compare.tree.sh -s $sp -g $x > $y/results-score.sc
done

y=$(echo $o | sed -e 's/^.*\///')
tar czf $o/testDISTIQUE-all-BalME-fastme.tar.gz $tmpDIR/* 
rm -r $tmpDIR
