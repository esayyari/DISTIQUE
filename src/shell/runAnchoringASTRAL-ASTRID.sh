#!/bin/bash
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
set -x
show_help(){
cat << EOF
USAGE: ${0##*/} [-h] [-i DATA PATH] [-o OUTPATH] [-n NUMROUNDS] [-v VERSION]
Runs DISTIQUE Anchoring on given gene tree available at DATA PATH. The results will be written at OUTPATH. NUMROUNDS is number of rounds for sampling around each polytomy, while VERSION indicates which version of the code to use.
EOF
}
n=2
v=3
z=B
while getopts "hi:o:n:v:z:" opt; do
        case $opt in
        h)
                show_help
                exit 0
                ;;
        i)
                i=$OPTARG
                ;;
        o)
                o=$OPTARG
                ;;
	n)
		n=$OPTARG
		;;
	v)
		v=$OPTARG
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

if [ -z "$i" ]; then
        printf "Enter data path\n"
        show_help
        exit 1
fi
if [ -z "$o" ]; then
        printf "Enter output path\n"
        show_help
        exit 1
fi

WS_LOC_UTIL=$WS_HOME/ASTRID/
echo $WS_LOC_UTIL
tmpDIR=`mktemp -d`
echo $tmpDIR

gt1000half=`mktemp  ${tmpDIR}/genetrees1000.half.XXXXX` || exit 1
gt200half=`mktemp  ${tmpDIR}/genetrees200.half.XXXXX`  || exit 1
gt50half=`mktemp  ${tmpDIR}/genetrees50.half.XXXXXX`  || exit 1
gt1000true=`mktemp  ${tmpDIR}/genetrees1000.true.XXXXX` || exit 1
gt200true=`mktemp  ${tmpDIR}/genetrees200.true.XXXXX`  || exit 1
gt50true=`mktemp  ${tmpDIR}/genetrees50.true.XXXXX`   || exit 1

head -n 1000 $i/estimatedgenetre.halfresolved > $gt1000half
head -n 200  $i/estimatedgenetre.halfresolved > $gt200half
head -n 50   $i/estimatedgenetre.halfresolved > $gt50half
head -n 1000 $i/truegenetrees > $gt1000true
head -n 200  $i/truegenetrees > $gt200true
head -n 50   $i/truegenetrees > $gt50true
sp=$i/s_tree.trees

res50true=`mktemp -d $tmpDIR/50genes.true.XXXXX` || exit 1
res200true=`mktemp -d $tmpDIR/200genes.true.XXXXX` || exit 1
res1000true=`mktemp -d $tmpDIR/1000genes.true.XXXXX` || exit 1

res50half=`mktemp -d $tmpDIR/50genes.half.XXXXX` || exit 1
res200half=`mktemp -d $tmpDIR/200genes.half.XXXXX` || exit 1
res1000half=`mktemp -d $tmpDIR/1000genes.half.XXXXX` || exit 1

if [ -s $gt1000half ]; then
cd $res1000half
/usr/bin/time -p python $WS_LOC_UTIL/ASTRID -i $gt1000half -o $res1000half/distance.d_astrid_tree.nwk -m fastme2 -c CACHE.csv > $res1000half/results.log 2>&1
echo "working on $i/estimatedgenetre.halfresolved1000 has been finished!"
fi
if [ -s $gt200half ]; then
cd $res200half
/usr/bin/time -p python $WS_LOC_UTIL/ASTRID -i $gt200half  -o $res200half/distance.d_astrid_tree.nwk -m fastme2 -c CACHE.csv > $res200half/results.log 2>&1
echo "working on $i/estimatedgenetre.halfresolved200 has been finished!"
fi
if [ -s $gt50half ]; then
cd $res50half
/usr/bin/time -p python $WS_LOC_UTIL/ASTRID -i $gt50half   -o $res50half/distance.d_astrid_tree.nwk -m fastme2 -c CACHE.csv > $res50half/results.log 2>&1
echo "working on $i/estimatedgenetre.halfresolved50 has been finished!"
fi
if [ -s $res1000true ]; then
cd $res1000true
/usr/bin/time -p python $WS_LOC_UTIL/ASTRID -i $gt1000true -o $res1000true/distance.d_astrid_tree.nwk -m fastme2 -c CACHE.csv > $res1000true/results.log 2>&1
echo "working on $i/truegenetre1000 has been finished!"
fi
if [ -s $res200true ]; then
cd $res200true
/usr/bin/time -p python $WS_LOC_UTIL/ASTRID -i $gt200true  -o $res200true/distance.d_astrid_tree.nwk -m fastme2 -c CACHE.csv > $res200true/results.log 2>&1
echo "working on $i/truegenetre200 has been finished!"
fi
if [ -s $res50true ]; then
cd $res50true
/usr/bin/time -p python $WS_LOC_UTIL/ASTRID -i $gt50true   -o $res50true/distance.d_astrid_tree.nwk -m fastme2 -c CACHE.csv > $res50true/results.log 2>&1
echo "working on $i/truegenetre50 has been finished!"
fi
for x in `find $tmpDIR -name "distance.d_astrid_tree.nwk"`; do
	y=$(dirname $x)
	$WS_HOME/DISTIQUE/src/shell/compare.tree.sh -s $sp -g $x > $y/results-astrid.score
done
y=$(echo $o | sed -e 's/^.*\///')
tar czf $o/njst-astrid.tar.gz $tmpDIR/* 
rm -r $tmpDIR
