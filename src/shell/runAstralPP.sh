#!/bin/bash
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
set -x
show_help(){
cat << EOF
USAGE: ${0##*/} [-h] [-i DATA PATH] [-o OUTPATH] 
Running ASTRAL on input gene tree with pp instead of usual ASTRAL. The output will be at OUTPATH.
EOF
}
while getopts "hi:o:" opt; do
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
WS_LOC_UTIL=$WS_HOME/ASTRAL/
tmpDIR=`mktemp -d`
echo $tmpDIR

gt=`mktemp  ${tmpDIR}/genetre.XXXXX`   || exit 1
out=`mktemp ${tmpDIR}/astral-pp.nwk.XXXXX` || exit 1
nw_prune $i 1 2 > $gt 
java -Xmx3000M -jar $WS_LOC_UTIL/astral.4.10.4.jar  -i $gt -o $out -y -t 12  2>/dev/null
out=`mktemp ${tmpDIR}/astral.nwk.XXXXX` || exit 1
java -Xmx3000M -jar $WS_LOC_UTIL/astral.4.10.4.jar  -i $gt -o $out 2>/dev/null
out=`mktemp ${tmpDIR}/astrid.nwk.XXXXX` || exit 1
python $WS_LOC_UTIL/../ASTRID/ASTRID -i $gt -o $out -m fastme2 
tar czvf $o/results.tar.gz $tmpDIR/*
rm -r $tmpDIR
