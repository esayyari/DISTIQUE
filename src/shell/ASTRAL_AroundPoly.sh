#!/bin/bash
#set -x
show_help(){
cat << EOF
USAGE: ${0##*/} [-h] [-i SPECIES TREES] [-o OUTPUT FILENAME ] [-l FILE OUTPUT LOG] 
Finds MRL on given input SPECIES TREES, while OUTPUT FILENAME is the destination in which the results will be written. 
EOF
}
I=/dev/null
while getopts "hi:o:I:" opt; do
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
	l)
		I=$OUTARG
		;;
        '?')
                printf "Unknown input option"
                show_help
                ;;
        esac
done
if [ -z "$i" ]; then
        printf "Enter species trees file"
	show_help
        exit 1
fi

if [ -z "$o" ]; then
        printf "Enter an output filename"
	show_help
        exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export FASTMRP=$WS_HOME/mrpmatrix
raxml="$WS_HOME/standard-RAxML/raxmlHPC"
ot=$(dirname "$o")
outdir="$(cd "$ot" && pwd)"
cat $i > $outdir/trees.nwk
in=$outdir/trees.nwk
out=anchored_mrl_tree.nwk


#if [ -s $outdir/$out ]; then
#  echo "Ouput files exists. Refusing to rerun. "
#  exit 0;
#fi
cd $outdir
	#v=$(find $WS_HOME/ASTRAL/ -name "astral.*jar" )
	#vt=$(ls $v)
java -Xmx2500M -jar $WS_HOME/ASTRAL/astral.4.10.3.jar -i $in -o $o 
	

echo ASTRAL Done. Output at: $o


