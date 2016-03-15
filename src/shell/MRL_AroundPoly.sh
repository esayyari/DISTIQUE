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


if [ -s $outdir/$out ]; then
  echo "Ouput files exists. Refusing to rerun. "
  exit 0;
fi
cd $outdir
tmp=`mktemp $outdir/mrpmaptrix.$out.XXXXX`
java -jar $FASTMRP/mrp.jar $in $tmp  PHYLIP -randomize > $I 2>&1

#rm -f RAxML_*$out*
$raxml -m BINCAT -s $tmp -n $out -N 2 -p $RANDOM > $I 2>&1

#test "$?" != "0" && exit 1

mv RAxML_bestTree.$out* $o
mv RAxML_info.$out $outdir/mrl.log

rm RAxML_*$out* 


#echo MRL Done. Output at: $o


