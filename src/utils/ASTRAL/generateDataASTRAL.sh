#!/bin/bash
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
#set -x
show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-g GENETREES] [-t CONSENSUS THRESHOLD] [-m MISSING RATE] [-o OUTPATH] [ -s TRUE SPECIES TREE]
Generates a pool of branches with their posterior probabilities using Astral.

	-h	HELP		display help and exit
	-g  	GENETREES 	Path where gene trees available and estimated species trees available
	-t	THRESHOLD 	The consensus threshold to generate randomly resolved samples.
	-m 	MISSING RATE 	Number of missing taxa per each gene.
	-o 	OUTPUT FOLDER	The output folder to put the results there.
	-s 	SPECIES TREE 	The true species tree to compare the results.
EOF
}

if [ $# -lt 1 ]; then 
	show_help
	exit 1 
fi

t=0.2;
m=0;

while getopts "hg:t:m:o:s:" opt; do
        case $opt in
        h)
                show_help
                exit 0;
                ;;
        g)
                p=$OPTARG
                ;;
        t)
                t=$OPTARG
                ;;
        m)
                m=$OPTARG
                ;;
        o)
                out=$OPTARG
                ;;
	s)
		s=$OPTARG
		;;
	?)
		show_help
		exit 1;
		;;
	esac
done
#rm $out/*
version=`grep _versinon $WS_HOME/ASTRAL/main/phylonet/coalescent/CommandLine.java|grep String|sed -e "s/.*= .//g" -e "s/.;//g"`
#echo Version $version
astral_1000_half=astral/Best.tre
astral_1000_true=astral/Best.tre

gt=genetrees.gt
tgt="true-genetrees.gt"

TmpFolder=`mktemp -d`;

A1000half=`mktemp -p $TmpFolder astral_1000_half.nwk.XXXXX`;
A1000true=`mktemp -p $TmpFolder astral_1000_true.nwk.XXXXX`;
#echo $TmpFolder

astral_1000_halfStat=`mktemp -p $TmpFolder astral_1000_halfStat.XXXXX`
astral_1000_trueStat=`mktemp -p $TmpFolder astral_1000_trueStat.XXXXX`


gt1000=`mktemp -p $TmpFolder estimated1000.XXXXX`
tgt1000=`mktemp -p $TmpFolder true1000.XXXXX`
head -n 1000 $p/$gt >$gt1000
head -n 1000 $p/$tgt >$tgt1000

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt1000 -q $p/$astral_1000_half -t 4 > $A1000half  2>$astral_1000_halfStat ;

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt1000 -q $p/$astral_1000_true -t 4 > $A1000true  2>$astral_1000_trueStat ;


A_trueSp1000half=`mktemp -p $TmpFolder astral_1000_half_sp.nwk.XXXXX`;
A_trueSp1000true=`mktemp -p $TmpFolder astral_1000_true_sp.nwk.XXXXX`;

#echo $TmpFolder
astral_1000_half_trueSp_Stat=`mktemp -p $TmpFolder astral_1000_half_sp_Stat.XXXXX`

astral_1000_true_trueSp_Stat=`mktemp -p $TmpFolder astral_1000_true_sp_Stat.XXXXX`

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -q $s -i $gt1000 -t 4 > $A_trueSp1000half 2> $astral_1000_half_trueSp_Stat ;

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -q $s -i $tgt1000 -t 4 > $A_trueSp1000true 2> $astral_1000_true_trueSp_Stat ;


#echo "bipartition and quartetpartition info of astral and njst trees has been generated"


TmpSpStat=`mktemp -p $TmpFolder spTreeStat.XXXXX`;
TmpSpTree=`mktemp -p $TmpFolder spTreePP.XXXXX`;
tmptmp=`mktemp -p $TmpFolder`;
java -jar $WS_HOME/ASTRAL/astral.$version.jar -i $s -q $s -t 4 >> $TmpSpTree 2>>$TmpSpStat;
cat $TmpSpStat | grep "\{" > $tmptmp
awk 'NR %3 == 1' $tmptmp > $TmpSpStat
#echo "bipartition and quartetpartition info of species tree have been generated"
res_astral_1000_half=`mktemp -p $TmpFolder ppOfBranches_astral_1000_half.XXXXX`;
res_astral_1000_true=`mktemp -p $TmpFolder ppOfBranches_astral_1000_true.XXXXX`;

res_sp_1000_half=`mktemp -p $TmpFolder ppOfBranches_sp_1000_half.XXXXX`;
$DIR/extractPPofPoolOfBranches.py -i $astral_1000_halfStat -s $TmpSpStat -o $res_astral_1000_half
$DIR/extractPPofPoolOfBranches.py -i $astral_1000_trueStat -s $TmpSpStat -o $res_astral_1000_true
$DIR/extractPPofPoolOfBranches.py -i $astral_1000_half_trueSp_Stat -s $TmpSpStat -o $res_sp_1000_half
#cat $A_trueSp1000half | grep "^(" > $TmpFolder/sp_scored_half.nwk
#cat $A_trueSp1000true | grep "^(" > $TmpFolder/sp_scored_true.nwk

#$WS_HOME/global/src/shell/compareTrees $s $TmpFolder/sp_scored_half.nwk > $TmpFolder/astral.bl 
#echo "pp of branches computed"
cp $p/$gt $TmpFolder
cp $p/$tgt $TmpFolder
tar czf $out/ppAnalysis.tar.gz $TmpFolder

rm -r $TmpFolder
