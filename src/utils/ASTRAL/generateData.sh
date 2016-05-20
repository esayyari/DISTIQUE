#!/bin/bash
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
set -x
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
echo Version $version
astral_1000_half=astral-v474-p1-halfresolved.genes1000 
astral_200_half=astral-v474-p1-halfresolved.genes200
astral_50_half=astral-v474-p1-halfresolved.genes50 
astral_1000_true=astral-v474-p1-halfresolved.genes1000 
astral_200_true=astral-v474-p1-halfresolved.genes200 
astral_50_true=astral-v474-p1-halfresolved.genes50 

concat_1000=concatenatedtree.genes1000 
concat_200=concatenatedtree.genes200 
concat_50=concatenatedtree.genes50
 

njst_1000=njst.halfresolved.genes1000 
njst_200=njst.halfresolved.genes200 
njst_50=njst.halfresolved.genes50

concat_1000_true=concatenatedtree.genes1000
concat_200_true=concatenatedtree.genes200
concat_50_true=concatenatedtree.genes50


njst_1000_true=njst.halfresolved.genes1000
njst_200_true=njst.halfresolved.genes200
njst_50_true=njst.halfresolved.genes50

if [ "$m" == "40" ]; then
gt=estimatedgenetre.halfresolved-missing40
tgt=truegenetrees-missing40
elif [ "$m" == "80" ]; then
gt=estimatedgenetre.halfresolved-missing80
tgt=truegenetrees-missing80
else
gt=estimatedgenetre.halfresolved
tgt=truegenetrees
fi

TmpFolder=`mktemp -d`;

A1000half=`mktemp $TmpFolder/astral_1000_half.nwk.XXXXX`;
A200half=`mktemp $TmpFolder/astral_200_half.nwk.XXXXX`;
A50half=`mktemp $TmpFolder/astral_50_half.nwk.XXXXX`;

A1000true=`mktemp $TmpFolder/astral_1000_true.nwk.XXXXX`;
A200true=`mktemp $TmpFolder/astral_200_true.nwk.XXXXX`;
A50true=`mktemp $TmpFolder/astral_50_true.nwk.XXXXX`;

N1000half=`mktemp $TmpFolder/njst_1000.nwk.XXXXX`;
N200half=`mktemp $TmpFolder/njst_200.nwk.XXXXX`;
N50half=`mktemp $TmpFolder/njst_50.nwk.XXXXX`;

C1000half=`mktemp $TmpFolder/concat_1000.nwk.XXXXX`;
C200half=`mktemp $TmpFolder/concat_200.nwk.XXXXX`;
C50half=`mktemp $TmpFolder/concat_50.nwk.XXXXX`;

N1000true=`mktemp $TmpFolder/njst_1000_true.nwk.XXXXX`;
N200true=`mktemp $TmpFolder/njst_200_true.nwk.XXXXX`;
N50true=`mktemp $TmpFolder/njst_50_true.nwk.XXXXX`;

C1000true=`mktemp  $TmpFolder/concat_1000_true.nwk.XXXXX`;
C200true=`mktemp  $TmpFolder/concat_200_true.nwk.XXXXX`;
C50true=`mktemp  $TmpFolder/concat_50_true.nwk.XXXXX`;
echo $TmpFolder

astral_1000_halfStat=`mktemp  $TmpFolder/astral_1000_halfStat.XXXXX`
astral_200_halfStat=`mktemp  $TmpFolder/astral_200_halfStat.XXXXX`
astral_50_halfStat=`mktemp  $TmpFolder/astral_50_halfStat.XXXXX`

astral_1000_trueStat=`mktemp  $TmpFolder/astral_1000_trueStat.XXXXX`
astral_200_trueStat=`mktemp  $TmpFolder/astral_200_trueStat.XXXXX`
astral_50_trueStat=`mktemp  $TmpFolder/astral_50_trueStat.XXXXX`

njst_1000_halfStat=`mktemp  $TmpFolder/njst_1000Stat.XXXXX`
njst_200_halfStat=`mktemp  $TmpFolder/njst_200Stat.XXXXX`
njst_50_halfStat=`mktemp  $TmpFolder/njst_50Stat.XXXXX`


concat_1000_halfStat=`mktemp  $TmpFolder/concat_1000Stat.XXXXX`
concat_200_halfStat=`mktemp  $TmpFolder/concat_200Stat.XXXXX`
concat_50_halfStat=`mktemp  $TmpFolder/concat_50Stat.XXXXX`

njst_1000_trueStat=`mktemp  $TmpFolder/njst_1000_trueStat.XXXXX`
njst_200_trueStat=`mktemp  $TmpFolder/njst_200_trueStat.XXXXX`
njst_50_trueStat=`mktemp  $TmpFolder/njst_50_trueStat.XXXXX`


concat_1000_trueStat=`mktemp  $TmpFolder/concat_1000_trueStat.XXXXX`
concat_200_trueStat=`mktemp  $TmpFolder/concat_200_trueStat.XXXXX`
concat_50_trueStat=`mktemp  $TmpFolder/concat_50_trueStat.XXXXX`


gt1000=`mktemp  $TmpFolder/estimated1000.XXXXX`
gt200=`mktemp  $TmpFolder/estimated200.XXXXX`
gt50=`mktemp  $TmpFolder/estimated50.XXXXX`

tgt1000=`mktemp  $TmpFolder/true1000.XXXXX`
tgt200=`mktemp  $TmpFolder/true200.XXXXX`
tgt50=`mktemp  $TmpFolder/true50.XXXXX`


head -n 1000 $p/$gt >$gt1000
head -n 200 $gt1000 > $gt200
head -n 50 $gt1000 > $gt50



head -n 1000 $p/$tgt > $tgt1000
head -n 200 $tgt1000 > $tgt200
head -n 50 $tgt1000 > $tgt50

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt1000 -q $p/$astral_1000_half -t 6 > $A1000half  2>$astral_1000_halfStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt200 -q $p/$astral_200_half -t 6 >   $A200half 2>$astral_200_halfStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt50 -q $p/$astral_50_half -t 6 > $A50half 2>$astral_50_halfStat ;

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt1000 -q $p/$astral_1000_true -t 6 > $A1000true 2> $astral_1000_trueStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt200 -q $p/$astral_200_true -t 6 > $A200true 2> $astral_200_trueStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt50 -q $p/$astral_50_true -t 6 > $A50true 2> $astral_50_trueStat ;

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt1000 -q $p/$njst_1000 -t 6 > $N1000half 2> $njst_1000_halfStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt200 -q $p/$njst_200 -t 6 > $N200half 2> $njst_200_halfStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt50 -q $p/$njst_50 -t 6 > $N50half 2> $njst_50_halfStat ;

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt1000 -q $p/$concat_1000 -t 6 > $C1000half 2> $concat_1000_halfStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt200 -q $p/$concat_200 -t 6 > $C200half 2> $concat_200_halfStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt50 -q $p/$concat_50 -t 6 > $C50half 2> $concat_50_halfStat ;


java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt1000 -q $p/$njst_1000_true -t 6 > $N1000true 2> $njst_1000_trueStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt200 -q $p/$njst_200_true -t 6 > $N200true 2> $njst_200_trueStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt50 -q $p/$njst_50_true -t 6 > $N50true 2> $njst_50_trueStat ;

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt1000 -q $p/$concat_1000_true -t 6 > $C1000true 2> $concat_1000_trueStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt200 -q $p/$concat_200_true -t 6 > $C200true 2> $concat_200_trueStat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tgt50 -q $p/$concat_50_true -t 6 > $C50true 2> $concat_50_trueStat ;



A_trueSp200half=`mktemp  $TmpFolder/astral_200_half_sp.nwk.XXXXX`;
A_trueSp1000half=`mktemp  $TmpFolder/astral_1000_half_sp.nwk.XXXXX`;
A_trueSp50half=`mktemp  $TmpFolder/astral_50_half_sp.nwk.XXXXX`;
A_trueSp1000true=`mktemp  $TmpFolder/astral_1000_true_sp.nwk.XXXXX`;
A_trueSp200true=`mktemp  $TmpFolder/astral_200_true_sp.nwk.XXXXX`;
A_trueSp50true=`mktemp  $TmpFolder/astral_50_true_sp.nwk.XXXXX`;
#N_trueSp1000half=`mktemp  $TmpFolder njst_1000_sp.nwk.XXXXX`;
#N_trueSp200half=`mktemp  $TmpFolder njst_200_sp.nwk.XXXXX`;
#N_trueSp50half=`mktemp  $TmpFolder njst_50_sp.nwk.XXXXX`;
#C_trueSp1000half=`mktemp  $TmpFolder concat_1000_sp.nwk.XXXXX`;
#C_trueSp200half=`mktemp  $TmpFolder concat_200_sp.nwk.XXXXX`;
#C_trueSp50half=`mktemp  $TmpFolder concat_50_sp.nwk.XXXXX`;
echo $TmpFolder
astral_1000_half_trueSp_Stat=`mktemp  $TmpFolder/astral_1000_half_sp_Stat.XXXXX`
astral_200_half_trueSp_Stat=`mktemp  $TmpFolder/astral_200_half_sp_Stat.XXXXX`
astral_50_half_trueSp_Stat=`mktemp  $TmpFolder/astral_50_half_sp_Stat.XXXXX`

astral_1000_true_trueSp_Stat=`mktemp  $TmpFolder/astral_1000_true_sp_Stat.XXXXX`
astral_200_true_trueSp_Stat=`mktemp  $TmpFolder/astral_200_true_sp_Stat.XXXXX`
astral_50_true_trueSp_Stat=`mktemp  $TmpFolder/astral_50_true_sp_Stat.XXXXX`
#njst_1000_half_trueSp_Stat=`mktemp  $TmpFolder njst_1000_sp_Stat.XXXXX`
#njst_200_half_trueSp_Stat=`mktemp  $TmpFolder njst_200_sp_Stat.XXXXX`
#njst_50_half_trueSp_Stat=`mktemp  $TmpFolder njst_50_sp_Stat.XXXXX`

#concat_1000_half_trueSp_Stat=`mktemp  $TmpFolder concat_1000_sp_Stat.XXXXX`
#concat_200_half_trueSp_Stat=`mktemp  $TmpFolder concat_200_sp_Stat.XXXXX`
#concat_50_half_trueSp_Stat=`mktemp  $TmpFolder concat_50_sp_Stat.XXXXX`



java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -q $s -i $gt1000 -t 6 > $A_trueSp1000half 2> $astral_1000_half_trueSp_Stat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -q $s -i $gt200 -t 6 > $A_trueSp200half 2> $astral_200_half_trueSp_Stat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -q $s -i $gt50  -t 6 > $A_trueSp50half 2> $astral_50_half_trueSp_Stat ;

java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -q $s -i $tgt1000 -t 6 > $A_trueSp1000true 2> $astral_1000_true_trueSp_Stat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -q $s -i $tgt200 -t 6 > $A_trueSp200true 2> $astral_200_true_trueSp_Stat ;
java -Xmx2000M -jar $WS_HOME/ASTRAL/astral.$version.jar -q $s -i $tgt50  -t 6 > $A_trueSp50true 2> $astral_50_true_trueSp_Stat ;





echo "bipartition and quartetpartition info of astral and njst trees has been generated"


TmpSpStat=`mktemp -p $TmpFolder spTreeStat.XXXXX`;
TmpSpTree=`mktemp -p $TmpFolder spTreePP.XXXXX`;
tmptmp=`mktemp -p $TmpFolder`;
java -jar $WS_HOME/ASTRAL/astral.$version.jar -i $s -q $s -t 6 >> $TmpSpTree 2>>$TmpSpStat;
cat $TmpSpStat | grep "\{" > $tmptmp
awk 'NR %3 == 1' $tmptmp > $TmpSpStat
echo "bipartition and quartetpartition info of species tree have been generated"
res_astral_1000_half=`mktemp  $TmpFolder/ppOfBranches_astral_1000_half.XXXXX`;
res_astral_200_half=`mktemp  $TmpFolder/ppOfBranches_astral_200_half.XXXXX`;
res_astral_50_half=`mktemp  $TmpFolder/ppOfBranches_astral_50_half.XXXXX`;

res_astral_1000_true=`mktemp  $TmpFolder/ppOfBranches_astral_1000_true.XXXXX`;
res_astral_200_true=`mktemp  $TmpFolder/ppOfBranches_astral_200_true.XXXXX`;
res_astral_50_true=`mktemp  $TmpFolder/ppOfBranches_astral_50_true.XXXXX`;

res_njst_1000_half=`mktemp  $TmpFolder/ppOfBranches_njst_1000_half.XXXXX`;
res_njst_200_half=`mktemp  $TmpFolder/ppOfBranches_njst_200_half.XXXXX`;
res_njst_50_half=`mktemp  $TmpFolder/ppOfBranches_njst_50_half.XXXXX`;

res_concat_1000_half=`mktemp  $TmpFolder/ppOfBranches_concat_1000_half.XXXXX`;
res_concat_200_half=`mktemp  $TmpFolder/ppOfBranches_concat_200_half.XXXXX`;
res_concat_50_half=`mktemp  $TmpFolder/ppOfBranches_concat_50_half.XXXXX`;

res_sp_1000_half=`mktemp  $TmpFolder/ppOfBranches_sp_1000_half.XXXXX`;
res_sp_200_half=`mktemp  $TmpFolder/ppOfBranches_sp_200_half.XXXXX`;
res_sp_50_half=`mktemp  $TmpFolder/ppOfBranches_sp_50_half.XXXXX`;

res_sp_1000_true=`mktemp  $TmpFolder/ppOfBranches_sp_1000_true.XXXXX`;
res_sp_200_true=`mktemp  $TmpFolder/ppOfBranches_sp_200_true.XXXXX`;
res_sp_50_true=`mktemp  $TmpFolder/ppOfBranches_sp_50_true.XXXXX`;

res_njst_1000_true=`mktemp  $TmpFolder/ppOfBranches_njst_1000_true.XXXXX`;
res_njst_200_true=`mktemp  $TmpFolder/ppOfBranches_njst_200_true.XXXXX`;
res_njst_50_true=`mktemp  $TmpFolder/ppOfBranches_njst_50_true.XXXXX`;

res_concat_1000_true=`mktemp  $TmpFolder/ppOfBranches_concat_1000_true.XXXXX`;
res_concat_200_true=`mktemp  $TmpFolder/ppOfBranches_concat_200_true.XXXXX`;
res_concat_50_true=`mktemp  $TmpFolder/ppOfBranches_concat_50_true.XXXXX`;

$DIR/extractPPofPoolOfBranches.py -i $astral_1000_halfStat -s $TmpSpStat -o $res_astral_1000_half
$DIR/extractPPofPoolOfBranches.py -i $astral_200_halfStat -s $TmpSpStat -o $res_astral_200_half
$DIR/extractPPofPoolOfBranches.py -i $astral_50_halfStat -s $TmpSpStat -o $res_astral_50_half

$DIR/extractPPofPoolOfBranches.py -i $astral_1000_trueStat -s $TmpSpStat -o $res_astral_1000_true
$DIR/extractPPofPoolOfBranches.py -i $astral_200_trueStat -s $TmpSpStat -o $res_astral_200_true
$DIR/extractPPofPoolOfBranches.py -i $astral_50_trueStat -s $TmpSpStat -o $res_astral_50_true

$DIR/extractPPofPoolOfBranches.py -i $njst_1000_halfStat -s $TmpSpStat -o $res_njst_1000_half
$DIR/extractPPofPoolOfBranches.py -i $njst_200_halfStat -s $TmpSpStat -o $res_njst_200_half
$DIR/extractPPofPoolOfBranches.py -i $njst_50_halfStat -s $TmpSpStat -o $res_njst_50_half


$DIR/extractPPofPoolOfBranches.py -i $concat_1000_halfStat -s $TmpSpStat -o $res_concat_1000_half
$DIR/extractPPofPoolOfBranches.py -i $concat_200_halfStat -s $TmpSpStat -o $res_concat_200_half
$DIR/extractPPofPoolOfBranches.py -i $concat_50_halfStat -s $TmpSpStat -o $res_concat_50_half

$DIR/extractPPofPoolOfBranches.py -i $njst_1000_trueStat -s $TmpSpStat -o $res_njst_1000_true
$DIR/extractPPofPoolOfBranches.py -i $njst_200_trueStat -s $TmpSpStat -o $res_njst_200_true
$DIR/extractPPofPoolOfBranches.py -i $njst_50_trueStat -s $TmpSpStat -o $res_njst_50_true


$DIR/extractPPofPoolOfBranches.py -i $concat_1000_trueStat -s $TmpSpStat -o $res_concat_1000_true
$DIR/extractPPofPoolOfBranches.py -i $concat_200_trueStat -s $TmpSpStat -o $res_concat_200_true
$DIR/extractPPofPoolOfBranches.py -i $concat_50_trueStat -s $TmpSpStat -o $res_concat_50_true

$DIR/extractPPofPoolOfBranches.py -i $astral_1000_half_trueSp_Stat -s $TmpSpStat -o $res_sp_1000_half
$DIR/extractPPofPoolOfBranches.py -i $astral_200_half_trueSp_Stat -s $TmpSpStat -o $res_sp_200_half
$DIR/extractPPofPoolOfBranches.py -i $astral_50_half_trueSp_Stat -s $TmpSpStat -o $res_sp_50_half

$DIR/extractPPofPoolOfBranches.py -i $astral_1000_true_trueSp_Stat -s $TmpSpStat -o $res_sp_1000_true
$DIR/extractPPofPoolOfBranches.py -i $astral_200_true_trueSp_Stat -s $TmpSpStat -o $res_sp_200_true
$DIR/extractPPofPoolOfBranches.py -i $astral_50_true_trueSp_Stat -s $TmpSpStat -o $res_sp_50_true


echo "pp of branches computed"
cp $p/$gt $TmpFolder
cp $p/$tgt $TmpFolder
tar czvf $out/astral-BUG-missing"$m".tar.gz $TmpFolder

rm -r $TmpFolder
