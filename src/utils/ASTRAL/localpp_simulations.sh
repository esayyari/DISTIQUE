#!/bin/bash
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
set -x
show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-p FILEDIRECTORY] [-g GENETREES FILENAME ] [ -o OUTPUT FOLDER ] [ -x TRUE SPECIES TREE ] [ -s ESTIMATED SPECIES TREE ]
Generates a pool of branches with their posterior probabilities using Astral.

	-h	HELP			     display help and exit
	-p  	FILEDIRECTORY 		     where estimated, and true gene trees, and estimated and true species trees are available 
	-g 	GENETREES FILENAME estimated gene trees filename
	-o 	OUTPUT FOLDER	             The output folder to put the results there.
	-x	TRUE SPECIES TREE	     The true species tree filename.
	-s 	ESTIMATED SPECIES TREE 	     The estimated species tree filename.
EOF
}

if [ $# -ne 5 ]; then 
	show_help
	exit 1 
fi


while getopts "hp:g:o:x:s:" opt; do
        case $opt in
        h)
                show_help
                exit 0;
                ;;
        p)
                p=$OPTARG
                ;;
        g)
                g=$OPTARG
                ;;
        o)
                o=$OPTARG
                ;;
	s)
		s=$OPTARG
		;;
	x)
		x=$OPTARG
		;;
	?)
		show_help
		exit 1;
		;;
	esac
done

version=`grep _versinon $WS_HOME/ASTRAL/main/phylonet/coalescent/CommandLine.java|grep String|sed -e "s/.*= .//g" -e "s/.;//g"`
echo Version $version

TmpFolder=`mktemp -d`;
sp=$TmpFolder/species_tree.trees
spStat=$TmpFolder/species_tree.stat
spScored=$TmpFolder/scored_species_tree.trees
gt=$TmpFolder/gene_trees.trees
tsp=$TmpFolder/true_species_tree.trees
tspStat=$TmpFolder/true_species_tree.stat
tspScored=$TmpFolder/scored_true_species_tree.trees
tmptmp=`mktemp`

cp $p/$s $sp
cp $p/$g $gt
cp $p/$x $tsp

java -Xmx4000M -jar $WS_HOME/ASTRAL/astral.$version.jar -i $gt -q $sp -t 6 > $spScored  2>$spStat ;
test "$?" -ne 0 && echo "error raised in computing bipartitions and quadripartitions info of estimated species tree" && cp $spStat $o && exit 1
echo "bipartitions and quartetpartition info of estimated species tree have been generated"

java -jar $WS_HOME/ASTRAL/astral.$version.jar -i $tsp -q $tsp -t 6 >> $tspScored 2>>$tspStat;
test "$?" -ne 0 && echo "error raised in computing bipartition and quadripartitions info of true species tree" && cp $tspStat $o && exit 1
cat $tspStat | grep "\{" > $tmptmp
awk 'NR %3 == 1' $tmptmp > $tspStat
echo "bipartitions and quadripartitions info of true species tree have been generated"

res=$TmpFolder/ppOfBranches_species_tree
$DIR/extractPPofPoolOfBranches.py -i $spStat -s $tspStat -o $res
echo "pp of branches computed"

tar czvf $out/ppAnalysis-$version.tar.gz $TmpFolder

rm -r $TmpFolder
rm $tmptmp
