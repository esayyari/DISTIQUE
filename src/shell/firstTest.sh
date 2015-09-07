#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/setup.sh
source $DIR/setupData.sh

show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-s SUMMARIZING  METHOD] [-e SUMMARIING  OPTION] [-d DISTANCE METHOD] [-p PSEUDO COUNT] 
[-o DISTANCE OPTION] [-c SPECIES TREE FILE] [GENE TREES FILE]
Estimate species tree using ************ method. GNEE TREES FILE is the file of gene trees concatenated 
and SPECIES TREE FILE is the species tree file.

	-c	SPECIES TREE FILE	indicates to compare the tree provided by our method with the SPECIES TREE FILE.
	-h      			display help and exit
	-s 	SUMMARIZING METHOD	the method to produce estimated species tree based on distance matrix.
					We use fastME (fm) and PhyDstar (pd). Note that thsse two pacakges should be installed
					beforehand. PhyDstar should be located under WS_HOME path (described at setup.sh) 
					under the folder PhyDstar.
	-e	SUMMARIZING OPTION	Indicates wich method to use to infer species tree based on distance matrix if you
					use PhyDstar. Options are NJ, BioNJ, MVR, and UNJ. For more details look at: 
					http://www.atgc-montpellier.fr/phyd/usersguide.php
	-d 	DISTANCE METHOD		the method to compute distances. The methods are: min, prod, minavg, minmed
	-o 	DISTANCE OPTION		If you choose the methods minavg or minmed you could define the percentile of which
					the average or median will be computed as the distance
	-p 	PSEUDO COUNT		The pseudo count to avoid empty bin problem. This should be a positive integer, 
					which indicates	power of 10 i.e. 10^-ps.
EOF
}
m=fm
d=mvr
pc=8
dm=prod
do=1

if [ $# -lt 1 ]; then
	show_help
	exit 0;
fi  
while getopts "hc:s:e:d:o:p:" opt; do
	case $opt in 
	h) 
		show_help
		exit 0;
		;;
	c)
		sp_tree=$OPTARG
		echo $sp_tree
		;;
	s) 
		m=$OPTARG
		;;
	e) 
		d=$OPTARG
		;;
	o)
		do=$OPTARG
		;;
	p)  	
		pc=$OPTARG
		;;
	d)	
		dm=$OPTARG
		;;
	'?')
		show_help
		exit 1;
		;;
	esac
done
shift "$((OPTIND-1))"
FILE_NAME=$1
WS_RES=$(dirname "${FILE_NAME}")
NM=$(basename "${FILE_NAME}")
printf  "start computing quartets\n"
$DIR/quartet.freq.sh -f $FILE_NAME -m $dm -p $do -c $pc >$WS_RES/$NM"_d" ;
echo $m
printf "start building species tree\n"
case $m in
fm)
	fastme -i $WS_RES/$NM"_d">$WS_RES/log.inof;
	;;
	
pd)
	java -jar $WS_HOME/PhyDstar/PhyDstar.jar -d $d -i $WS_RES/$NM"_d"
	;;
esac

if [ ! -z "$sp_tree" ]; then
	if [ "$m" == "pd" ]; then
		case $d in
		nj)
			echo "`head -n 1 $WS_RES/$NM"_d_nj.t"`"> $WS_RES/$NM"_d_nj.t"
			echo "comparing files `echo -n $sp_tree` and `echo -n $WS_RES/$NM"_d_nj.t"`"
			$WS_GLB_SH/compareTrees.missingBranch $sp_tree $WS_RES/$NM"_d_nj.t">$WS_RES/$NM"_d_nj.t_res_stat.txt"
			;;
		bionj)
			echo "`head -n 1 $WS_RES/$NM"_d_bionj.t"`"> $WS_RES/$NM"_d_bionj.t"
			echo "comparing files `echo -n $sp_tree` and `echo -n $WS_RES/$NM"_d_bionj.t"`"
			$WS_GLB_SH/compareTrees.missingBranch $sp_tree $WS_RES/$NM"_d_bionj.t">$WS_RES/$NM"_d_bionj.t_res_stat.txt"
			;;	
		mvr)
			echo "`head -n 1 $WS_RES/$NM"_d_mvr.t"`"> $WS_RES/$NM"_d_mvr.t"
			echo "comparing files `echo -n $sp_tree` and `echo -n $WS_RES/$NM"_d_mvr.t"`"
			$WS_GLB_SH/compareTrees.missingBranch $sp_tree $WS_RES/$NM"_d_mvr.t">$WS_RES/$NM"_d_mrv.t_res_stat.txt"
			;;
		unj)
			echo "`head -n 1 $WS_RES/$NM"_d_unj.t"`"> $WS_RES/$NM"_d_unj.t"
			echo "comparing files `echo -n $sp_tree` and `echo -n $WS_RES/$NM"_d_unj.t"`"
			$WS_GLB_SH/compareTrees.missingBranch $sp_tree $WS_RES/$NM"_d_unj.t">$WS_RES/$NM"_d_unj.t_res_stat.txt"
			;;
		?)
			printf "Unknown option. -d should be either nj, or bionj, or mvr, or unj."
			exit 1;
			;;
		esac
	else
 				
		echo "`head -n 1 $WS_RES/$NM"_d_fastme_tree.nwk"`"> $WS_RES/$NM"_d_fastme_tree.nwk"
		echo "comparing files `echo -n $sp_tree` and `echo -n $WS_RES/$NM"_d_fastme_tree.nwk"`"
		$WS_GLB_SH/compareTrees.missingBranch $sp_tree $WS_RES/$NM"_d_fastme_tree.nwk">$WS_RES/$NM"_d_fastme_tree.nwk_res_stat.txt"
	fi
fi


