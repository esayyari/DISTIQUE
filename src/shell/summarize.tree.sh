#!/bin/bash

#set -x
#set -e
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/setup.sh
source $DIR/setupData.sh

show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-s SUMMARIZING  METHOD] [-e SUMMARIING  OPTION] [-d DISTANCE METHOD] [-p PSEUDO COUNT] 
[-o DISTANCE OPTION]  [-r OUTPUT FOLDER ] [-w OVERWRITE FLAG] [GENE TREES FILE]
Estimate species tree using ************ method. GNEE TREES FILE is the file of gene trees concatenated.

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
	-r	OUTPUT FOLDER		The output folder to put the results
	-w 	OVERWRITE FLAG		flag that indicates to overwrite to the OUTPUT FOLDER that already exists or not
EOF
}
m=fm
d=mvr
pc=8
dm=prod
do=1
w=0
if [ $# -lt 1 ]; then
	show_help
	exit 0;
fi  
while getopts "hr:s:e:d:o:p:w:" opt; do
	case $opt in 
	h) 
		show_help
		exit 0;
		;;
	r)
		res_root=$OPTARG
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
	w)	
		w=$OPTARG
		;;
	'?')
		show_help
		exit 1;
		;;
	esac
done
if [ -z $res_root ]; then 
	printf "Enter the output directory\n"
	show_help;
	exit 1
fi

if [ -d $res_root -a $w -eq 0 ]; then
	printf "Enter another output folder. It already exists.\n"
	show_help;
	exit 1
elif [ ! -d $res_root ]; then
	mkdir $res_root
fi

shift "$((OPTIND-1))"
FILE_NAME=$1

printf  "start computing quartets\n"
$DIR/quartet.freq.sh -w $w -f $FILE_NAME -m $dm -p $do -c $pc -o $res_root  >$res_root/distance.d ;
printf "start building species tree\n"
case $m in
fm)
	$WS_HOME/"fastme-2.1.4"/src/fastme -i $res_root/distance.d >$res_root/log.inof;
#	$WS_HOME/"fastme-2.1.4"/src/fastme -i $res_root/distance.d -o estimated_species_tree.tre >$res_root/log.inof;	
	;;
	
pd)
	java -jar $WS_HOME/PhyDstar/PhyDstar.jar -d $d -i $res_root/distance.d
#	$res_root/distance.d_"$d".t > estimated_species_tree.tre
	;;
esac



