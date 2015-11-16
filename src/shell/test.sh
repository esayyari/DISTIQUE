#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/setup.sh

log() {
  printf "\nCalled as: $0: $cmdargs\n\n"
  printf "Start time: "; /bin/date
  printf "Running as user: "; /usr/bin/id
  printf "Running on node: "; /bin/hostname
  printf "Node IP address: "; /bin/hostname -I
  printparams
  printf "\nEnvironment:\n\n"
  printenv | sort
}
show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-g GENE TREES FILE] [-s TRUE SPECIES TREE] [-a ANCHORS] [-f QUARTET TABLE FILE] 
Estimate species tree using DISTIQUE-ANCHORING method. GNEE TREES FILE is the file of gene trees concatenated.

	-g 	GENE TREES FILE		The location of gene trees file.
	-s 	TRUE SPECIES TREE	The location of true species tree for comparison
	-a 	ANCHORS			The anchors
	-o 	OUTPUT FOLDER		The folder to put the results
	-f 	QUARTET TABLE FILE	The file of pre-computed quartet table
EOF
}
if [ $# -lt 1 ]; then
	show_help
	exit 0;
fi  
while getopts "hg:s:a:o:f:" opt; do
	case $opt in 
	h) 
		show_help
		exit 0;
		;;
	g)
		g=$OPTARG
		;;
	s) 
		s=$OPTARG
		;;
	a) 
		a=$OPTARG
		;;
	o)
		o=$OPTARG
		;;
	f)
		f=$OPTARG
		;;
	esac
done
if [ ! -z $f ]; then
python $WS_LOC_PUTIL/testAnchoring.py  "-"g $g "-"s $s  "-"o $o "-"a $a "-"f $f
else
python $WS_LOC_PUTIL/testAnchoring.py  "-"g $g "-"s $s  "-"o $o "-"a $a;
fi

