#!/bin/bash

DIR=$( dirname "${0}")
source $DIR/setup.sh
show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-m DISTANCE METHOD] [-c PSEUDO COUNT] [-p PERCENTILE] [-f FILENAME] [-o OUTPUT FOLDER] [-w OVERWRITE FLAG] 
Computes the distances using the python function distance.py on top of the gene trees file FILENAME:

	-h	HELP		display help and exit
	-c  	PSEUDO COUNT 	pseudo count that is used by the distance.py to handle empty bin problem
				It should be a positive number that defines the power of 10, ie 1e-pc
	-m	DISTANCE METHOD The method to compute distances between taxa. It could be min, prod, minavg,
				and minmed
	-f 	FILENAME 	The gene trees filename
	-p	PERCENTILE	The percentile of the data that in methods minavg, and minmed will be used to 
				compute distance. The default is 1.
	-o 	OUTPUT FOLDER	The output folder to put the results there.i
	-w 	OVERWRITE FLAG	The flag that indicates to overwrite the output folder if exists or not. 
EOF
}

if [ $# -lt 1 ]; then 
	show_help
	exit 1 
fi
m=prod
c=8
p=1
w=0
while getopts "hc:m:f:p:o:w:" opt; do
        case $opt in
        h)
                show_help
                exit 0;
                ;;
        c)
                c=$OPTARG
                ;;
        f)
                filename=$OPTARG
                ;;
        p)
                p=$OPTARG
                ;;
        m)
                m=$OPTARG
                ;;
	o)
		o=$OPTARG
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
shift "$((OPTIND-1))"

if [ -z $o ]; then
printf "Enter the output directory"
exit 1
fi

if [ -d $o ] && [ $w -eq 0 ]; then 
printf "The output directory already exists, enter another output folder"
exit 1
elif [ ! -d $o ]; then
mkdir $o
fi
tmp=`mktemp` 
qfile="quartets.q"
for x in `cat $filename`; do 
  echo -n "$x" >$tmp; 
  $WS_GLB_BIN/quart_bin fancy printQuartets $tmp;
done |sed -e 's/^.*, //'| sed -e 's/://' > $o/quartets_original.q;
cat $o/quartets_original.q  |  python $WS_LOC_PUTIL/quartetTable.py> $o/$qfile
 
#python $WS_LOC_PUTIL/distance.py -m $m -c $c -p $p -f $o/$qfile
rm $tmp
