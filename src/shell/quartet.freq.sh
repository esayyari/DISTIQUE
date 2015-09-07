#!/bin/bash

DIR=$( dirname "${0}")
source $DIR/setup.sh
show_help() {
cat << EOF
USAGE: ${0##*/} [-h] [-m DISTANCE METHOD] [-c PSEUDO COUNT] [-p PERCENTILE] [-f FILENAME] 
Computes the distances using the python function distance.py on top of the gene trees file FILENAME:

	-h	HELP		display help and exit
	-c  	PSEUDO COUNT 	pseudo count that is used by the distance.py to handle empty bin problem
				It should be a positive number that defines the power of 10, ie 1e-pc
	-m	DISTANCE METHOD The method to compute distances between taxa. It could be min, prod, minavg,
				and minmed
	-f 	FILENAME 	The gene trees filename
	-p	PERCENTILE	The percentile of the data that in methods minavg, and minmed will be used to 
				compute distance. The default is 1.
EOF
}

if [ $# -lt 1 ]; then 
	show_help 
fi
m=prod
c=8
p=1

while getopts "hc:m:f:p:" opt; do
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
        '?')
                show_help
                exit 1;
                ;;
        esac
done
shift "$((OPTIND-1))"


tmp=`mktemp`
NM=$(basename "${filename}")
WS_RES=$( dirname "${filename}")
for x in `cat $filename`; do 
  echo -n "$x" >$tmp; 
  $WS_GLB_BIN/quart_bin fancy printQuartets $tmp;
done |sed -e "s/^.*: //" | python $WS_LOC_PUTIL/quartetTable.py > $WS_RES/$NM"_qttable.q"
 
python $WS_LOC_PUTIL/distance.py -m $m -c $c -p $p -f $WS_RES/$NM"_qttable.q" >$WS_RES/$NM"_d" 

rm $tmp;
