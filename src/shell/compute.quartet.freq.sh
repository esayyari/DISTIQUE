#!/bin/bash

DIR=$( dirname "${0}")
source $DIR/setup.sh
show_help() {
cat << EOF
USAGE: ${0##*/} [-h]  [-f FILENAME] [-o OUTPUT FOLDER] [-w OVERWRITE FLAG] 
Computes the quartet tables  on top of the gene trees file FILENAME:

	-h	HELP		display help and exit
	-f 	FILENAME 	The gene trees filename
	-o 	OUTPUT FOLDER	The output folder to put the results there.i
	-w 	OVERWRITE FLAG	The flag that indicates to overwrite the output folder if exists or not. 
EOF
}

if [ $# -lt 1 ]; then 
	show_help
	exit 1 
fi
w=0
while getopts "h:f:o:w:" opt; do
        case $opt in
        h)
                show_help
                exit 0;
                ;;
        f)
                filename=$OPTARG
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
qfile="quartet_tmp.q"
for x in `cat $filename | sed -e 's/\[&U\]//'`; do 
  echo -n "$x" >$tmp; 
  $WS_GLB_BIN/quart_bin fancy printQuartets $tmp;
done |sed -e 's/^.*, //'| sed -e 's/://'  |  python $WS_LOC_PUTIL/quartetTable.py> $o/$qfile
 
rm $tmp
