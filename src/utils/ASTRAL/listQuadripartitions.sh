#!/bin/bash

function show_help(){
	cat << EOF 
	USAGE: ${0##*/} [-h] [-s SPECIES TREE] [-o OUTPUT]
	-h	HELP show help	and exit
	-o 	OUTPUT FOLDER	The output folder to put the results there.
	-s 	SPECIES TREE 	The true species tree to compare the results.
EOF
}

if [ $# -lt 1 ]; then 
	show_help
	exit 1 
fi


while getopts "ho:s:" opt; do
        case $opt in
        h)
                show_help
                exit 0;
                ;;
        s)
                s=$OPTARG
                ;;
        o)
                o=$OPTARG
                ;;
	esac
done
tmp=`mktemp`
astral=$(basename $(find $WS_HOME/ASTRAL/ -name "astral*.jar" ))
java -jar $WS_HOME/ASTRAL/$astral -i $s -q $s -t 7  > $tmp 2>&1
grep "^{" $tmp > $o/quadripartitions.q
