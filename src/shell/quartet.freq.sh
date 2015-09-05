#!/bin/bash

test $# == 1 || { echo USAGE: $0 file_with_newick_trees; exit 1; }



DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/setup.sh

tmp=`mktemp`
NM=$(basename "${1}")
WS_RES=$( dirname "${1}")
for x in `cat $1`; do 
  echo -n "$x" >$tmp; 
  $WS_GLB_BIN/quart_bin fancy printQuartets $tmp;
done |sed -e "s/^.*: //" | python $WS_LOC_PUTIL/quartetTable.py > $WS_RES/$NM"_qttable.q"
printf $WS_RES/$NM"_qttable.q\n"
python $WS_LOC_PUTIL/distance.py -m min -c 1e-12 -p 1 -f $WS_RES/$NM"_qttable.q" >$WS_RES/$NM"_d" 
 

rm $tmp;
