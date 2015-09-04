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
done |sed -e "s/^.*: //" | python $WS_LOC_PUTIL/quartetTable.py>$WS_RES/$NM"_qrtTable.q"
cat $WS_RES/$NM"_qrtTable.q" | python $WS_LOC_PUTIL/distance.py>$WS_RES/$NM"_d" 


rm $tmp;
