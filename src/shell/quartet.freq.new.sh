#!/bin/bash

test $# == 1 || { echo USAGE: $0 file_with_newick_trees; exit 1; }

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/setup.sh

tmp=`mktemp`

for x in `cat $1`; do 
  echo -n "$x" >$tmp; 
  $WS_GLB_BIN/quart_bin fancy printQuartets $tmp;
done |sed -e "s/^.*: //"> /Users/Erfan/Documents/Reasearch/global/Erfan/Test/test.txt;
cat /Users/Erfan/Documents/Reasearch/global/Erfan/Test/test.txt | python $WS_GLB_PUTIL/summarize.quartets.py ; 

rm $tmp;
