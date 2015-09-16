#!/bin/bash
test $# == 2 || { echo USAGE: $0 path_to_species_trees dictionary_file; exit 1; }
DIR=$( dirname "${BASH_SOURCE[0]}" )
source $DIR/setup.sh
for data_name in `find $1 -name "s_tree.trees"`; do
path=$(dirname "${data_name}")
printf "$path\n"
printf "$data_name\n"
#python $WS_LOC_PUTIL/mapsequences.py $data_name $2 $path/true_sp_tree.nwk
done
