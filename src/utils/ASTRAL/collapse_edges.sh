#!/bin/bash
source $WS_HOME/global/src/shell/setup.sh

for x in `find $path -name ""`; do

$WS_GLB_PUTIL/remove_edges_from_tree.py    $x 75 "$x-col"

done
