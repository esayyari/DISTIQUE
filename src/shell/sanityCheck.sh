#!/bin/bash

export DIR=`pwd`
source $DIR/setup.sh
source $WS_HOME/global/src/shell/setup.sh
sp=11
$WS_GLB_BIN/quart_bin genTree $sp > $WS_LOC_TS/testSanity/test.nwk;
for i in {1..100}; do
	echo "">>$WS_LOC_TS/testSanity/test.nwk;
	echo "`head -n 1 $WS_LOC_TS/testSanity/test.nwk`">> $WS_LOC_TS/testSanity/test.nwk;
done
for i in {1..900}; do
	echo "">>$WS_LOC_TS/testSanity/test.nwk;
	$WS_GLB_BIN/quart_bin genTree $sp >> $WS_LOC_TS/testSanity/test.nwk;
done

$DIR/summarize.tree.sh -w 1 -s fm -d prod -r $WS_LOC_TS/testSanity  $WS_LOC_TS/testSanity/test.nwk;
echo "";
echo "`head -n 1 $WS_LOC_TS/testSanity/distance.d_fastme_tree.nwk`";

echo "";
echo "`head -n 1 $WS_LOC_TS/testSanity/test.nwk`">$WS_LOC_TS/testSanity/ref.nwk;

echo "`head -n 1 $WS_LOC_TS/testSanity/test.nwk`";
echo "`tail -n 2 $WS_LOC_TS/testSanity/test.nwk`";
echo "`head -n 1 $WS_LOC_TS/testSanity/distance.d_fastme_tree.nwk`"> $WS_LOC_TS/testSanity/distance.d_fastme_tree.nwk

$WS_GLB_SH/compareTrees.missingBranch $WS_LOC_TS/testSanity/ref.nwk  $WS_LOC_TS/testSanity/distance.d_fastme_tree.nwk
