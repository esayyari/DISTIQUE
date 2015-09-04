#!/bin/bash

export DIR=`pwd`
source $DIR/setup.sh
sp=10
$WS_GLB_BIN/quart_bin genTree $sp > $WS_GLB_TS/test.nwk;
for i in {1..100}; do
	
	echo "">>$WS_GLB_TS/test.nwk;
	echo "`head -n 1 $WS_GLB_TS/test.nwk`">> $WS_GLB_TS/test.nwk;
done
for i in {1..400}; do
	echo "">>$WS_GLB_TS/test.nwk;
	$WS_GLB_BIN/quart_bin genTree $sp >> $WS_GLB_TS/test.nwk;
done

$DIR/quartet.freq.sh $WS_GLB_TS/test.nwk;
fastme -v 0 -i $WS_GLB_RES/quartetDistance.txt;
sed -ie 's/:[[:digit:]]*\+\.[[:digit:]]*\|:-[[:digit:]]*\+\.[[:digit:]]*//g' $WS_GLB_RES/quartetDistance.txt_fastme_tree.nwk;

echo "";
echo "`head -n 1 $WS_GLB_RES/quartetDistance.txt_fastme_tree.nwk`";

echo "";
echo "`head -n 1 $WS_GLB_TS/test.nwk`">$WS_GLB_TS/ref.nwk;

echo "`head -n 1 $WS_GLB_TS/test.nwk`";
echo "`tail -n 2 $WS_GLB_TS/test.nwk`";
echo "`head -n 1 $WS_GLB_RES/quartetDistance.txt_fastme_tree.nwk`"> $WS_GLB_RES/quartetDistance.txt_fastme_tree.nwk

$DIR/compareTrees.missingBranch $WS_GLB_TS/ref.nwk  $WS_GLB_RES/quartetDistance.txt_fastme_tree.nwk;
