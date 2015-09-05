#!/bin/bash

export DIR=`pwd`
source $DIR/setup.sh
source $WS_HOME/global/src/shell/setup.sh
sp=6
$WS_LOC_BIN/quart_bin genTree $sp > $WS_LOC_TS/test.nwk;
for i in {1..100}; do
	
	echo "">>$WS_LOC_TS/test.nwk;
	echo "`head -n 1 $WS_LOC_TS/test.nwk`">> $WS_LOC_TS/test.nwk;
done
#for i in {1..400}; do
#	echo "">>$WS_LOC_TS/test.nwk;
#	$WS_LOC_BIN/quart_bin genTree $sp >> $WS_LOC_TS/test.nwk;
#done

$DIR/quartet.freq.sh $WS_LOC_TS/test.nwk;
fastme -v 0 -i $WS_LOC_TS/test.nwk_d
sed -ie 's/:[[:digit:]]*\+\.[[:digit:]]*\|:-[[:digit:]]*\+\.[[:digit:]]*//g' $WS_LOC_TS/test.nwk_d_fastme_tree.nwk

echo "";
echo "`head -n 1 $WS_LOC_TS/test.nwk_d_fastme_tree.nwk`";

echo "";
echo "`head -n 1 $WS_LOC_TS/test.nwk`">$WS_LOC_TS/ref.nwk;

echo "`head -n 1 $WS_LOC_TS/test.nwk`";
echo "`tail -n 2 $WS_LOC_TS/test.nwk`";
echo "`head -n 1 $WS_LOC_TS/test.nwk_d_fastme_tree.nwk`"> $WS_LOC_TS/test.nwk_d_fastme.tree.nwk

$WS_GLB_SH/compareTrees.missingBranch $WS_LOC_TS/ref.nwk  $WS_LOC_TS/test.nwk_d_fastme.tree.nwk
