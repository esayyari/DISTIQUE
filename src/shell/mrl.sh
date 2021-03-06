#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export FASTMRP=$WS_HOME/mrpmatrix
raxml="$WS_HOME/standard-RAxML/raxmlHPC"

outdir="$(cd "$1" && pwd)"
for x in `cat $outdir/distance-*.nwk`; do   echo $x | sed -e 's/\[&U\]//g'; done>$outdir/trees.nwk
for x in `find $outdir -name "distance-*.nwk"`; do y=$(cat $x | sed -e 's/\[&U\]//g'); echo -n $y; echo; done > $outdir/trees.nwk
in=$outdir/trees.nwk
out=anchored_mrl_tree.nwk
met=$2
if [ "$met" == "" ]; then
	met="mrl"
fi
#if [ -s $outdir/$out ]; then
 # echo "Ouput files exists. Refusing to rerun. "
  #exit 0;
#fi
cd $outdir
tmp=`mktemp $outdir/mrpmaptrix.$out.XXXXX`
if [ "$met" == "mrl" ]; then 
java -jar $FASTMRP/mrp.jar $in $tmp  PHYLIP -randomize

#rm -f RAxML_*$out*
$raxml -m BINCAT -s $tmp -n $out -N 2 -p $RANDOM > /dev/null 2>&1

#test "$?" != "0" && exit 1

mv RAxML_bestTree.$out* $out
mv RAxML_info.$out mrl.log

rm RAxML_*$out* 
else
	java -jar -Xmx3000M $WS_HOME/ASTRAL/astral.4.10.3.jar -i $in -o $out > /dev/null 2>&1
fi

echo MRL Done. Output at: $outdir/$out


