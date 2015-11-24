#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export FASTMRP=~/Documents/Research/mrpmatrix
export raxml="raxmlHPC"
outdir=$1
for x in `cat $outdir/distance-*.nwk`; do   echo $x | sed -e 's/\[&U\]//g'; done>$outdir/trees.nwk
for x in `find $outdir -name "distance-*.nwk"`; do y=$(cat $x | sed -e 's/\[&U\]//g'); echo -n $y; echo; done > $outdir/trees.nwk
in=$outdir/trees.nwk
out=anchored_fastme_tree.nwk


cd $outdir

if [ -s $out ]; then
  echo "Ouput files exists. Refusing to rerun. "
  exit 0;
fi

tmp=`mktemp -p . mrpmaptrix.$out.XXXXX`
java -jar $FASTMRP/mrp.jar $in $tmp  PHYLIP -randomize

rm -f RAxML_*$out*

$raxml -m BINCAT -s $tmp -n $out -N 2 -p $RANDOM

test "$?" != "0" && exit 1

mv RAxML_bestTree.$out* $out
mv RAxML_info.$out mrl.log

rm RAxML_*$out* $tmp

echo MRL Done. Output at: $out


