#!/bin/bash

tar -xzf ./repository.tar.gz
shell_path=./repository/DISTIQ/src/shell

mv ./distiq.sh $shell_path
cd $shell_path

DIR=$shell_path
source $DIR/setup.sh
s=pd
e=nj
c=8
d=prod
o=1
true_sp=$2
gt=$1
outputdir=$3
outputfile=$4


mkdir -p $outputdir
./summarize.tree.sh -w 1 -s $s -e $e  -d $d -o $o -r $outputdir $gt

if [ "$s" -eq "pd" ]; then
$DIR/compare.tree.sh -s $true_sp -g $outputdir/distance.d_"$e".t> $outputdir/res.txt

elif [ "$s" -eq "fm" ]; then
$DIR/compare.tree.sh -s $true_sp -g $outputdir/distance.d_fastme_tree.nwk> $outputdir/res.txt
fi

tar czf $outputfile $outputdir
