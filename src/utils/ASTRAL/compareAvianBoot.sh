#!/bin/bash

test $# == 2 || { echo USAGE: $0 bootstrap_genetrees species_tree; exit 1; }
DIR=$WS_HOME/global/src/shell
cdir=`pwd`;
source $DIR/setup.sh
outFold=$(dirname "$1")
mkdir -p $cdir/$outFold
out=$cdir/$outFold
echo $out
x=$(basename "$1"| sed -e 's/\..*//')
tmpFile=`mktemp -p $out $x.XXXXX`;
outFile=`mktemp -p $out $x-comp-sp.XXXXX`;
gunzip -c $1 > $tmpFile
$WS_GLB_SH/compareTrees.missingBranch $2 $tmpFile > $outFile;
trueTmpFile=$(basename "$tmpFile" | sed -e 's/\..*$//' )
trueTmpFolder=$(dirname "$tmpFile")
mv $tmpFile $trueTmpFolder/$trueTmpFile

trueOutFile=$(basename "$outFile" | sed -e 's/\..*$//' )
trueOutFolder=$(dirname "$outFile")
mv $outFile $trueOutFolder/$trueOutFile

