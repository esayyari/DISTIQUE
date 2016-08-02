#!/bin/bash
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
set -x
show_help(){
cat << EOF
USAGE: ${0##*/} [-h] [-i DATA PATH] [-g QUARTET PATH] [-o OUTPATH] [-n NUMROUNDS] [-v VERSION]
Runs DISTIQUE Anchoring on all replicates available under DATA PATH, based on available quartet files. The results will be written at OUTPATH. NUMROUNDS is number of rounds for sampling around each polytomy, while VERSION indicates which version of the code to use.
EOF
}
n=2
v=3
while getopts "hi:g:o:n:v:s:" opt; do
        case $opt in
        h)
                show_help
                exit 0
                ;;
        i)
                i=$OPTARG
                ;;
        o)
                o=$OPTARG
                ;;
        g)
		q=$OPTARG
                ;;
	n)
		n=$OPTARG
		;;
	v)
		v=$OPTARG
		;;
	s)
		sp=$OPTARG
		;;
        '?')
                printf "Unknown input option"
                show_help
                ;;
        esac
done

if [ -z "$i" ]; then
        printf "Enter data path\n"
        show_help
        exit 1
fi
if [ -z "$o" ]; then
        printf "Enter output path\n"
        show_help
        exit 1
fi
if [ -z "$q" ]; then
      printf "Enter quartet file path\n"  
      show_help
      exit 1
fi

WS_LOC_UTIL=$WS_HOME/DISTIQUE/src/utils/
WS_LOC_SHEL=$WS_HOME/DISTIQUE/src/shell/
tmpDIR=`mktemp -d`
echo $tmpDIR

for r in `find $i -maxdepth 1 -type d -name "R*"`; do
    R=$(echo $r | sed -e 's/.*\///')
    tmptmpDIR=`mktemp -d $tmpDIR/tmp$R.XXXXXXX`
    echo $tmptmpDIR
      if [ "$v" == "0" ]; then
	 l=freq
	 a=mean
	 m=prod
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_BalME-$l-$a-"$m".XXXXX`
         /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z B -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
       	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_OLSME-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z O -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
         tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-BIONJ-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z I -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-NJ-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z N -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_BalME2-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z B2  -l $l -a $a -m $m> $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_OLSME2-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py-g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z O2 -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-PhyDstar-MVR-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u phydstar -z MVR -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-PhyDstar-BioNJ-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u phydstar -z BioNJ -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-cons-add-freq-prod-mean.XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -l freq -m prod -a mean -u fastme -z D> $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-cons-max-freq-min-mean.XXXXX`
       	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -l freq -m min  -a mean -u fastme -z D> $tmptmptmpDIR/results.log 2>&1
 
	 l=log
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_BalME-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z B -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_OLSME-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z O -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-BIONJ-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z I -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-NJ-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z N -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_BalME2-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z B2  -l $l -a $a -m $m> $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_OLSME2-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py-g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z O2 -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-PhyDstar-MVR-$l-$a-"$m".XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u phydstar -z MVR -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-PhyDstar-BioNJ-$l-$a-"$m".XXXXX`
         /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u phydstar -z BioNJ -l $l -a $a -m $m > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-cons-add-$l-prod-mean.XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -l $l -m $m -a $a -u fastme -z D > $tmptmptmpDIR/results.log 2>&1
	 tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-cons-max-$l-min-mean.XXXXX`
	 /usr/bin/time -p $WS_LOC_UTIL/distique-2.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -l $l -m min  -a $a -u fastme -z D> $tmptmptmpDIR/results.log 2>&1
    fi
done
for x in `find $tmpDIR -name "distance.d_distique_tree.nwk*"`; do
	k=$(dirname $x)
	$WS_LOC_SHEL/compare.tree.sh -s $sp -g $x > $k/results-score.sc
done
find $tmpDIR -name "distance.d_distique_tree.nwk*" > $tmpDIR/listToTar
find $tmpDIR -name "results-score.sc" >> $tmpDIR/listToTar
find $tmpDIR -name "results.log" >> $tmpDIR/listToTar

y=$(echo $o | sed -e 's/^.*\///')
tar czf $o/testAnchoring-V$v-$n.tar.gz -T $tmpDIR/listToTar 
rm -r $tmpDIR
