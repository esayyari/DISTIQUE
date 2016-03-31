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
while getopts "hi:g:o:n:v:" opt; do
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
tmpDIR=`mktemp -d`
echo $tmpDIR

for r in `find $i -maxdepth 1 -type d -name "R*"`; do
	R=$(echo $r | sed -e 's/.*\///')
	tmptmpDIR=`mktemp -d $tmpDIR/tmp$R.XXXXXXX`
	echo $tmptmpDIR
#	for ((C=1;C<=5;C++)); do
#	{
		#tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp$C.XXXXXX`
	if [ "$v" == "0" ]; then
		tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_BalME.XXXXX`
        	/usr/bin/time -p $WS_LOC_UTIL/distique.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z B  > $tmptmptmpDIR/results.log 2>&1
       		tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_OLSME.XXXXX`
	        /usr/bin/time -p $WS_LOC_UTIL/distique.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z O  > $tmptmptmpDIR/results.log 2>&1
        	tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-BIONJ.XXXXX`
	        /usr/bin/time -p $WS_LOC_UTIL/distique.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z I  > $tmptmptmpDIR/results.log 2>&1
        	tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-NJ.XXXXX`
	        /usr/bin/time -p $WS_LOC_UTIL/distique.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z N  > $tmptmptmpDIR/results.log 2>&1
		tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_BalME2.XXXXX`
                /usr/bin/time -p $WS_LOC_UTIL/distique.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z B2  > $tmptmptmpDIR/results.log 2>&1
                tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_OLSME2.XXXXX`
                /usr/bin/time -p $WS_LOC_UTIL/distique.py-g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u fastme -z O2  > $tmptmptmpDIR/results.log 2>&1
                tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-PhyDstar-MVR.XXXXX`
                /usr/bin/time -p $WS_LOC_UTIL/distique.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u phydstar -z MVR  > $tmptmptmpDIR/results.log 2>&1
                tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-PhyDstar-BioNJ.XXXXX`
                /usr/bin/time -p $WS_LOC_UTIL/distique.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -u phydstar -z BioNJ  > $tmptmptmpDIR/results.log 2>&1
	fi
	tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_BalME.XXXXX`
	/usr/bin/time -p $WS_LOC_UTIL/testAnchoring-v$v.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -n $n -u fastme -z B  > $tmptmptmpDIR/results.log 2>&1
	tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_OLSME.XXXXX`
        /usr/bin/time -p $WS_LOC_UTIL/testAnchoring-v$v.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -n $n -u fastme -z O  > $tmptmptmpDIR/results.log 2>&1
	tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-BIONJ.XXXXX`
        /usr/bin/time -p $WS_LOC_UTIL/testAnchoring-v$v.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -n $n -u fastme -z I  > $tmptmptmpDIR/results.log 2>&1
	tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-NJ.XXXXX`
        /usr/bin/time -p $WS_LOC_UTIL/testAnchoring-v$v.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -n $n -u fastme -z N  > $tmptmptmpDIR/results.log 2>&1
	tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_BalME2.XXXXX`
        /usr/bin/time -p $WS_LOC_UTIL/testAnchoring-v$v.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -n $n -u fastme -z B2  > $tmptmptmpDIR/results.log 2>&1
        tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-fastme-TaxAdd_OLSME2.XXXXX`
        /usr/bin/time -p $WS_LOC_UTIL/testAnchoring-v$v.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -n $n -u fastme -z O2  > $tmptmptmpDIR/results.log 2>&1
        tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-PhyDstar-MVR.XXXXX`
        /usr/bin/time -p $WS_LOC_UTIL/testAnchoring-v$v.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -n $n -u phydstar -z MVR  > $tmptmptmpDIR/results.log 2>&1
        tmptmptmpDIR=`mktemp -d $tmptmpDIR/tmp-PhyDstar-BioNJ.XXXXX`
        /usr/bin/time -p $WS_LOC_UTIL/testAnchoring-v$v.py -g $r/genetrees.gt -o $tmptmptmpDIR -f $q/$R/quartets.q -n $n -u phydstar -z BioNJ  > $tmptmptmpDIR/results.log 2>&1
#	}
	done
#done

y=$(echo $o | sed -e 's/^.*\///')
tar czf $o/testAnchoring-V$v-$y-$n.tar.gz $tmpDIR/* 
rm -r $tmpDIR
