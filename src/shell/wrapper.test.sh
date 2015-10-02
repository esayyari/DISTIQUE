#!/bin/bash
set -e
DIR=$( pwd )
#b=( 50 100 200 )
b=( 10 25 50 100 200 )
g=( 100 500 1000 )
res=~/Documents/Research/results/
dataroot=~/Documents/Research/data/11-taxon/
model=( model.10.1800000.0.000000111  model.10.200000.0.000001000  model.10.5400000.0.000000037  model.10.600000.0.000000333 )
#model=(  model.10.5400000.0.000000037 ) 
#model=model.10.600000.0.000000333
out=~/Documents/Research/results/tests/min
for mov in "${model[@]}"; do
        for bv in "${b[@]}"; do
                for gv in "${g[@]}"; do
                         for i in {1..50}; do
                                mkdir -p $res/11-taxon/$mov/$i/relabeled_shortened_data_$bv"_"subset_$gv/$1/$i
                                resulttmp=$res/11-taxon/$mov/$i/relabeled_shortened_data_$bv"_"subset_$gv/$1/$i
                                $DIR/summarize.tree.sh -s fm -d $1 -p 4 -r $resulttmp  -w 1 $dataroot/$mov/$i/relabeled_shortened_data_$bv"_"subset_$gv/all_fasttree.genetrees;
				$DIR/compare.tree.sh -s $dataroot/$mov/$i/S_relabeled_tree.trees -g $resulttmp/distance.d_fastme_tree.nwk > $out/11-taxon_$mov"-"11-fm_$1"-"$bv"-"$gv"-"$i"-"results.txt
                                 printf "$i gene is $gv base is $bv model is $mov\n $resulttmp\n $dataroot/$mov/$i/relabeled_shortened_data_$bv"_"subset_$gv/all_fasttree.genetrees\n";
                        done
                done
        done
done

#grep -r "" | sed -e 's/\.11\./ 11 /' | sed -e 's/fm_prod\./fm_prod /' | sed -e 's/\./-/g' | sed -e 's/-/\./' | sed -e 's/-/\./' | sed -e 's/-/\./' | sed -e 's/-/\./' |  sed -e 's/-/ /' | sed -e 's/-/ /' | sed -e 's/-/ /' |sed -e 's/results-txt:[0-9]* [0-9]* //g' | sed -e 's/-/ /'  |sed -e 's/--/ -/' | sed -e 's/-/\./'| sed -e 's/ /\./' |sed -e 's/ /,/g'| sed -e 's/11.taxon/11-taxon/' > newres.csv
