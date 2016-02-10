#!/bin/bash

DIR=$WS_HOME/DISTIQUE/src/utils
path=$WS_HOME/oasis/data/11-taxon/
for x in `ls $path`; do
        for z in `ls $path/$x`; do
                for y in `ls $path/$x/$z`; do
                mkdir -p ~/Documents/Research/oasis/outputs/new_out/11-taxon-$x-$z-$y-corrected
                $DIR/distique.py -g $path/$x/$z/$y/genetrees.gt -o ~/Documents/Research/oasis/outputs/new_out/11-taxon-$x-$z-$y-corrected -f ~/Documents/Research/oasis/outputs/11-taxon/$x/$z/$y/distique-2/prod/prod/quartet_tmp.q -l log;
                echo "working on 11-taxon-$x-$z-$y has been finished!"
        done
done
done
