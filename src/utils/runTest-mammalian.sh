#!/bin/bash

DIR=$WS_HOME/DISTIQUE/src/utils
path=$WS_HOME/oasis/data/mammalian/
for x in `ls $path`; do
        for z in `ls $path/$x`; do
                mkdir -p ~/Documents/Research/oasis/outputs/new_out/mammalian-$x-$z
                $DIR/distique.py -g $path/$x/$z/genetrees.gt -o ~/Documents/Research/oasis/outputs/new_out/mammalian-$x-$z -f ~/Documents/Research/oasis/outputs/mammalian/$x/$z/distique-2/prod/prod/quartet_tmp.q -l log;
                echo "working on mammalian-$x-$z has been finished!"
        done
done

