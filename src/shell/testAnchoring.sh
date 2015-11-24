#!/bin/bash

rm res.txt; rm test/*; for i in {1..20}; do ./testAnchoring.py -g ~/Documents/Research/oasis/data/mammalian/data-mammalian/0.2X-3200-500/R$i/genetrees.gt -o test -n 12 -t 0.5; ../shell/mrl.sh ~/Documents/Research/DISTIQUE/src/utils/test/; ../shell/compare.tree.sh -s ~/Documents/Research/oasis/data/mammalian/mammalian-model-species.tre -g test/anchored_fastme_tree.nwk >> res.txt; rm test/*; done
