#!/usr/bin/env python
import dendropy
import sys
import os
import numpy as np
import itertools
import subprocess
import printTools as pr
import anchoredTableTools as atbs
import toolsTreeTaxa as tstt
import timer as tm
from optparse import OptionParser
import tableManipulationToolsAnchoring as tbsa
import tempfile
import warnings
import testDependencies as tD
usage = "usage: %prog [options]"
parser = OptionParser(usage)

parser.add_option("-f",dest="freq",
        help="The path to frequency table path")
parser.add_option("-a",dest="anch",
                  help = "anchors")
parser.add_option("-t",dest="list",
                  help = "list of taxa around polytomy")
parser.add_option("-m",dest="met",
        help="The method to summerize quartet results around each node, freq, or log, Default is log", default="log")
parser.add_option("-l",dest="fillmethod",
        help="The method to fill empty cells in distance tables, const, rand, or normConst. Default is const", default="const")
parser.add_option("-o",dest="outp",
        help="The output folder to write the output file")
(options,args) = parser.parse_args()
filename = options.freq
anch = list((options.anch).split(","))
taxa_list_f = options.list
method = options.fillmethod
met = options.met
op = options.outp

src_fpath = os.path.expanduser(os.path.expandvars(filename))
if not os.path.exists(src_fpath):
    sys.stderr.write('Not found: "%s"' % src_fpath)
frq = tD.readAnchoredFrqFromFile(src_fpath)

src_fpath = os.path.expanduser(os.path.expandvars(taxa_list_f))
if not os.path.exists(src_fpath):
    sys.stderr.write('Not found: "%s"' % src_fpath)
taxa_list = tD.readTaxaList(src_fpath)
method = "mean"
quartTable = tbsa.findTrueAverageTableAnchoringAddDistances(frq,anch,taxa_list,method,met)
out_path = op
out_path = out_path+"/quartTable_test.qt"

tD.printQuartetTableAveraged(quartTable,out_path)

