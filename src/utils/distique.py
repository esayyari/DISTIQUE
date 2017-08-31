#!/usr/bin/env python
import os
from optparse import OptionParser
from analyze import analyze
from reader import asllOptions


usage = "usage: %prog [options]"
parser = OptionParser(usage)


parser.add_option("-a", "--averagemethod", dest="av",
                  help="The average method to find the average quartet table. Options are geometric mean (gmean), averaging (mean), or root mean square (otherwise). Default is mean.",
                  default="mean")
parser.add_option("-c", "--annotation", dest="annotation",
                 help="The annotion file for the multi-individual inference.", default=None)
parser.add_option("-d", dest="debug",
                  help="The flag for indicating that this run is for debugging!", default=False)
parser.add_option("-e", "--method", dest="am",
                  help="The averaging method for finding average quartet table", default="mean")
parser.add_option("-f", "--file", dest="filename", type="string",
                  help="read quartet table from FILENAME")
parser.add_option("-g", "--gene", dest="gt", type="string",
                  help="read genetrees from FILENAME")
parser.add_option("-i", "--initialtree", dest="initTree", type = str,
                  help="The initial backbone species tree with polytomies to be resolved instead of the consensus", default=None)
parser.add_option("-l", dest="fillmethod",
                  help="The method to fill empty cells in distance tables, const, rand, or normConst. Default is const",
                  default="const")
parser.add_option("-m", "--distmethod", dest="method", type=str,
                  help="The method to compute the distance of taxa. Options are prod or min. The default is prod.",
                  default="prod")
parser.add_option("-n", "--numStep", dest="num",
                  help="The number of anchors, default is 3", default=2)
parser.add_option("-o", "--output", dest="out", type="string",
                  help="the PATH to write the generated files")
parser.add_option("-p", dest="met",
                  help="The method to summerize quartet results around each node, freq, or log, Default is freq",
                  default="freq")
parser.add_option("-r", dest="outlier",
                  help="The strategy for outlier removal. The options are pairwise1, pairwise2, consensus10, or consensus3. Default is None",
                  default="consensus3")
parser.add_option("-s", "--sp", dest="sp",
                  help="species tree")
parser.add_option("-t", "--strategy", dest="strat", type="int",
                  help="The version of DISTIQUE to be run 1 (all-paris, prod), 2 (all-paris, max), 3 (Distance-sum), and 4 (Tree-sum), default is DISTIQUE Distance-SUM (3)",
                  default=3)
parser.add_option("-u", dest="sumProg",
                  help="The summerize method program to find species tree from distance matrix. The options are ninja, fastme, phydstar. Default is fastme ",
                  default="fastme")
parser.add_option("-v", "--verbose", dest="verbose",
                  help="Verbose", default=1)
parser.add_option("-x", dest="summary",
                  help="The summary method that will be used to summarize inferred species trees. Default is mrl",
                  default="mrl")
parser.add_option("-y", "--threshold", dest="thr", type=float,
                  help="the minimum frequency that consensus will use. Default is 0.5", default=0.5)
parser.add_option("-z", dest="sumProgOption",
                  help="The distance method to build the tree. If sumProg is set to fastme the options are TaxAdd_(B)alME (-s) (Default), TaxAdd_(B2)alME (-n), (D) default of fastme, TaxAdd_(O)LSME (-s), TaxAdd_(O2)LSME (-n), B(I)ONJ, (N)J. The default in this case is TaxAdd_(B)alME. if the  sumProg is set to phydstar, the options are BioNJ, MVR, and NJ. The default is TaxAdd_(B)alME.",
                  default="B")

(options, args) = parser.parse_args()
opt = asllOptions(options)
analyzer = analyze(opt)
analyzer.distique()

