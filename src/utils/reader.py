import os
import toolsTreeTaxa as tstt
import dendropy
import tableManipulationToolsAnchoring as tbsa
import timer as tm
import copy
import tempfile
import sys
class asllOptions:
    def __init__(self,options):
        return self.reader(options)

    def reader(self,options):
        self.strat = options.strat

        self.filename = options.filename
        self.gt = options.gt
        self.outpath = options.out
        self.thr = options.thr
        self.av = options.av
        self.sumProg = options.sumProg
        self.sumProgOption = options.sumProgOption
        self.verbose = options.verbose
        self.debugFlag = (options.debug == "1")
        self.sp = options.sp
        self.num = int(options.num)
        self.summary = options.summary
        self.met = options.met
        self.strategy = options.outlier
        self.fillmethod = options.fillmethod
        self.initTree = options.initTree
        self.annotation = options.annotation
        if self.verbose:
            print self.strategy
        if self.strategy is not None:
            self.removeOutliers = True
        else:
            self.removeOutliers = False

        self.am = options.am
        self.randomSample = True
        if options.filename:
            self.readFromFile = True
            self.frqT = tbsa.readTable(self.filename)
        else:
            self.readFromFile = False

        if self.sumProg == "phydstar" and self.sumProgOption == "":
            self.sumProgOption = "BioNJ"
        elif self.sumProg == "phydstar":
            self.sumProgOption = options.sumProgOption

        if options.filename:
            self.readFromFile = True
        else:
            self.readFromFile = False
        if (not options.gt or not options.out):
            sys.exit("Please enter genetrees file, and output folder location")

        self.src_fpath = os.path.expanduser(os.path.expandvars(self.gt))

        self.trees = dendropy.TreeList.get_from_path(self.src_fpath, 'newick')

        (self.converted_labels, self.new_labels) = tstt.changeLabelsToNumbers(self.trees, self.verbose)

        if self.initTree is not None:
            if self.verbose:
                print "reading the initial tree with polytomies"
            self.con_tree = dendropy.Tree.get(path=self.initTree,schema="newick")
        else:
            if self.verbose:
                print "computing the consensus tree with the threshold " + str(self.thr)
                print "time to compute consensus is: "
                tm.tic()
            self.con_tree = self.trees.consensus(min_freq=self.thr)
            if self.verbose:
                tm.toc()

        self.ftmpt = tempfile.mkstemp(suffix='.nwk', prefix="consensusTree", dir=self.outpath, text=None)
        self.con_tree2 = copy.deepcopy(self.con_tree)
        tstt.changeLabelsToNames(self.con_tree2, self.new_labels, self.verbose)

        self.con_tree2.write(path=self.ftmpt[1], schema="newick", suppress_rooting=True)

        os.close(self.ftmpt[0])
        tstt.labelNodes(self.con_tree)

        (self.to_resolve, self.maxPolyOrder) = tstt.findPolytomies(self.con_tree)
        self.taxa = list()
        if self.annotation is not None:
            self.multiInd = True
        else:
            self.multiInd = False
        if self.multiInd:
            allLines = open(self.annotation,'r').readlines()
            self.mapping = list()
        idx = 0
        self.mapSpeciesToIdx = dict()
        for e in self.con_tree.leaf_nodes():
            self.taxa.append(e.taxon.label)
            self.mapSpeciesToIdx[e.taxon.label] = idx
            idx += 1
        if self.multiInd:
            for idx in range(0,len(self.con_tree.leaf_nodes())):
                self.mapping.append([])
        if self.multiInd:
            for line in allLines:
                line = line.strip('\n')
                line = line.strip()
                listLine = line.split('\t')
                species = listLine[1]
                spIdx = self.mapSpeciesToIdx[species]
                self.mapping[spIdx].append(listLine[0])


        self.n = len(self.con_tree.leaf_nodes())
        self.num_taxa = len(self.trees[0].leaf_nodes())
        self.method = "prod"
        self.WS_LOC_SHELL = os.environ['WS_HOME'] + '/DISTIQUE/src/shell'
        self.WS_LOC_FM = os.environ['WS_HOME'] + '/DISTIQUE/bin'
        self.mapSpeciesdict = dict()








