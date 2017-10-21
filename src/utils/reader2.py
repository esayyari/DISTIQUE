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

        self.trees = dendropy.TreeList.get_from_path(self.src_fpath, 'newick', preserve_underscores=True)
        self.taxa = list()
        if self.annotation is not None:
            self.multiInd = True
        else:
            self.multiInd = False
        self.mapping = list()
        self.mapSpeciesToIdx = dict()
        self.mapIdxTospecies = list()
        self.mapIndToSp = dict()
        self.method = "prod"
        self.WS_LOC_SHELL = os.environ['WS_HOME'] + '/DISTIQUE/src/shell'
        self.WS_LOC_FM = os.environ['WS_HOME'] + '/DISTIQUE/bin'
        self.num_taxa = len(self.trees[0].leaf_nodes())
        self._mapIndToSpTrueLabels = dict()
        self._indTaxa = list
        self.geneIndTable = list()
        self._listIndForSpecies = list()
        self.reader()


    def reader(self):
        (self.converted_labels, self.new_labels2) = tstt.changeLabelsToNumbers(self.trees, self.verbose)
        self.makeBackBoneTree()
        (self.to_resolve, self.maxPolyOrder) = tstt.findPolytomies(self.con_tree)
        self.makeMappings()
        self.simplifyGenes()


    def mapNames(self):
        print self.converted_labels
        for node in self.con_tree.leaf_node_iter():
            node.taxon.label = self.converted_labels[str(node.taxon.label)]
        return

    def changeLabelsToNumbers(self,tree):
        converted_labels = dict()
        new_labels = dict()
        i = 0
        for taxon in tree.taxon_namespace:
            if taxon.label is not None:
                converted_labels[str(taxon.label)] = "l"+str(i)
                new_labels["l"+str(i)] = str(taxon.label)
                i += 1
        for node in tree.leaf_node_iter():
            node.taxon.label = converted_labels[str(node.taxon.label)]
        return (converted_labels,new_labels)



    def simplifyGenes(self):
        self._trees = copy.deepcopy(self.trees)

        for tree in self._trees:
            for leaf in tree.leaf_node_iter():
                leaf.spidx = self._mapIndToSpTrueLabels[leaf.taxon.label]

        for _, tree in enumerate(self._trees):
            tree.reroot_at_edge(tree.leaf_nodes()[0].edge, update_bipartitions=False)
            for node in tree.postorder_node_iter():
                node.desc_paths = set()
                if node.is_leaf():
                    node.desc_paths.add(node.spidx)
                else:
                    children = node.child_nodes()
                    flag = self.checkIfAllLeaf(children)

                    for child in children:
                        node.desc_paths = node.desc_paths.union(child.desc_paths)
                    if len(node.desc_paths) == 1:
                        c = len(children)
                        taxon = children[0].taxon.label
                        for chIdx in range(len(children)):
                            node.remove_child(children[chIdx])
                            node.taxon = dendropy.Taxon(label=taxon)
                            node.label = node.taxon.label
                            node.spidx = next(iter(node.desc_paths))
                            node.mult = c
                    elif flag:
                        h = dict()
                        h2 = set()
                        for child in children:
                            if child.spidx not in h:
                                h[child.spidx] = 0
                                h2.add(child)
                            else:
                                h[child.spidx] += 1
                        for child in children:
                            if child in h2:
                                child.mult = h[child.spidx]
                                continue
                            else:
                                node.remove_child(child, suppress_unifurcations=True)

        self.makeGeneTable()

    def checkIfAllLeaf(self,children):
        for child in children:
            if not child.is_leaf():
                return False
        return True

    def findAllNodes(self, tree, label):
        filter_fn = lambda n: n.taxon is not None and n.taxon.labe == label
        nodes = tree.find_nodes(filter_fn=filter_fn)
        return nodes

    def makeMappings(self):
        idx = 0
        if self.multiInd:
            self.readMappingFile()
        for e in self.con_tree.leaf_nodes():
            self.mapSpeciesToIdx[e.taxon.label] = idx
            self.mapIdxTospecies.append(e.taxon.label)
            idx += 1

        for idx in range(0, len(self.con_tree.leaf_nodes())):
            self.mapping.append([])

        if self.multiInd:
            for species in self.mapSpeciesToIdx:
                spIdx = self.mapSpeciesToIdx[species]
                self.mapping[spIdx].append(self.converted_labels[ind] for ind in self._listIndForSpecies[spIdx])
                self.mapIndToSp[self.converted_labels[listLine[0]]] = spIdx
                self.taxa.append(self.converted_labels[listLine[0]])

        else:
            for tx in self.con_tree.leaf_nodes():
                self.taxa.append(tx.taxon.label)
                species = tx.taxon.label
                spIdx = self.mapSpeciesToIdx[species]
                self.mapIndToSp[self.converted_labels[species]] = spIdx
                self.mapping[spIdx].append(species)

    def makeBackBoneTree(self):
        if self.initTree is not None:
            if self.verbose:
                print "reading the initial tree with polytomies"
            self.con_tree = dendropy.Tree.get(path=self.initTree,schema="newick",preserve_underscores=True)
        else:
            if self.verbose:
                print "computing the consensus tree with the threshold " + str(self.thr)
                print "time to compute consensus is: "
                tm.tic()
            self.con_tree = self.trees.consensus(min_freq=self.thr)
            if self.verbose:
                tm.toc()

        self.con_tree.write(path=self.outpath + "/consensus_tree.trees", schema="newick", suppress_rooting=True)

        (self.species_converted_labels, self.new_labels) = self.changeLabelsToNumbers(self.con_tree)
        tstt.labelNodes(self.con_tree)
        self.n = len(self.con_tree.leaf_nodes())


    def makeGeneTable(self):
        for gIdx in range(len(self._trees)):
            self.geneIndTable.append(list())
            for _ in range(len(self.new_labels)):
                self.geneIndTable[gIdx].append(list())

        for gIdx in range(len(self._trees)):
            for leaf in self._trees[gIdx].leaf_nodes():
                self.geneIndTable[gIdx][leaf.spidx].append(leaf)

    def readMappingFile(self):
        if self.multiInd:
            allLines = open(self.annotation,'r').readlines()
        for _ in range(len(self.mapSpeciesToIdx())):
            self._listIndForSpecies.append(list())

        if self.multiInd:
            for line in allLines:
                line = line.strip('\n')
                line = line.strip()
                listLine = line.split()
                species = listLine[1]
                spIdx = self.mapSpeciesToIdx[self.species_converted_labels[species]]
                self._mapIndToSpTrueLabels[listLine[0]] = spIdx
                self._indTaxa.append(listLine[0])
                self._listIndForSpecies[spIdx].append(listLine[0])
