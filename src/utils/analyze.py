import os

import numpy as np
import tableManipulationTools as tbs
import tempfile
import timer as tm
import subprocess
import printTools as pr
import anchoredTableTools as atbs
import toolsTreeTaxa as tstt
import tableManipulationToolsAnchoring as tbsa
import dendropy
import gc
import prodDistance as pd
import random
random.seed(a=121089923)

class analyze:
    def __init__(self,opt):
        self.opt = opt

    def distique(self):
        if self.opt.strat == 1:
            self.prodmethod()
        elif self.opt.strat == 2:
            self.minmethod()
        elif self.opt.strat == 3:
            if self.opt.verbose:
                tm.tic()
            self.distancesum()
            if self.opt.verbose:
                print "Species Tree infered after :"
                tm.toc()

        elif self.opt.strat == 4:
            if self.opt.verbose:
                tm.tic()
            self.treesum()
            if self.opt.verbose:
                print "Species Tree infered after :"
                tm.toc()
    def prodmethod(self):
        method = "prod"

        if self.opt.verbose:
            print "The summary program is: " + self.opt.sumProg
            print "The self.option for this summary program is: " + self.opt.sumProg
            print "Number of taxa is: " + str(self.opt.n)
            print "the number of polytomies is: " + str(len(self.opt.to_resolve))
            print "the maximum order of polytomies is: " + str(self.opt.maxPolyOrder)
            print "computing the total quartet table"
        if self.opt.readFromFile:
            frq = tbs.readTable(self.opt.filename)
        else:
            if self.opt.verbose:
                tm.tic()
            frq = tbs.findQuartetTable(self.opt.trees, self.opt.taxa, 0, self.opt.tempFold, self.opt.verbose)
            if self.opt.verbose:
                print "time to find quartet lists: "
                tm.toc()
        if self.opt.verbose:
            tm.tic()

        for e in self.opt.con_tree.postorder_node_iter():

            if e in self.opt.to_resolve:

                self.opt.to_resolve[e]
                (self.opt.taxa_list, taxa_inv) = tstt.getTaxaList(self.opt.to_resolve[e])



                quartTable = tbs.findTrueAverageTable(frq, self.opt.taxa_list, self.opt.av, self.opt.met)

                if self.opt.verbose:
                    print "computing distance table using the method: " + str(method)

                ftmp3 = tempfile.mkstemp(suffix='.d', prefix="distancet.d", dir=self.opt.tempFold, text=None)
                tbs.distanceTable(quartTable, method, ftmp3[1], self.opt.met)
                os.close(ftmp3[0])
                ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix="distance.d_fastme_tree.nwk", dir=self.opt.tempFold, text=None)
                FNULL = open(os.devnull, 'w')
                subprocess.call(
                    [self.opt.WS_LOC_FM + "/fastme", "-i", ftmp3[1], "-w", "none", "-o", ftmp4[1], "-I", "/dev/null"],
                    stdout=FNULL, stderr=subprocess.STDOUT)
                os.close(ftmp4[0])


                tstt.resolvePolytomy(ftmp4[1], e, self.opt.verbose)

        if self.opt.verbose:
            print "resolving polytomies takes about: "
            tm.toc()

        outfile = self.opt.outpath + "/distique_all_pairs_prod.nwk"

        if self.opt.verbose:
            print "writing the resulting tree as: " + outfile

        tstt.changeLabelsToNames(self.opt.con_tree, self.opt.new_labels, self.opt.verbose)
        self.opt.con_tree.write(path=self.opt.outfile, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)

    def minmethod(self):
        method = "min"

        if self.opt.verbose:
            print "computing the total quartet table"
        if self.opt.readFromFile:
            frq = tbs.readTable(self.opt.filename)
        else:
            tm.tic()
            frq = tbs.findQuartetTable(self.opt.trees, self.opt.taxa, 0, self.opt.tempFold, self.opt.verbose)

            print "time to find quartet lists: "
            tm.toc()
        if self.opt.verbose:
            tm.tic()
        for e in self.opt.con_tree.postorder_node_iter():
            if e in self.opt.to_resolve:
                val = self.opt.to_resolve[e]
                (self.opt.taxa_list, taxa_inv) = tstt.getTaxaList(self.opt.to_resolve[e])

                quartTable = tbs.findTrueAverageTable(frq, self.opt.taxa_list, self.opt.av, self.opt.met)
                if self.opt.verbose:
                    print "computing distance table using the method: " + str(method)
                ftmp3 = tempfile.mkstemp(suffix='.d', prefix="distancet.d", dir=self.opt.tempFold, text=None)
                tbs.distanceTable(quartTable, method, ftmp3[1], self.opt.met)
                os.close(ftmp3[0])
                ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix="distance.d_fastme_tree.nwk", dir=self.opt.tempFold, text=None)
                FNULL = open(os.devnull, 'w')
                subprocess.call(
                    [self.opt.WS_LOC_FM + "/fastme", "-i", ftmp3[1], "-w", "none", "-o", ftmp4[1], "-I", "/dev/null"],
                    stdout=FNULL, stderr=subprocess.STDOUT)
                os.close(ftmp4[0])


                res = tstt.resolvePolytomy(ftmp4[1], e, self.opt.verbose)
        if self.opt.verbose:
            print "resolving polytomies takes about: "
            tm.toc()

        outfile = self.opt.outpath + "/distique_all_pairs_max.nwk"
        if self.opt.verbose:
            print "writing the resulting tree as: " + outfile
        tstt.changeLabelsToNames(self.opt.con_tree, self.opt.new_labels, self.opt.verbose)
        self.opt.con_tree.write(path=outfile, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)

    def distancesum(self):
        (to_resolvettt, _) = tstt.findPolytomies_with_names(self.opt.con_tree)
        (ac, acSmall, smallAnchs, numAnchors, n) = self.pickAnchors(to_resolvettt)
        skippedPoly = self.printInfo(n, numAnchors, smallAnchs)
        listPoly = self.opt.to_resolve.keys()

        self.initializeTables(listPoly)
        if len(ac) != 0:
            self.computeAllFinalDistanceTablesBigPoly(listPoly, skippedPoly, ac)
        self.computeAllFinalDistanceTablesSmallPoly(listPoly,skippedPoly,acSmall,smallAnchs)

        (normalizedD,normalizedC) = self.normalizeTables(listPoly)

        self.writeAnchorsTofile( ac, acSmall)
        self.resolveAllPolytomies(listPoly, normalizedD)
        tstt.prune_tree_trivial_nodes(self.opt.con_tree)

        outfile = self.opt.outpath + "/distique_distance-sum.nwk"

        if self.opt.verbose:
            print "writing the resulting tree as: " + outfile

        tstt.changeLabelsToNames(self.opt.con_tree, self.opt.new_labels, self.opt.verbose)
        self.opt.con_tree.write(path=outfile, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)



    def treesum(self):


        (to_resolvettt, _) = tstt.findPolytomies_with_names(self.opt.con_tree)
        (ac, acSmall, smallAnchs, numAnchors, n) = self.pickAnchors(to_resolvettt)
        if self.opt.verbose:
            tm.tic()
        skippedPoly = self.printInfo(n,numAnchors,smallAnchs)
        self.writeAnchorsTofile( ac, acSmall)
        self.TreeList = dict()

        if ac is not None:
            self.updateTreeListBigPolyAllPolytomies(ac, skippedPoly)
        if self.opt.verbose:
            print "Start finding resolution for polytomies with degree smaller than 6"
        self.updateTreeListSmallAllPolytomies(skippedPoly,acSmall)
        if self.opt.removeOutliers:
            if self.opt.verbose:
                print "removeing outliers"
            for e in self.TreeList:
                tstt.remove_outliers(self.TreeList[e], self.opt.strategy, self.opt.tempFold, e, self.opt.summary)

        for e in self.opt.con_tree.postorder_node_iter():
            if e in self.opt.to_resolve:
                ftmp4 = tstt.findMRL(self.TreeList[e.label], e.label, self.opt.tempFold, self.opt.summary)

                atbs.resolvePolytomy(ftmp4, e, self.opt.verbose)
        tstt.prune_tree_trivial_nodes(self.opt.con_tree)

        outfile = self.opt.outpath + "/distique_tree-sum.nwk"
        if self.opt.verbose:
            print "writing the resulting tree as: " + outfile
        tstt.changeLabelsToNames(self.opt.con_tree, self.opt.new_labels, self.opt.verbose)
        self.opt.con_tree.write(path=outfile, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)
        if self.opt.verbose:
            print "The overll time to estimate tree is: "
            tm.toc()



    def resolveAllPolytomies(self,listPoly, normalizedD):
        for e in self.opt.con_tree.postorder_node_iter():
            if e in self.opt.to_resolve:
                z = listPoly.index(e)
                keyDict = sorted(list(np.unique((" ".join(normalizedD[z].keys())).split(" "))))
                fileDistance = "distancet-" + str(e.label) + ".d"
                ftmp3 = tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=self.opt.tempFold, text=None)
                pr.printDistanceTableToFile(normalizedD[z], keyDict, ftmp3[1])
                print "writing distance table to " + str(ftmp3[1])
                os.close(ftmp3[0])
                ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix=fileDistance + "_fastme_tree.nwk", dir=self.opt.tempFold,
                                         text=None)
                tstt.buildTreeFromDistanceMatrix(ftmp3[1], ftmp4[1], self.opt.sumProg, self.opt.sumProgOption)
                os.close(ftmp4[0])
                atbs.resolvePolytomy(ftmp4[1], e, self.opt.verbose)

    def writeAnchorsTofile(self, ac, acSmall):
        fileAnch = "listAnchors"
        ftmp3 = tempfile.mkstemp(suffix='.txt', prefix=fileAnch, dir=self.opt.tempFold, text=None)
        f = open(ftmp3[1], 'w')
        for anch in ac:
            anch = sorted(list(anch))
            anch_temp = "|".join(anch)
            print >> f, anch_temp
        for key in acSmall.keys():
            achList = acSmall[key]
            for ac in achList:
                for anch in ac:
                    anch = sorted(list(anch))
                    anch_temp = "|".join(anch)
                    print >> f, anch_temp

        f.close()
        os.close(ftmp3[0])

    def printInfo(self, n, numAnchors, smallAnchs):
            skippedPoly = set()
            print "Number of taxa is: " + str(n)
            print "Number of polytomies is: " + str(len(self.opt.to_resolve))
            for k in self.opt.to_resolve:
                print "The size of polytomy around node " + k.label + " is: " + str(len(self.opt.to_resolve[k].keys()))
                if len(self.opt.to_resolve[k].keys()) < 6:
                    skippedPoly.add(k)
            print "Filling method is: " + self.opt.fillmethod
            print "Distance method is: " + self.opt.method
            print "Averaging distances around polytomies using method: " + self.opt.met
            print "The number of anchors are: " + str(numAnchors)
            print "The number of anchors for small polytomies is: " + str(len(smallAnchs))
            return skippedPoly
    def computeFrqandQuartTables(self,anch,e,N,taxa_list,taxa_inv,mults):
        if self.opt.readFromFile:
            frq = atbs.findAnchoredDistanceTableFromFile(anch, self.opt.frqT, self.opt.taxa, self.opt.tempFold,self.opt.mapSpeciesToIdx,self.opt.mapping)
        else:
            anch1 = taxa_inv[anch[0]]
            anch2 = taxa_inv[anch[1]]
            A = sorted([anch1, anch2])
            frq = atbs.findAnchoredDistanceTableOverallp(e, N, anch, taxa_list, taxa_inv, self.opt.trees, self.opt.taxa,
                                                         self.opt.tempFold,
                                                         self.opt.debugFlag,self.opt.mapSpeciesToIdx,self.opt.mapping,A,self.opt.geneIndTable,mults)

        if self.opt.verbose:
            print "time elapsing time for counting number of quartets for anch " + anch[0] + " " + anch[1]

        if self.opt.readFromFile:
            quartTable = tbsa.findTrueAverageTableAnchoringAddDistancesOverallFromFile(frq, anch, taxa_list,
                                                                                       N,
                                                                                       self.opt.method, self.opt.met)
        else:
            anch1 = taxa_inv[anch[0]]
            anch2 = taxa_inv[anch[1]]
            anch = sorted([anch1,anch2])
            quartTable = tbsa.findTrueAverageTableAnchoringAddDistancesOverall(frq, anch, taxa_list, N,
                                                                               self.opt.method,
                                                                               self.opt.met)
        return (quartTable, frq)

    def addDistances(self,z,Dtmp,Ctmp):
        if len(self.distance_tables[z]) == 0:
            if self.opt.fillmethod == "normConst":
                Dmax = max(np.abs(Dtmp.values())) * 1.
                if Dmax == 0:
                    Dmax = 1.
                for kttDtmp in Dtmp:
                    Dtmp[kttDtmp] = Dtmp[kttDtmp] / Dmax
                self.distance_tables[z] = Dtmp
                self.count_distance_table[z] = Ctmp
            else:
                self.distance_tables[z] = Dtmp
                self.count_distance_table[z] = Ctmp
        else:
            atbs.addDistanceAnchores(self.distance_tables[z], Dtmp, self.count_distance_table[z], Ctmp, self.opt.fillmethod)
        return

    def pickAnchors(self,to_resolvettt):

        acSmall = dict()
        if self.opt.randomSample:
            ac = tstt.pickAnchors(self.opt.taxa, self.opt.to_resolve, self.opt.num, self.opt.debugFlag,self.opt.seedNum)
        numAnchors = len(ac)
        smallAnchs = set()
        for e in sorted(to_resolvettt.keys()):
            if len(to_resolvettt[e].keys()) < 6:
                val = to_resolvettt[e]
                (taxa_list, taxa_inv) = tstt.getTaxaList(val)
                acSmall[e] = tstt.chooseAnchoresAll(taxa_list, 1, self.opt.debugFlag)
                for acList in acSmall[e]:
                    for anch in acList:
                        smallAnchs.add(anch)
                        if len(ac) == 0:
                            continue
                        ac.append(anch)

        n = len(self.opt.con_tree.leaf_nodes())
        ac = set(ac)
        if len(ac) == 0:
            for e in acSmall:
                numAnchors += len(acSmall[e][0]) * 1
        else:
            numAnchors = len(ac)
        return (ac, acSmall, smallAnchs, numAnchors,n)

    def computeDistancesForOneAnchor(self,anch,taxa_list,taxa_inv, z,e,mults):
        anch = sorted(list(anch))

        N1 = taxa_inv[anch[0]]
        N2 = taxa_inv[anch[1]]
        if N1 == N2:
            contFlag = False
            return(None, None, contFlag)
        N = {N1, N2}
        (quartTable, frq) = self.computeFrqandQuartTables(anch, e, N, taxa_list, taxa_inv,mults)
        [Dtmp, Ctmp] = atbs.anchoredDistanceFromFrqAddDistances(quartTable, sorted(list(N)), taxa_list)
        atbs.fillEmptyElementsDistanceTable(Dtmp, Ctmp, self.opt.fillmethod)
        contFlag = True
        self.addDistances(z, Dtmp, Ctmp)
        del quartTable
        del frq
        return (contFlag)


    def computeDistancesForAllAnchors(self,ac,taxa_list,taxa_inv,z,e,count,mults):
        for anch in ac:
            A = anch
            (contFlag) = self.computeDistancesForOneAnchor(A, taxa_list,taxa_inv,z,e,mults)
            if not contFlag:
                continue

            if self.opt.verbose:
                print "Computing distance table using anchors " + anch[0] + " and " + anch[
                 1] + " has been finished!"
                print "The anchor " + str(count) + " out of " + str(len(ac)) + " anchors has been finished!"
            count += 1
            if (count % 100) == 0:
                gc.collect()

        return (count)



    def computeAllFinalDistanceTablesBigPoly(self,listPoly,skippedPoly,ac):
        count = 1
        for z in range(len(listPoly)):
            e = listPoly[z]
            if e in skippedPoly:
                continue
            val = self.opt.to_resolve[e]
            (taxa_list, taxa_inv) = tstt.getTaxaList2(val,self.opt.mapping,self.opt.mapSpeciesToIdx)
            mults = self.preprocessEmptyQuartTables(len(self.opt.trees), taxa_list)

            (count) = self.computeDistancesForAllAnchors(ac, taxa_list,taxa_inv,z,e, count, mults)
        return

    def computeAllFinalDistanceTablesSmallPoly(self,listPoly,skippedPoly,acSmall,smallAnchs):
        countR = 1
        skippedPolyList = list(skippedPoly)
        for e in skippedPolyList:
            z = listPoly.index(e)
            val = self.opt.to_resolve[e]
            (taxa_list, taxa_inv) = tstt.getTaxaList2(val,self.opt.mapping,self.opt.mapSpeciesToIdx)
            countT = 1
            mults = self.preprocessEmptyQuartTables(len(self.opt.trees), taxa_list)


            for achList in acSmall[e.label]:

                (countT) = self.computeDistancesForAllAnchors(achList,taxa_list,taxa_inv,z, e,countT,mults)
            print "Number of small polytomies finished is " + str(countR) + " out of " + str(len(acSmall.keys()))
            countR += 1
        return

    def initializeTables(self,listPoly):
        count_distance_table = dict()

        distance_tables = []

        for z in range(len(listPoly)):
            distance_tables.append([])
        self.count_distance_table = count_distance_table
        self.distance_tables = distance_tables

    def normalizeTables(self,listPoly):
        normalizedD = self.distance_tables
        normalizedC = self.count_distance_table

        for z in range(len(listPoly)):
            atbs.normalizeDistanceTable(normalizedD[z], normalizedC[z])
        return (normalizedD,normalizedC)

    def computeFrqandQuartTablesTreesum(self,anch,e,N1,N2,taxa_list,taxa_inv):
        N = {N1,N2}
        if self.opt.readFromFile:
            frq = atbs.findAnchoredDistanceTableFromFile(anch, self.opt.frqT, self.opt.taxa, self.opt.tempFold,self.opt.mapSpeciesToIdx,self.opt.mapping)

        else:
            frq = atbs.findAnchoredDistanceTableOverallp(e, N, anch, taxa_list, taxa_inv, self.opt.trees, self.opt.taxa,
                                                         self.opt.tempFold,
                                                         self.opt.debugFlag,self.opt.mapSpeciesToIdx,self.opt.mapping)

        if self.opt.readFromFile:
            quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesOverallFromFile(frq, anch,
                                                                                           taxa_list,
                                                                                           N1, N2, self.opt.am,
                                                                                           self.opt.met)
        else:
            quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesOverall(frq, anch, taxa_list, N1,
                                                                                   N2,
                                                                                   self.opt.am, self.opt.met)
        return (quartTable,frq)

    def computeFrqandQuartTablesTreesumSkippedPoly(self, anch, e, N1, N2, taxa_list, taxa_inv, quartTable):

        N = {N1, N2}
        if self.opt.readFromFile:
            frq = atbs.findAnchoredDistanceTableFromFile(anch, self.opt.frqT, self.opt.taxa,self.opt.tempFold,self.opt.mapSpeciesToIdx,self.opt.mapping)

        else:

            frq = atbs.findAnchoredDistanceTableOverallp(e, N, anch, taxa_list, taxa_inv, self.opt.trees, self.opt.taxa,
                                                         self.opt.tempFold,
                                                         self.opt.debugFlag,self.opt.mapSpeciesToIdx,self.opt.mapping)

        if self.opt.readFromFile:
            quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomiesOverallFromFile(frq, quartTable, anch, taxa_list,
                                                                                        N1, N2)

        else:
            quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomiesOverall(frq, quartTable, anch, taxa_list, N1, N2,self.opt.met)

        return (quartTable,frq)

    def makeTreeInTreesum(self,quartTable,anch,e):
        D = atbs.anchoredDistanceFromFrq(quartTable, anch)
        keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
        fileDistance = "distancet-" + str(anch[0]) + "-" + str(anch[1]) + ".d"
        ftmp3 = tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=self.opt.tempFold, text=None)
        pr.printDistanceTableToFile(D, keyDict, ftmp3[1])
        os.close(ftmp3[0])
        ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix=fileDistance + "_fastme_tree.nwk", dir=self.opt.tempFold,
                                 text=None)
        tstt.buildTreeFromDistanceMatrix(ftmp3[1], ftmp4[1], self.opt.sumProg, self.opt.sumProgOption)
        os.close(ftmp4[0])
        tree_tmp = dendropy.Tree.get(path=ftmp4[1], schema='newick', rooting="force-unrooted")
        self.TreeList[e.label].append(tree_tmp)

    def makeTreeInTreesumSkipped(self,quartTable,i,e):
        Frq = atbs.anchoredDistanceFromFrqSmallPolytomies(quartTable, self.opt.am, self.opt.met)
        D = pd.prodDistance(Frq, self.opt.met)
        keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
        fileDistance = "distancet-anchList-" + str(i) + ".d"
        ftmp3 = tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=self.opt.tempFold, text=None)
        pr.printDistanceTableToFile(D, keyDict, ftmp3[1])
        os.close(ftmp3[0])
        ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix=fileDistance + "_fastme_tree.nwk", dir=self.opt.tempFold,
                                 text=None)
        tstt.buildTreeFromDistanceMatrix(ftmp3[1], ftmp4[1], self.opt.sumProg, self.opt.sumProgOption)
        os.close(ftmp4[0])
        tree_tmp = dendropy.Tree.get(path=ftmp4[1], schema='newick')
        self.TreeList[e.label].append(tree_tmp)

    def updateTreeListSmallPolyForOnePolytomy(self,acSmall,e,count):
        quartTable = dict()
        if e.label not in self.TreeList:
            self.TreeList[e.label] = dendropy.TreeList()

        i = 0
        val = self.opt.to_resolve[e]
        (taxa_list, taxa_inv) = tstt.getTaxaList(val)
        for achList in acSmall[e.label]:

            if self.opt.verbose:
                tm.tic()
            for anch in achList:
                anch = sorted(list(anch))
                N1 = taxa_inv[anch[0]]
                N2 = taxa_inv[anch[1]]
                if N1 == N2:
                    continue
                (quartTableTmp, frq) = self.computeFrqandQuartTablesTreesumSkippedPoly(anch, e, N1, N2, taxa_list,
                                                                                       taxa_inv, quartTable)
                quartTable = quartTableTmp
                count += 1
            self.makeTreeInTreesumSkipped(quartTable, i, e)
            i += 1
            if self.opt.verbose:
                tm.toc()

    def updateTreeListBigPolyForOnePolytomy(self,e,skippedPoly,anch):
        contFlag = False

        if e.label not in self.TreeList:
            self.TreeList[e.label] = dendropy.TreeList()
        if e in skippedPoly:
            contFlag = True
            return contFlag
        val = self.opt.to_resolve[e]
        (taxa_list, taxa_inv) = tstt.getTaxaList(val)
        N1 = taxa_inv[anch[0]]
        N2 = taxa_inv[anch[1]]
        if N1 == N2:
            contFlag = True
            return contFlag
        (quartTable, frq) = self.computeFrqandQuartTablesTreesum(anch, e, N1, N2, taxa_list, taxa_inv)


        self.makeTreeInTreesum(quartTable, anch, e)
        return contFlag

    def updateTreeListBigPolyAllPolytomies(self,ac,skippedPoly):
        count = 1
        for anch in ac:

            anch = sorted(list(anch))
            for e in self.opt.to_resolve:
                contFlag = self.updateTreeListBigPolyForOnePolytomy(e, skippedPoly, anch)
                if contFlag:
                    continue

            count += 1

    def updateTreeListSmallAllPolytomies(self,skippedPoly,acSmall):
        count = 1
        for e in skippedPoly:
            self.updateTreeListSmallPolyForOnePolytomy(acSmall, e, count)

    def getSpfromInd(self, indSp):
        return self.opt.mapIndToSp[indSp]

    def getIndlistFromSp(self, species):
        spIdx = self.opt.mapSpeciesToIdx[species]
        return self.mapping[spIdx]

    def addOne(self, indSp):
        self.counts[self.getSpfromInd(indSp)] += 1
        return

    def taxaToIndlist(self,taxa_list ):
        return

    def preprocessEmptyQuartTables(self,k,taxa_list):
        mults = dict()

        for gIdx in range(0, k):
            for key in taxa_list.keys():
                if key not in mults:
                    mults[key] = 0
                for sp in taxa_list[key]:
                    for l in range(0, len(self.opt.geneIndTable[gIdx][self.opt.mapSpeciesToIdx[sp]])):
                        mults[key] += self.opt.geneIndTable[gIdx][self.opt.mapSpeciesToIdx[sp]][l].mult
        return mults