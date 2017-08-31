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


class analyze:
    def __init__(self,opt):
        self.opt = opt

    def distique(self):
        if self.opt.strat == 1:
            self.prodmethod()
        elif self.opt.strat == 2:
            self.minmethod()
        elif self.opt.strat == 3:
            self.distancesum()
        elif self.opt.strat == 4:
            self.treesum()

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
            frq = tbs.findQuartetTable(self.opt.trees, self.opt.taxa, 0, self.opt.outpath, self.opt.verbose)
            if self.opt.verbose:
                print "time to find quartet lists: "
                tm.toc()
        if self.opt.verbose:
            tm.tic()

        for e in self.opt.con_tree.postorder_node_iter():

            if e in self.opt.to_resolve:

                self.opt.to_resolve[e]
                (self.opt.taxa_list, taxa_inv) = tstt.getTaxaList(self.opt.to_resolve[e])

                if self.opt.verbose:
                    print "computing the partial quartet table"

                quartTable = tbs.findTrueAverageTable(frq, self.opt.taxa_list, self.opt.av, self.opt.met)

                if self.opt.verbose:
                    print "computing distance table using the method: " + str(method)

                ftmp3 = tempfile.mkstemp(suffix='.d', prefix="distancet.d", dir=self.opt.outpath, text=None)
                tbs.distanceTable(quartTable, method, ftmp3[1], self.opt.met)
                os.close(ftmp3[0])
                ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix="distance.d_fastme_tree.nwk", dir=self.opt.outpath, text=None)
                FNULL = open(os.devnull, 'w')
                subprocess.call(
                    [self.opt.WS_LOC_FM + "/fastme", "-i", ftmp3[1], "-w", "none", "-o", ftmp4[1], "-I", "/dev/null"],
                    stdout=FNULL, stderr=subprocess.STDOUT)
                os.close(ftmp4[0])

                if self.opt.verbose:
                    print "starting to resolve polytomy"

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
            frq = tbs.findQuartetTable(self.opt.trees, self.opt.taxa, 0, self.opt.outpath, self.opt.verbose)

            print "time to find quartet lists: "
            tm.toc()
        tm.tic()
        for e in self.opt.con_tree.postorder_node_iter():
            if e in self.opt.to_resolve:
                val = self.opt.to_resolve[e]
                (self.opt.taxa_list, taxa_inv) = tstt.getTaxaList(self.opt.to_resolve[e])
                if self.opt.verbose:
                    print "computing the partial quartet table"
                quartTable = tbs.findTrueAverageTable(frq, self.opt.taxa_list, self.opt.av, self.opt.met)
                if self.opt.verbose:
                    print "computing distance table using the method: " + str(method)
                ftmp3 = tempfile.mkstemp(suffix='.d', prefix="distancet.d", dir=self.opt.outpath, text=None)
                tbs.distanceTable(quartTable, method, ftmp3[1], self.opt.met)
                os.close(ftmp3[0])
                ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix="distance.d_fastme_tree.nwk", dir=self.opt.outpath, text=None)
                FNULL = open(os.devnull, 'w')
                subprocess.call(
                    [self.opt.WS_LOC_FM + "/fastme", "-i", ftmp3[1], "-w", "none", "-o", ftmp4[1], "-I", "/dev/null"],
                    stdout=FNULL, stderr=subprocess.STDOUT)
                os.close(ftmp4[0])

                if self.opt.verbose:
                    print "starting to resolve polytomy"
                res = tstt.resolvePolytomy(ftmp4[1], e, self.opt.verbose)
        print "resolving polytomies takes about: "
        tm.toc()

        outfile = self.opt.outpath + "/distique_all_pairs_max.nwk"
        if self.opt.verbose:
            print "writing the resulting tree as: " + outfile
        tstt.changeLabelsToNames(self.opt.con_tree, self.opt.new_labels, self.opt.verbose)
        self.opt.con_tree.write(path=outfile, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)

    def distancesum(self):
        if self.opt.verbose:
            tm.tic()

        skippedPoly = set()

        (to_resolvettt, _) = tstt.findPolytomies_with_names(self.opt.con_tree)
        (ac, acSmall, smallAnchs, numAnchors, n) = self.pickAnchors(to_resolvettt)

        if self.opt.verbose:
            self.printInfo(n, skippedPoly, numAnchors, smallAnchs)

        count_distance_table = dict()
        distance_tables = []
        listPoly = self.opt.to_resolve.keys()
        for z in range(len(listPoly)):
            distance_tables.append([])
        self.count_distance_table = count_distance_table
        self.distance_tables = distance_tables
        if len(ac) != 0:
            self.computeAllFinalDistanceTablesBigPoly(listPoly, skippedPoly, ac)

        self.computeAllFinalDistanceTablesSmallPoly(listPoly,skippedPoly,acSmall)

        normalizedD = self.distance_tables
        normalizedC = self.count_distance_table
        for z in range(len(listPoly)):
            atbs.normalizeDistanceTable(normalizedD[z], normalizedC[z])

        self.writeAnchorsTofile( ac, acSmall)

        for e in self.opt.con_tree.postorder_node_iter():
            if e in self.opt.to_resolve:
                z = listPoly.index(e)
                keyDict = sorted(list(np.unique((" ".join(normalizedD[z].keys())).split(" "))))
                fileDistance = "distancet-" + str(e.label) + ".d"
                ftmp3 = tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=self.opt.outpath, text=None)
                pr.printDistanceTableToFile(normalizedD[z], keyDict, ftmp3[1])
                print "writing distance table to " + str(ftmp3[1])
                os.close(ftmp3[0])
                ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix=fileDistance + "_fastme_tree.nwk", dir=self.opt.outpath,
                                         text=None)
                tstt.buildTreeFromDistanceMatrix(ftmp3[1], ftmp4[1], self.opt.sumProg, self.opt.sumProgOption)
                os.close(ftmp4[0])
                if self.opt.verbose:
                    print "starting to resolve polytomy"
                    atbs.resolvePolytomy(ftmp4[1], e, self.opt.verbose)
        tstt.prune_tree_trivial_nodes(self.opt.con_tree)
        outfile = self.opt.outpath + "/distique_distance-sum.nwk"
        if self.opt.verbose:
            print "writing the resulting tree as: " + outfile
        tstt.changeLabelsToNames(self.opt.con_tree, self.opt.new_labels, self.opt.verbose)
        self.opt.con_tree.write(path=outfile, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)
        print "The overall time to infer the species tree is: "
        tm.toc()



    def treesum(self):
        skippedPoly = set()

        (to_resolvettt, _) = tstt.findPolytomies_with_names(self.opt.con_tree)
        (ac, acSmall, smallAnchs, numAnchors, n) = self.pickAnchors(to_resolvettt)
        numSmallAnchors = 0
        if len(ac) == 0:
            for e in acSmall:
                numAnchors += len(acSmall[e][0]) * 1
        else:
            for e in acSmall:
                numSmallAnchors += len(acSmall[e][0])
            numAnchors = len(ac)

        if self.opt.verbose:
            print "Number of taxa is: " + str(n)
            print "Number of polytomies is: " + str(len(self.opt.to_resolve))
            for k in self.opt.to_resolve:
                print "The size of polytomy around node " + k.label + " is: " + str(len(self.opt.to_resolve[k].keys()))
                if len(self.opt.to_resolve[k].keys()) < 6:
                    skippedPoly.add(k)

            print "Averaging distances around polytomies using method: " + self.opt.met
            print "The number of anchores are: " + str(numAnchors)
        TreeList = dict()
        if ac is not None:
            count = 1
            for anch in ac:
                if self.opt.verbose:
                    tm.tic()
                anch = sorted(list(anch))
                print anch
                for e in self.opt.to_resolve:
                    if e in skippedPoly:
                        continue
                    val = self.opt.to_resolve[e]
                    (taxa_list, taxa_inv) = tstt.getTaxaList(val)
                    N1 = taxa_inv[anch[0]]
                    N2 = taxa_inv[anch[1]]
                    if N1 == N2:
                        continue
                    N = {N1, N2}
                    if self.opt.verbose:
                        print "The size of polytomy is: " + str(len(taxa_list))
                    tm.tic()

                    if self.opt.readFromFile:
                        frq = atbs.findAnchoredDistanceTableFromFile(anch, self.opt.frqT, self.opt.taxa, self.opt.outpath)
                    else:
                        frq = atbs.findAnchoredDistanceTableOverallp(e, N, anch, taxa_list, taxa_inv, self.opt.trees, self.opt.taxa,
                                                                     self.opt.outpath,
                                                                     self.opt.debugFlag)

                    if self.opt.verbose:
                        print "computing the partial quartet table"
                    if self.opt.readFromFile:
                        quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesOverallFromFile(frq, anch,
                                                                                                       taxa_list,
                                                                                                       N1, N2, self.opt.am, self.opt.met)
                    else:
                        quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesOverall(frq, anch, taxa_list, N1,
                                                                                               N2,
                                                                                        self.opt.am, self.opt.met)
                    if e.label not in TreeList:
                        TreeList[e.label] = dendropy.TreeList()
                    if self.opt.verbose:
                        print "computing distance table using the method: " + str(self.opt.am)
                    D = atbs.anchoredDistanceFromFrq(quartTable, anch)

                    keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
                    fileDistance = "distancet-" + str(anch[0]) + "-" + str(anch[1]) + ".d"
                    ftmp3 = tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=self.opt.outpath, text=None)
                    pr.printDistanceTableToFile(D, keyDict, ftmp3[1])
                    os.close(ftmp3[0])
                    ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix=fileDistance + "_fastme_tree.nwk", dir=self.opt.outpath,
                                             text=None)
                    tstt.buildTreeFromDistanceMatrix(ftmp3[1], ftmp4[1], self.opt.sumProg, self.opt.sumProgOption)
                    os.close(ftmp4[0])
                    tree_tmp = dendropy.Tree.get(path=ftmp4[1], schema='newick', rooting="force-unrooted")
                    TreeList[e.label].append(tree_tmp)
                if self.opt.verbose:
                    tm.toc()
                    print "The anchor " + str(count) + " out of " + str(len(ac)) + " anchors has been finished!"
                count += 1

        if self.opt.verbose:
            print "Start finding resolution for polytomies with degree smaller than 6"
        count = 1
        for e in skippedPoly:
            i = 0
            val = self.opt.to_resolve[e]
            (taxa_list, taxa_inv) = tstt.getTaxaList(val)
            for achList in acSmall[e.label]:
                quartTable = dict()
                if self.opt.verbose:
                    tm.tic()
                for anch in achList:
                    anch = sorted(list(anch))
                    if self.opt.verbose:
                        print anch
                        print "time elapsing  for counting number of quartets for this anchors is: "
                        tm.toc()
                    N1 = taxa_inv[anch[0]]
                    N2 = taxa_inv[anch[1]]
                    if N1 == N2:
                        continue
                    N = {N1, N2}
                    if self.opt.verbose:
                        print "The size of polytomy is: " + str(len(taxa_list))
                    tm.tic()

                    if self.opt.readFromFile:
                        frq = atbs.findAnchoredDistanceTableFromFile(anch, self.opt.frqT, self.opt.taxa, self.opt.outpath)
                    else:
                        frq = atbs.findAnchoredDistanceTableOverallp(e, N, anch, taxa_list, taxa_inv, self.opt.trees, self.opt.taxa,
                                                                     self.opt.outpath,
                                                                     self.opt.debugFlag)
                    if e.label not in TreeList:
                        TreeList[e.label] = dendropy.TreeList()


                    if self.opt.verbose:
                        print "computing the partial quartet table"
                    if self.opt.readFromFile:
                        quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomiesOverallFromFile(frq,
                                                                                                         anch,
                                                                                                         taxa_list,
                                                                                                         N1, N2)
                    else:
                        quartTable = tbsa.findTrueAverageTableAnchoringOnDifferentSidesSmallPolytomiesOverall(frq, anch,
                                                                                                 taxa_list, N1, N2,
                                                                                                 self.opt.met)
                    if self.opt.verbose:
                        print "The anchor " + str(count) + " out of " + str(
                            numSmallAnchors) + " anchors has been finished!"
                    count += 1
                if self.opt.verbose:
                    print "computing distance table using the method: " + str(self.opt.am)
                Frq = atbs.anchoredDistanceFromFrqSmallPolytomies(quartTable, self.opt.am, self.opt.met)
                D = pd.prodDistance(Frq, self.opt.met)
                keyDict = sorted(list(np.unique((" ".join(D.keys())).split(" "))))
                fileDistance = "distancet-anchList-" + str(i) + ".d"
                ftmp3 = tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=self.opt.outpath, text=None)
                pr.printDistanceTableToFile(D, keyDict, ftmp3[1])
                os.close(ftmp3[0])
                ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix=fileDistance + "_fastme_tree.nwk", dir=self.opt.outpath,
                                         text=None)
                tstt.buildTreeFromDistanceMatrix(ftmp3[1], ftmp4[1], self.opt.sumProg, self.opt.sumProgOption)
                os.close(ftmp4[0])
                tree_tmp = dendropy.Tree.get(path=ftmp4[1], schema='newick')
                if e.label in TreeList:
                    TreeList[e.label].append(tree_tmp)
                else:
                    TreeList[e.label] = dendropy.TreeList()
                    TreeList[e.label].append(tree_tmp)

                if self.opt.verbose:
                    tm.toc()
        if self.opt.removeOutliers and self.opt.verbose:
            print "removeing outliers"

        for e in TreeList:
            if self.opt.removeOutliers:
                tstt.remove_outliers(TreeList[e], self.opt.strategy, self.opt.outpath, e, self.opt.summary)
        for e in self.opt.con_tree.postorder_node_iter():
            if e in self.opt.to_resolve:
                ftmp4 = tstt.findMRL(TreeList[e.label], e.label, self.opt.outpath, self.opt.summary)
                if self.opt.verbose:
                    print "starting to resolve polytomy"
                    atbs.resolvePolytomy(ftmp4, e, self.opt.verbose)
        tstt.prune_tree_trivial_nodes(self.opt.con_tree)

        outfile = self.opt.outpath + "/distique_tree-sum.nwk"
        if self.opt.verbose:
            print "writing the resulting tree as: " + outfile
        tstt.changeLabelsToNames(self.opt.con_tree, self.opt.new_labels, self.opt.verbose)
        self.opt.con_tree.write(path=outfile, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)



    def resolveAllPolytomies(self,listPoly, normalizedD):
        for e in self.opt.con_tree.postorder_node_iter():
            if e in self.opt.to_resolve:
                z = listPoly.index(e)
                keyDict = sorted(list(np.unique((" ".join(normalizedD[z].keys())).split(" "))))
                fileDistance = "distancet-" + str(e.label) + ".d"
                ftmp3 = tempfile.mkstemp(suffix='.d', prefix=fileDistance, dir=self.opt.outpath, text=None)
                pr.printDistanceTableToFile(normalizedD[z], keyDict, ftmp3[1])
                print "writing distance table to " + str(ftmp3[1])
                os.close(ftmp3[0])
                ftmp4 = tempfile.mkstemp(suffix='.nwk', prefix=fileDistance + "_fastme_tree.nwk", dir=self.opt.outpath,
                                         text=None)
                tstt.buildTreeFromDistanceMatrix(ftmp3[1], ftmp4[1], self.opt.sumProg, self.opt.sumProg)
                os.close(ftmp4[0])
                if self.opt.verbose:
                    print "starting to resolve polytomy"
                    atbs.resolvePolytomy(ftmp4[1], e, self.opt.verbose)

    def writeAnchorsTofile(self, ac, acSmall):
        fileAnch = "listAnchors"
        ftmp3 = tempfile.mkstemp(suffix='.txt', prefix=fileAnch, dir=self.opt.outpath, text=None)
        f = open(ftmp3[1], 'w')
        for anch in ac:
            anch = sorted(list(anch))
            anch_temp = "|".join(anch)
            print >> f, anch_temp

        for achList in acSmall:
            for anch in achList:
                anch = sorted(list(anch))
                anch_temp = "|".join(anch)
                print >> f, anch_temp

        f.close()
        os.close(ftmp3[0])

    def printInfo(self, n, skippedPoly, numAnchors, smallAnchs):
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
    def computeFrqandQuartTables(self,anch,e,N,taxa_list,taxa_inv):
        if self.opt.readFromFile:
            frq = atbs.findAnchoredDistanceTableFromFile(anch, self.opt.frqT, self.opt.taxa, self.opt.outpath)
        else:
            frq = atbs.findAnchoredDistanceTableOverallp(e, N, anch, taxa_list, taxa_inv, self.opt.trees, self.opt.taxa,
                                                         self.opt.outpath,
                                                         self.opt.debugFlag)

        print "time elapsing time for counting number of quartets for anch " + anch[0] + " " + anch[1]
        tm.toc()

        tm.tic()
        if self.opt.readFromFile:
            quartTable = tbsa.findTrueAverageTableAnchoringAddDistancesOverallFromFile(frq, anch, taxa_list,
                                                                                       N,
                                                                                       self.opt.method, self.opt.met)
        else:
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
        numAnchors = 0
        if self.opt.randomSample:
            ac = tstt.pickAnchors(self.opt.taxa, self.opt.to_resolve, self.opt.num, self.opt.debugFlag)
        smallAnchs = set()
        for e in to_resolvettt:
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

    def computeDistancesForOneAnchor(self,anch,taxa_list,taxa_inv, z,e):
        anch = sorted(list(anch))

        N1 = taxa_inv[anch[0]]
        N2 = taxa_inv[anch[1]]
        if N1 == N2:
            contFlag = False
            return(None, None, contFlag)
        N = {N1, N2}
        if self.opt.verbose:
            print "The size of polytomy is: " + str(len(taxa_list))
        (quartTable, frq) = self.computeFrqandQuartTables(anch, e, N, taxa_list, taxa_inv)
        [Dtmp, Ctmp] = atbs.anchoredDistanceFromFrqAddDistances(quartTable, anch, taxa_list)
        atbs.fillEmptyElementsDistanceTable(Dtmp, Ctmp, self.opt.fillmethod)
        contFlag = True
        self.addDistances(z, Dtmp, Ctmp)
        del quartTable
        del frq
        return (contFlag)


    def computeDistancesForAllAnchors(self,ac,taxa_list,taxa_inv,z,e,count):
        for anch in ac:
            (contFlag) = self.computeDistancesForOneAnchor(anch, taxa_list,taxa_inv,z,e)
            if not contFlag:
                continue

            if self.opt.verbose:
                print "Computing distance table using anchors " + anch[0] + " and " + anch[
                    1] + " has been finished!"
                print "The anchor " + str(count) + " out of " + str(len(ac)) + " anchors has been finished!"
            count += 1
            if (count % 50) == 0:
                gc.collect()

        return (count)



    def computeAllFinalDistanceTablesBigPoly(self,listPoly,skippedPoly,ac):
        count = 1
        for z in range(len(listPoly)):
            e = listPoly[z]
            if e in skippedPoly:
                continue
            val = self.opt.to_resolve[e]
            (taxa_list, taxa_inv) = tstt.getTaxaList(val)
            (count) = self.computeDistancesForAllAnchors(ac, taxa_list,taxa_inv,z,e, count)
        return

    def computeAllFinalDistanceTablesSmallPoly(self,listPoly,skippedPoly,acSmall):
        countT = 1
        skippedPolyList = list(skippedPoly)
        for e in skippedPolyList:
            z = listPoly.index(e)
            val = self.opt.to_resolve[e]
            (taxa_list, taxa_inv) = tstt.getTaxaList(val)
            for achList in acSmall[e.label]:
                (countT) = self.computeDistancesForAllAnchors(achList,taxa_list,taxa_inv,z, e,countT)
        return