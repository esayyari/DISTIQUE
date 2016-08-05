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

def readAnchoredFrqFromFile(out_path):
    frq = dict()
    f = open(out_path, 'r')
    frq = dict()
    for line in f:
        k=line.split()
        frq[k[0]] = [float(k[1]),float(k[2])]
    f.close()
    return frq
def readTaxaList(out_path):
    taxa_list = dict()
    f = open(out_path, 'r')
    for line in f:
        k=line.split()
        taxa_list[k[0]] = k[1:len(k)]
    f.close()
    return taxa_list

def printQuartTable(quartTable,out_path):
    out_path = os.path.expanduser(os.path.expandvars(out_path))
    orig_stdout = sys.stdout
    f = open(out_path,'w')
    sys.stdout = f
    for k in quartTable:
        print k,str(quartTable[k][0]),quartTable[k][1]
    sys.stdout = orig_stdout
    f.close()
    return
def printQuartetTableAveraged(quartTable,out_path):
    out_path = os.path.expanduser(os.path.expandvars(out_path))
    orig_stdout = sys.stdout
    f = open(out_path,'w')
    sys.stdout = f
    for k in quartTable:
        print k,quartTable[k]

    sys.stdout = orig_stdout
    f.close()
    return
def printTaxaList(taxaList,out_path):
    out_path = os.path.expanduser(os.path.expandvars(out_path))
    orig_stdout = sys.stdout
    f = open(out_path,'w')
    sys.stdout = f
    for k in taxaList:
        print k,
        for v in taxaList[k]:
            print v,
        print 
    f.close()
    sys.stdout = orig_stdout
    return
def readAnchoredQuartTableFromFile(out_path):
    frq = dict()
    f = open(out_path, 'r')
    frq = dict()
    for line in f:
        k=line.split()
        frq[k[0]] = float(k[1])
    f.close()
    return frq
def printDistanceTable(Dtmp,Ctmp,Dpath,Cpath):
    out_pathD = os.path.expanduser(os.path.expandvars(Dpath))
    out_pathC = os.path.expanduser(os.path.expandvars(Cpath))
    fD = open(out_pathD,'w')
    fC = open(out_pathC,'w')
    for key in Dtmp:
        print >> fD,key,Dtmp[key]
        print >> fC,key,Ctmp[key]
    fD.close()
    fC.close()
    return