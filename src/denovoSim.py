

# runs the denovo sim.

# the whole enchilada

#################
# DENOVO-PORTION#
#################

# simulate data

   # output expression in time steps, file per time step

# run denovo gene sets

   # outputs a list of trees.

# run denovo extractor

   # output a filtered expr for each tree passing a filter

# run statistical modeller.

   # output model results... predicted output, and summary



import numpy, time
import simGroups as ss
import filterSignal as fs
import denovoGeneSets as dg
import treeFilterAndExtration as cf
import model as mm
import analysis as an
import extractSubGraphs as es
import setScoring as scr


def runDenovoSim(datadir, Nf):
    # defaults
    ngenes = 120
    nparts = 6
    nsamples = 40
    filteredPrefix = "filtered_"
    denovoPrefix = "denovo_trees"
    levelThresh = 3
    topNTrees = 20
    crossVal = 5
    deltad = 6.0

    print('\nworking in ' + datadir)

    numpy.random.seed(seed=int(time.time()))

    # first simulate the data
    x = ss.runSim_DisjointSets(datadir, ngenes=ngenes, nparts=nparts, nsamples=nsamples, deltad=deltad)

    # filter the data
    y = fs.filterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])

    s = es.allSubgraphs(x[0],x[1],30,200)

    # recover the trees
    z = dg.denovoGeneSets(filelist=y[0], dirs=x[0], outputprefix=denovoPrefix, adjmat=x[1])

    # filter trees and extract data for modeling
    trees, genes, means = cf.treeFilterAndEx(dirs=x[1], treefile=z[0], filterfiles=y[0], levelThresh=levelThresh, topNTrees=topNTrees)

    # run models
    m = mm.rfModel(dirs=x[0], exprfile=x[2], pheno=x[3], genes=genes, cvs=crossVal)

    # compare model results to simulation.
    an.analysis(predacc=m, genes=genes, trees=trees, means=means, dirs=x[0], setfile=x[4])



def runDenovoSimReuseData(datadir, Nf):
    #
    # here let's not rebuild the network and data.
    #
    # defaults
    ngenes = 120
    nparts = 6
    nsamples = 10
    filteredPrefix = "filtered_"
    denovoPrefix = "denovo_trees"
    levelThresh = 3
    topNTrees = 20
    crossVal = 5
    deltad = 6.0
    Nf = int(Nf)

    print('\nworking in ' + datadir)

    numpy.random.seed(seed=int(time.time()))

    # first simulate the data
    #x = ss.runSim_DisjointSets(datadir, ngenes=ngenes, nparts=nparts, nsamples=nsamples, deltad=deltad)
    x = [datadir,"scorematrix.tsv", 'exprdat.tsv', 'phenotype.tsv', "setmatrix.tsv"]

    # filter the data
    #y = fs.filterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])
    y = ['filtered_files_list.txt']

    # building empirical subgraph distribution
    #s = es.allSubgraphs(x[0],x[1],60,200)
    s = 'all_subgraphs.txt'

    # recover the trees
    z = dg.denovoGeneSets(filelist=y[0], dirs=x[0], outputprefix=denovoPrefix, adjmat=x[1])

    # filter trees and extract data for modeling
    trees, genes, means = cf.treeFilterAndEx(dirs=x[1], treefile=z[0], filterfiles=y[0], levelThresh=levelThresh, topNTrees=topNTrees)

    # run models
    m = mm.rfModel(dirs=x[0], exprfile=x[2], pheno=x[3], genes=genes, cvs=crossVal)

    # score the gene sets.
    out = scr.setScoringDenovo(dir=datadir, Nf=Nf, exprfile=x[2], subgraphfile=s, filterfiles=y[0], genes=genes)

    # compare model results to simulation.
    g = an.analysis(predacc=m, genes=genes, trees=trees, means=means, dirs=x[0], setfile=x[4], setscores=out)

    return(out)
