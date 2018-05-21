

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
import simulateData_DisjointSets as ss
import filterSignal as fs
import blobFinder as dg
import treeFilterAndExtration as cf
import randomForest_model as mm
import analysisDenovoSimulation as an
import extractSubGraphs as es
import setScoring as scr
import ssGSEA as gsea

def runDenovoSim(datadir, Nf, subgraphFile, filterType):
    # defaults

    ngenes = 80    # number of nodes in the network
    nparts = 4     # number of sets in simulation
    nsamples = 24  # number of samples simulated
    filteredPrefix = "filtered_"  # file prefix for filtered files
    crossVal = 8   # random forest cross validation folds
    deltad = 4.0   # boost in the expression for target set
    Nf: int = int(Nf)   # number of scale levels for filtering
    numberSubGraphs = 100  # if generating subgraphs
    maxSubGraphSize = 25   # max size of subgraphs
    denovoPrefix = "denovo_trees"
    levelThresh = 4 # number of scale levels to transverse
    topNTrees = 20  # number of trees to return


    print('running full denovo sim')
    print('\nworking in ' + datadir)

    numpy.random.seed(seed=int(time.time()))

    # first simulate the data
    x = ss.runSim_DisjointSets(datadir, ngenes=ngenes, nparts=nparts, nsamples=nsamples, deltad=deltad)

    # filter the data
    if filterType == 'heat' or filterType == '':
        print("applying heat filter")
        y = fs.heatFilterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])
    elif filterType == 'mexi':
        y = fs.mexFilterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])

    # load the subgraphs (or generate them)
    if subgraphFile == '':
        s = es.allSubgraphs(x[0],x[1],maxSubGraphSize,numberSubGraphs)
    else:
        s = subgraphFile

    # recover the trees
    z = dg.blobs(filelist=y[0], dirs=x[0], outputprefix=denovoPrefix, adjmat=x[1])

    # filter trees and extract data for modeling
    trees, genes, means, levels, nlevs = cf.treeFilterAndEx(dirs=x[1], treefile=z[0], filterfiles=y[0], levelThresh=levelThresh, topNTrees=topNTrees)

    # score the gene sets.
    out,samps = scr.setScoringDenovoMultiScale(dir=datadir, Nf=Nf, exprfile=x[2], subgraphfile=s, filterfiles=y[0], genes=genes, levels=levels)

    # run random forest using gene set scores
    score, cvscores, clf, featImp = mm.rfModelSetScores(dirs=datadir, inputs=out, pheno=x[3], genes=genes, cvs=crossVal)

    # then run the ssGSEA code on the same data... add to the analysis
    #ssgsea = gsea.scoreSets(datadir, genes, x[2], 1.0)

    # compare model results to simulation.
    g = an.analysisDenovo(predacc=score, genes=genes, levels=nlevs, trees=trees, means=means, dirs=x[0], setfile=x[4], setscores=out, setsamples=samps, featureImp=featImp)

    return(out)


def runDenovoSimRerun(datadir, Nf, subgraphFile, filterType):
    #
    # here let's not rebuild the network and data.
    #
    # defaults
    ngenes = 80    # number of nodes in the network
    nparts = 4     # number of sets in simulation
    nsamples = 10  # number of samples simulated
    filteredPrefix = "filtered_"  # file prefix for filtered files
    crossVal = 10   # random forest cross validation folds
    deltad = 5.0   # boost in the expression for target set
    Nf = int(Nf)   # number of scale levels for filtering
    numberSubGraphs = 100  # if generating subgraphs
    maxSubGraphSize = 25   # max size of subgraphs
    denovoPrefix = "denovo_trees"
    levelThresh = 3 # number of scale levels to transverse
    topNTrees = 20  # number of trees to return

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
    s = subgraphFile

    # recover the trees
    z = dg.denovoGeneSets(filelist=y[0], dirs=x[0], outputprefix=denovoPrefix, adjmat=x[1])

    # filter trees and extract data for modeling
    trees, genes, means, levels = cf.treeFilterAndEx(dirs=x[1], treefile=z[0], filterfiles=y[0], levelThresh=levelThresh, topNTrees=topNTrees)

    # score the gene sets.
    out,samps = scr.setScoringDenovoMultiScale(dir=datadir, Nf=Nf, exprfile=x[2], subgraphfile=s, filterfiles=y[0], genes=genes, levels=levels)

    # run models
    score, cvscores, clf, featImp = mm.rfModelSetScores(dirs=x[0], inputs=out, pheno=x[3], genes=genes, cvs=crossVal)

    # compare model results to simulation.
    g = an.analysisDenovo(predacc=score, genes=genes, trees=trees, means=means, dirs=x[0], setfile=x[4], setscores=out, setsamples=samps, featureImp=featImp)

    return(out)
