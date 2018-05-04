
import numpy as np
import time
import simulateData_DisjointSets as ss
import filterSignal as fs
import randomForest_model as mm
import analysisStandardSimulation as an
import extractSubGraphs as es
import setScoring as scr


def buildListOfGenesFromSetMat(datadir, filename):
    setmat = open(datadir+filename,'r').read().strip().split('\n')
    setlist = []
    for si in setmat:
        bits = np.array(si.split('\t'))
        idx = np.where(bits == '1')
        setlist.append(idx[0])
    return(setlist)


def runStandard(datadir, Nf, subgraphFile):
    # defaults
    ngenes = 80    # number of nodes in the network
    nparts = 4     # number of sets in simulation
    nsamples = 10  # number of samples simulated
    filteredPrefix = "filtered_"  # file prefix for filtered files
    crossVal = 5   # random forest cross validation folds
    deltad = 5.0   # boost in the expression for target set
    Nf = int(Nf)   # number of scale levels for filtering
    numberSubGraphs = 100  # if generating subgraphs
    maxSubGraphSize = 25   # max size of subgraphs

    print('Running Standard Sim')
    print('\n    working in ' + datadir)

    np.random.seed(seed=int(time.time()))

    # first simulate the data
    x = ss.runSim_DisjointSets(datadir, ngenes=ngenes, nparts=nparts, nsamples=nsamples, deltad=deltad)

    # filter the data
    y = fs.filterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])

    # build the subgraph sets
    if subgraphFile == '':
        s = es.allSubgraphs(x[0],x[1],maxSubGraphSize,numberSubGraphs)
    else:
        s = subgraphFile

    # get a list of genes for each set in the setmatrix
    genes = buildListOfGenesFromSetMat(datadir, 'setmatrix.tsv')

    # run models
    m = mm.rfModel(dirs=x[0], exprfile=x[2], pheno=x[3], genes=genes, cvs=crossVal)

    # score the gene sets.
    out, samps = scr.setScoringStandard(dir=datadir, Nf=Nf, exprfile=x[2], subgraphfile=s, filterfiles=y[0], genes=genes)

    # compare model results to simulation.
    g = an.analysis(predacc=m, genes=genes, dirs=x[0], setfile=x[4], setscores=out, setsamples=samps)

    return(out)
