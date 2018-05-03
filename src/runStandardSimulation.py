
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


def runStandard(datadir, Nf):
    # defaults
    ngenes = 80
    nparts = 4
    nsamples = 10
    filteredPrefix = "filtered_"
    crossVal = 5
    deltad = 5.0
    Nf = int(Nf)
    numberSubGraphs = 100
    maxSubGraphSize = 25

    print('running full denovo sim')
    print('\nworking in ' + datadir)

    np.random.seed(seed=int(time.time()))

    # first simulate the data
    x = ss.runSim_DisjointSets(datadir, ngenes=ngenes, nparts=nparts, nsamples=nsamples, deltad=deltad)

    # filter the data
    y = fs.filterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])

    # build the subgraph sets
    # s = es.allSubgraphs(x[0],x[1],maxSubGraphSize,numberSubGraphs)
    s = 'all_subgraphs.txt'

    # get a list of genes for each set in the setmatrix
    genes = buildListOfGenesFromSetMat(datadir, 'setmatrix.tsv')

    # run models
    m = mm.rfModel(dirs=x[0], exprfile=x[2], pheno=x[3], genes=genes, cvs=crossVal)

    # score the gene sets.
    out, samps = scr.setScoringDenovo(dir=datadir, Nf=Nf, exprfile=x[2], subgraphfile=s, filterfiles=y[0], genes=genes)

    # compare model results to simulation.
    g = an.analysis(predacc=m, genes=genes, dirs=x[0], setfile=x[4], setscores=out, setsamples=samps)

    return(out)
