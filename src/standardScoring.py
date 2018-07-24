
import numpy as np
import time
import simulateData_DisjointSets as ss
import filterSignal as fs
import randomForest_model as mm
import analysisStandardSimulation as an
import makeGraphs as es
import setScoring as scr
import ssGSEA as ssgsea
import sys


def buildListOfGenesFromSetMat(datadir, filename):
    setmat = open(datadir+filename,'r').read().strip().split('\n')
    setlist = []
    for si in setmat:
        bits = np.array(si.split('\t'))
        idx = np.where(bits == '1')
        setlist.append(idx[0])
    return(setlist)


def runStandard(datadir, Nf, filename, filterType, cores):
    # defaults
    np.random.seed(seed=int(time.time()))

    # read the genes,
    genes = open(datadir+filename,'r').read().strip().split()
    ngenes = len(genes)   # number of nodes in the network

    # read the data
    # each sample becomes a vector,
    # of genes as in order genes above
    nsamples = len(data)  # number of samples simulated

    crossVal = 8   # random forest cross validation folds
    Nf = int(Nf)   # number of scale levels for filtering

    # filter the data
    if filterType == 'heat':
        print('running heat filter..')
        y = fs.heatFilterData(exprfile=exprfile, dirs=datadir, outputprefix='_filtered.tsv', Nf=Nf, adjmat=datadir+filename+'_adjmat.tsv.gz')
    else:
        y = fs.mexFilterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])

    # build the subgraph sets
    if subgraphFile == '':
        s = es.allSubgraphs(x[0],x[1],maxSubGraphSize,numberSubGraphs,int(cores))
    else:
        s = subgraphFile

    # get a list of genes for each set in the setmatrix
    genes = buildListOfGenesFromSetMat(datadir, 'setmatrix.tsv')  # also could specify what gene sets wanted...
    # or proc a gmt file

    ssgseaScores = ssgsea.scoreSets(dirs=x[0], geneSets=genes, exprfile=x[2], omega=2)

    # score the gene sets.

    out, samps = scr.setScoringStandardMultiScaleZscoreV2(dir=datadir, Nf=Nf, subgraphfile=s, filterfiles=y[0], genes=genes, cores=int(cores))

    # run models
    score, cvscores, clf, featImp = mm.rfModelSetScores(dirs=x[0], inputs=out, pheno=x[3], genes=genes, cvs=crossVal)

    # run models
    gseascore, gseacvscores, gseaclf, gseafeatImp = mm.rfModelSetScores(dirs=x[0], inputs=ssgseaScores, pheno=x[3], genes=genes, cvs=crossVal)

    # compare model results to simulation.
    g = an.analysis(predacc=score, genes=genes, dirs=x[0], setfile=x[4], setscores=out, setsamples=samps, featureImp=featImp, gseaScore=gseascore)

    return(out)

