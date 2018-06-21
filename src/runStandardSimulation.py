
import numpy as np
import time
import simulateData_DisjointSets as ss
import filterSignal as fs
import randomForest_model as mm
import analysisStandardSimulation as an
import subGraphGenerator as es
import setScoring as scr
import ssGSEA as ssgsea


def buildListOfGenesFromSetMat(datadir, filename):
    setmat = open(datadir+filename,'r').read().strip().split('\n')
    setlist = []
    for si in setmat:
        bits = np.array(si.split('\t'))
        idx = np.where(bits == '1')
        setlist.append(idx[0])
    return(setlist)


def runStandard(datadir, Nf, subgraphFile, filterType):
    # defaults
    ngenes = 100   # number of nodes in the network
    nparts = 5     # number of sets in simulation
    nsamples = 32  # number of samples simulated
    filteredPrefix = "filtered_"  # file prefix for filtered files
    crossVal = 8   # random forest cross validation folds
    deltad = 2.0   # boost in the expression for target set
    Nf = int(Nf)   # number of scale levels for filtering
    numberSubGraphs = 300  # if generating subgraphs
    maxSubGraphSize = 30   # max size of subgraphs

    print('Running Standard Sim')
    print('\n    working in ' + datadir)

    np.random.seed(seed=int(time.time()))

    # first simulate the data
    x = ss.runSim_DisjointSets(datadir, ngenes=ngenes, nparts=nparts, nsamples=nsamples, deltad=deltad)

    # filter the data
    if filterType == 'heat':
        print('running heat filter..')
        y = fs.heatFilterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])
    else:
        y = fs.mexFilterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])

    # build the subgraph sets
    if subgraphFile == '':
        s = es.allSubgraphs(x[0],x[1],maxSubGraphSize,numberSubGraphs)
    else:
        s = subgraphFile

    # get a list of genes for each set in the setmatrix
    genes = buildListOfGenesFromSetMat(datadir, 'setmatrix.tsv')  # also could specify what gene sets wanted...
    # or proc a gmt file

    ssgseaScores = ssgsea.scoreSets(dirs=x[0], geneSets=genes, exprfile=x[2], omega=2)

    # score the gene sets.
    out, samps = scr.setScoringStandardMultiScale(dir=datadir, Nf=Nf, subgraphfile=s, filterfiles=y[0], genes=genes)

    # run models
    score, cvscores, clf, featImp = mm.rfModelSetScores(dirs=x[0], inputs=out, pheno=x[3], genes=genes, cvs=crossVal)

    # run models
    gseascore, gseacvscores, gseaclf, gseafeatImp = mm.rfModelSetScores(dirs=x[0], inputs=ssgseaScores, pheno=x[3], genes=genes, cvs=crossVal)

    # compare model results to simulation.
    g = an.analysis(predacc=score, genes=genes, dirs=x[0], setfile=x[4], setscores=out, setsamples=samps, featureImp=featImp, gseaScore=gseascore)

    return(out)


    # score the gene sets.
    # out,samps = scr.setScoringDenovoMultiScale(dir=datadir, Nf=Nf, exprfile=x[2], subgraphfile=s, filterfiles=y[0], genes=genes, levels=levels)

    # run random forest using gene set scores
    # score, cvscores, clf, featImp = mm.rfModelSetScores(dirs=datadir, inputs=out, pheno=x[3], genes=genes, cvs=crossVal)