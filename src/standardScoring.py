
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
import gzip


def buildListOfGenesFromGeneSets(datadir, genes, genesetfile):
    allgenesdec = [gi.decode('utf-8') for gi in genes]
    allgenesdict = {g:i for i,g in enumerate(allgenesdec)}
    setlist = []
    genesets = []

    for line in open(datadir+genesetfile, 'r').read().strip().split('\n'):
        bits = line.split('\t')
        genes = [bi for i,bi in enumerate(bits) if i > 1]
        genesets.append(bits[0])
        idx = []
        for gi in genes:
            if gi in allgenesdict:
                idx.append(allgenesdict[gi]) # for gene gi in this gene set, what index is it?
            else:
                print('Missing gene from building gene sets: ' + gi)
        setlist.append(idx)
    return(setlist)

def buildListOfSymbolsFromGeneSets(datadir, genes, genesetfile):
    setlist = []
    setnames = []
    for line in open(datadir+genesetfile, 'r').read().strip().split('\n'):
        bits = line.split('\t')
        genes = [bi for i,bi in enumerate(bits) if i > 1]
        setlist.append(genes)
        setnames.append(bits[0])
    return(setlist, setnames)


def runStandard(datadir, Nf, exprfile, filterType, cores, subgraphs, genefile, genesets, adjmat, phenofile, threshold, edgeThreshold):
    # defaults
    np.random.seed(seed=int(time.time()))

    # read the genes for the graph
    genes = gzip.open(datadir+genefile,'r').read().strip().split()
    ngenes = len(genes)   # number of nodes in the network
    geneidx = {i: gi for i, gi in enumerate(genes)}  # look up gene indices

    crossVal = 8   # random forest cross validation folds
    Nf = int(Nf)   # number of scale levels for filtering
    threshold = float(threshold)

    # filter the data
    if Nf == 1:
        # then don't do a decomposition
        y = fs.noFilterData(exprfile=exprfile, dirs=datadir, outputprefix='_filtered.tsv', Nf=Nf, adjmat=adjmat, allgenes=genes)
    elif Nf > 1: # and (filterType == '' or filterType == 'heat'):
        y = fs.heatFilterData(exprfile=exprfile, dirs=datadir, outputprefix='_filtered.tsv', Nf=Nf, adjmat=adjmat, allgenes=genes, edgeT=float(edgeThreshold))
    else:
        print("Options error: please check your number of scale levels and filter type")

    # get a list of genes for each set in the gene set file
    genesetidx = buildListOfGenesFromGeneSets(datadir, genes, genesets)
    genesetsymbols, setnames = buildListOfSymbolsFromGeneSets(datadir, genes, genesets)

    # generate score matrices
    ssgseaScores = ssgsea.parScoreSets(dirs=datadir, geneSets=genesetsymbols, exprfile=exprfile, omega=2, cores=int(cores))
    msgsScores, samps = scr.setScoringStandardMultiScaleZscoreV2(dir=datadir, Nf=Nf, subgraphfile=subgraphs, filterfiles=y[0], genes=genesetidx, cores=int(cores), threshold=threshold)

    y1lev = fs.noFilterData(exprfile=exprfile, dirs=datadir, outputprefix='_not_filtered.tsv', Nf=Nf, adjmat=adjmat, allgenes=genes)
    msgs1LevelScores, samps1Level = scr.setScoringStandardMultiScaleZscoreV2(dir=datadir, Nf=1, subgraphfile=subgraphs, filterfiles=y1lev[0], genes=genesetidx, cores=int(cores), threshold=threshold)

    an.writeOutputsGSO(datadir,samps,ssgseaScores,'ssgsea_scores.tsv')
    an.writeOutputsGSO(datadir,samps,msgsScores,  'msgs_scores.tsv')
    an.writeOutputsGSO(datadir,samps1Level,msgs1LevelScores, 'msgs_1level_scores.tsv')

    if phenofile != '':
        # run models
        score, cvscores, clf, featImp = mm.rfModelSetScores(dirs=datadir, inputs=msgsScores, pheno=phenofile, genes=genes, cvs=crossVal)

        # run models for 1 level
        score1level, cvscores1level, clf1level, featImp1level = mm.rfModelSetScores(dirs=datadir, inputs=msgs1LevelScores, pheno=phenofile,  genes=genes, cvs=crossVal)

        # run models for ssGSEA
        gseascore, gseacvscores, gseaclf, gseafeatImp = mm.rfModelSetScores(dirs=datadir, inputs=ssgseaScores, pheno=phenofile, genes=genes, cvs=crossVal)

        # compare model results to simulation.
        g = an.analysis(predacc=score, genes=genesetsymbols, dirs=datadir, setscores=msgsScores, setsamples=samps, featureImp=featImp, gseaScore=gseascore, level1Score=score1level)

    return(1)

