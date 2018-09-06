
import numpy as np
import time
import filterSignal as fs
import randomForest_model as mm
import analysisStandardSimulation as an
import setScoring as scr
import ssGSEA as ssgsea
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
            #else:
            #    print('Missing gene from building gene sets: ' + gi)
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


def runStandard(datadir, Nf, exprfile, filterType, cores, subgraphs, genefile, genesets, adjmat, phenofile, threshold, edgeThreshold, outputdir):
    # defaults
    np.random.seed(seed=int(time.time()))

    crossVal = 8   # random forest cross validation folds
    Nf = int(Nf)   # number of scale levels for filtering
    threshold = float(threshold)

    # read the genes for the graph
    genes = gzip.open(datadir+genefile,'r').read().strip().split()

    # get a list of genes for each set in the gene set file
    genesetidx = buildListOfGenesFromGeneSets(datadir, genes, genesets)
    genesetsymbols, setnames = buildListOfSymbolsFromGeneSets(datadir, genes, genesets)

    # filter the data
    if Nf == 1:
        # then don't do a decomposition
        y = fs.noFilterData(exprfile=exprfile, dirs=datadir, outputprefix='_filtered.tsv', Nf=Nf, adjmat=adjmat, allgenes=genes, outdir=outputdir)
    elif Nf > 1: # and (filterType == '' or filterType == 'heat'):
        y = fs.heatFilterData(exprfile=exprfile, dirs=datadir, outputprefix='_filtered.tsv', Nf=Nf, adjmat=adjmat, allgenes=genes, edgeT=float(edgeThreshold), outdir=outputdir)
    else:
        print("Options error: please check your number of scale levels and filter type")

    # make this an option.
    # and just unprocessed data
    y1lev = fs.noFilterData(exprfile=exprfile, dirs=datadir, outputprefix='_not_filtered.tsv', Nf=Nf, adjmat=adjmat, allgenes=genes, outdir=outputdir)


    # generate score matrices
    msgsScores, samps = scr.setScoringStandardMultiScaleZscoreV2(dir=datadir, outdir=outputdir, Nf=Nf, subgraphfile=subgraphs, filterfiles=y[0], genes=genesetidx, cores=int(cores), threshold=threshold)
    ssgseaScores = ssgsea.parScoreSets(dirs=datadir, geneSets=genesetsymbols, exprfile=exprfile, omega=2, cores=int(cores))
    msgs1LevelScores, samps1Level = scr.setScoringStandardMultiScaleZscoreV2(dir=datadir, outdir=outputdir, Nf=1, subgraphfile=subgraphs, filterfiles=y1lev[0], genes=genesetidx, cores=int(cores), threshold=threshold)

    # write out the score matrices
    an.writeOutputsGSO(datadir,samps,ssgseaScores,'ssgsea_scores.tsv',outdir=outputdir)
    an.writeOutputsGSO(datadir,samps,msgsScores,  'msgs_scores.tsv',outdir=outputdir)
    an.writeOutputsGSO(datadir,samps1Level,msgs1LevelScores, 'msgs_1level_scores.tsv',outdir=outputdir)

    # and rank the gene sets using the given phenotypes.
    if phenofile != '':
        # run models
        score, cvscores, clf, featImp = mm.rfModelSetScores(dirs=datadir, inputs=msgsScores, pheno=phenofile, genes=genes, cvs=crossVal)
        score1level, cvscores1level, clf1level, featImp1level = mm.rfModelSetScores(dirs=datadir, inputs=msgs1LevelScores, pheno=phenofile,  genes=genes, cvs=crossVal)
        gseascore, gseacvscores, gseaclf, gseafeatImp = mm.rfModelSetScores(dirs=datadir, inputs=ssgseaScores, pheno=phenofile, genes=genes, cvs=crossVal)

        # compare model results
        g = an.analysis(predacc=score, genes=genesetsymbols, featureImp=featImp, gseaScore=gseascore, level1Score=score1level, outdir=outputdir)

    return(1)

