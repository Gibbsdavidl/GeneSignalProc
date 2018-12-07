
import numpy as np
import time
import entropy as ent1
import randomForest_model as mm
import analysisStandardSimulation as an
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


def runEntropy(datadir, Nf, exprfile, filterType, cores, subgraphs, genefile, genesets, adjmat, phenofile, threshold, edgeThreshold, outputdir):
    # defaults
    np.random.seed(seed=int(time.time()))

    crossVal = 8   # random forest cross validation folds
    threshold = float(threshold)

    # read the genes for the graph
    genes = gzip.open(datadir+genefile,'r').read().strip().split()

    # get a list of genes for each set in the gene set file
    genesetidx = buildListOfGenesFromGeneSets(datadir, genes, genesets)
    genesetsymbols, setnames = buildListOfSymbolsFromGeneSets(datadir, genes, genesets)


    # generate score matrices
    entScores, samps    = ent1.scoreSets(dirs=datadir, geneSets=genesetsymbols, exprfile=exprfile, cores=int(cores))
    ssgseaScores = ssgsea.parScoreSets(dirs=datadir, geneSets=genesetsymbols, exprfile=exprfile, omega=2, cores=int(cores))

    # write out the score matrices
    an.writeOutputsGSO(datadir,samps,ssgseaScores,'ssgsea_scores.tsv',outdir=outputdir)
    an.writeOutputsGSO(datadir,samps,entScores,   'entropy_scores.tsv',outdir=outputdir)

    # and rank the gene sets using the given phenotypes.
    if phenofile != '':
        # run models
        entrfscore, cvscores, clf, featImp = mm.rfModelSetScores(dirs=datadir, inputs=entScores, pheno=phenofile, genes=genes, cvs=crossVal)
        gseascore, gseacvscores, gseaclf, gseafeatImp = mm.rfModelSetScores(dirs=datadir, inputs=ssgseaScores, pheno=phenofile, genes=genes, cvs=crossVal)

        # compare model results
        g = an.analysis3(predacc=entrfscore, genes=genesetsymbols, featureImp=featImp, gseaScore=gseascore, outdir=outputdir)

    return(1)

