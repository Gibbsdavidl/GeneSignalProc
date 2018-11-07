
# setScoring #

import numpy as np
import scipy as sp
import scipy.stats as stats
import makeGraphs as sg
import sys
from multiprocessing import Pool

def iciRule(genes, msr):

    ciList = []
    for i in range(0,len(msr)):
        chrvals = (msr[i]).strip().split('\t')
        vals = np.array([float(x) for i,x in enumerate(chrvals) if i in genes])
        est = np.mean(vals)
        sd  = np.std(vals)
        lu = [est-sd, est+sd]
        ciList.append(lu)

    # algorithm
    # goal: want largest window where max(L) <= min(U)
    # start with full size window... satified?
    # then trim window by one, and slide ... any location that's satisifed? take it.

    windowLength = len(ciList) - 1
    while windowLength > 1:
        for loci in range(0,len(ciList)-windowLength): # move the starting point
            window = [i for i in range(loci, loci+windowLength+1)]  # get the window
            lbar = max([x[0] for i,x in enumerate(ciList) if i in window])
            ubar = min([x[1] for i,x in enumerate(ciList) if i in window])
            if lbar <= ubar:
                return(window)
        windowLength -= 1
    return([])


def setoverlap(x,y):
    return(sum([1 for xi in x if xi in y]))


def zscore(gsExpr, x):
    # assume gsExpr and x have same length
    z = (np.mean(gsExpr) - np.mean(x)) / np.std(x)
    return(z)


def iqrsum(x):
    # return the sum of values within the IQR
    y = [xi for xi in x if xi >= 0]
    #a = np.percentile(y, 25)
    #b = np.percentile(y, 75)
    #res0 = np.sum([yi for yi in y if yi >= a and yi <= b])
    res0 = np.median(y)
    return(res0)


def getProp(x, subGraphSums, i):
    # get the proportion of subgraphs with lower IQR sums.
    res0 = 0.0
    n = float(len(subGraphSums) - 1)
    for j,sgs in enumerate(subGraphSums):
        if i != j and x >= sgs:
            res0 += 1.0
    return(res0/n)


def buildPrior(inputs, genes, sgs, t, levelSet):
    # for each scale level
    sumsList = []  # list of summed IQRs of each sampled subgraph
    priorList = [] # list of fit beta parameters for each scale level and gene set size
    for li in levelSet:
        # build the priors for each gene set of interest
        exprMat = inputs[li].strip().split('\t')  # the filtered data
        expr = stats.rankdata([float(x) for x in exprMat])  # convert to floads
        gssum = dict()
        prior = dict()
        for gs in genes:
            # for each gene set, get size, get subgraphs, get proportions, fit beta
            m = len(gs)
            if m not in prior:
                subgraphs = [sgi for sgi in sgs[m] if setoverlap(sgi, gs) < int(t * m)]
                subGraphSums = np.array([ iqrsum([expr[j] for j in gx if expr[j] >= 0.0]) for gx in subgraphs])
                propList = [getProp(x, subGraphSums, i) for i,x in enumerate(subGraphSums)]
                gssum[m] = subGraphSums
                prior[m] = stats.beta.fit(propList)
        sumsList.append(gssum)  # sums for this level
        priorList.append(prior) # priors for this level
    return( (sumsList, priorList) )


def sampleScoringEB( inputv ):

    (dir, outdir, sample, inputFiles, subgraphfile, genes, t) = inputv

    # read in the filtered file for sample..
    inputs = open(outdir + inputFiles[sample], 'r').read().strip().split("\n")
    sgs = sg.loadSubGraphs(dir, subgraphfile)

    if len(inputs) == 1:
        toplevel = 1
    else:
        # levelSet = iciRule(gs, inputs)
        toplevel = int(len(inputs) / 3)

    levelSet = range(0, toplevel)

    print('building prior')
    (sumsList, priorList) = buildPrior(inputs, genes, sgs, t, levelSet)

    sampRes = [] # list of scores across gene sets

    for i, gs in enumerate(genes):  # for each gene set
        m = len(gs)
        if m in sgs:
            propSummary = 1.0
            for li in levelSet:
                exprMat = inputs[li].strip().split('\t')  # the filtered data
                expr = stats.rankdata([float(x) for x in exprMat])      # convert to floads
                gsSum  = iqrsum([expr[j] for j in gs if expr[j] >= 0.0])  # for this gene set... only genes we have measured.
                gsHigher = np.sum( [1.0 for xi in sumsList[li][m] if gsSum >= xi] ) # for sumsList in level li and set size m
                sgsN = float(len( sumsList[li][m]))
                alpha = priorList[li][m][0]
                beta  = priorList[li][m][1]
                gsProp = (gsHigher + alpha) / (sgsN + alpha + beta)
                propSummary *= gsProp # multiply it in.
            sampRes.append(propSummary)  # could take max here too
        else:
            sampRes.append(0.0)
    return(sampRes)


def setScoringStandardMultiScaleZscoreV2(dir, outdir, Nf, filterfiles, subgraphfile, genes, cores, threshold):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    inputFiles = open(outdir+filterfiles, 'r').read().strip().split('\n')
    sampleList = []

    # make list of inputs ... what's needed for each sample
    inputs = []
    for sample in range(0,len(inputFiles)):
        sampleList.append(sample)
        inputs.append( (dir, outdir, sample, inputFiles, subgraphfile, genes, float(threshold))  )  # gather the inputs
    with Pool(cores) as p:
        outputs = p.map(sampleScoringEB, inputs)

    return( (outputs,sampleList) )

#######################################################################################################
#######################################################################################################
#######################################################################################################


def a__sampleScoringZV2( inputv ):

    (dir, outdir, sample, inputFiles, subgraphfile, genes, t) = inputv

    # read in the filtered file for sample..
    inputs = open(outdir + inputFiles[sample], 'r').read().strip().split("\n")
    sampRes = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)

    for i, gs in enumerate(genes):  # for each gene set

        #levelSet = iciRule(gs, inputs)
        if len(inputs) == 1:
            toplevel = 1
        else:
            toplevel = int(len(inputs) / 2)
        levelSet = range(0, toplevel)
        m = len(gs)
        zs = []
        if m in sgs:
            subgraphs = [sgi for sgi in sgs[m] if setoverlap(sgi, gs) < int(t * m)]
            for li in levelSet:
                exprMat = inputs[li].strip().split('\t')  # the filtered data
                expr = [float(x) for x in exprMat]        # convert to floads
                gsExpr = np.array([expr[j] for j in gs if expr[j] >= 0.0])  # for this gene set... only genes we have measured.
                subGraphExpr = np.array([ [expr[j] for j in gx if expr[j] >= 0.0] for gx in subgraphs])

                subGraphZs = [zscore(gsExpr, x) for x in subGraphExpr]  # list of z scores
                zs.append(np.sum(subGraphZs) / np.std(subGraphZs)) # each level has a z

                #gsMean = np.mean(gsExpr)
                #subgraphMean =  np.mean( [np.mean(x) for x in subGraphExpr] )
                #subgraphSD = np.std( [np.std(x) for x in subGraphExpr] )
                #zs.append( (gsMean - subgraphMean) / subgraphSD )

                #scoreList = [stats.ttest_ind(a=gsExpr, b=x, equal_var=False, nan_policy='omit')[0] for x in subGraphExpr]
                #zs.append(stats.ttest_1samp(scoreList, 0.0)[0]) # T tests not sensitive

                #x = [item for sublist in subGraphExpr for item in sublist]
                #zs.append(stats.ttest_ind(a=gsExpr, b=x, equal_var=False, nan_policy='omit')[0])  # terrible

            sampRes.append(np.sum(zs))  # could take max here too
        else:
            sampRes.append(0.0)
    return(sampRes)
