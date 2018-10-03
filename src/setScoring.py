
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
    if len(gsExpr) != len(x):
        print("ERROR: gene set size mismatch")
        return(0.0)
    z = (np.mean(gsExpr) - np.mean(x)) / np.std(x)
    return(z)


def sampleScoringZV2( inputv ):

    (dir, outdir, sample, inputFiles, subgraphfile, genes, t) = inputv

    # read in the filtered file for sample..
    inputs = open(outdir + inputFiles[sample], 'r').read().strip().split("\n")
    sampRes = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)

    for i, gs in enumerate(genes):  # for each gene set

        levelSet = iciRule(gs, inputs)
        #levelSet = range(0,len(inputs))
        m = len(gs)
        zs = []
        if m in sgs:
            subgraphs = [sgi for sgi in sgs[m] if setoverlap(sgi, gs) < int(t * m)]
            for li in levelSet:
                exprMat = inputs[li].strip().split('\t')  # the filtered data
                expr = [float(x) for x in exprMat]        # convert to floads
                gsExpr = np.array([expr[j] for j in gs if expr[j] >= 0.0])  # for this gene set... only genes we have measured.
                subGraphExpr = np.array([ [expr[j] for j in gx if expr[j] >= 0.0] for gx in subgraphs])

                subGraphZs = [zscore(gsExpr, x) for x in subGraphExpr]
                zs.append(sum(subGraphZs) / np.std(subGraphZs)) # each level has a z

                #gsMean = np.mean(gsExpr)
                #subgraphMean =  np.mean( [np.mean(x) for x in subGraphExpr] )
                #subgraphSD = np.std( [np.std(x) for x in subGraphExpr] )
                #zs.append( (gsMean - subgraphMean) / subgraphSD )

                #scoreList = [stats.ttest_ind(a=gsExpr, b=x, equal_var=False, nan_policy='omit')[0] for x in subGraphExpr]
                #zs.append(stats.ttest_1samp(scoreList, 0.0)[0]) # T tests not sensitive

                #x = [item for sublist in subGraphExpr for item in sublist]
                #zs.append(stats.ttest_ind(a=gsExpr, b=x, equal_var=False, nan_policy='omit')[0])  # terrible

            sampRes.append(np.max(zs))  # could take max here too
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
        outputs = p.map(sampleScoringZV2, inputs)

    return( (outputs,sampleList) )

#######################################################################################################
#######################################################################################################
#######################################################################################################
