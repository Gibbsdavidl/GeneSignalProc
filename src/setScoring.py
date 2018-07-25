
# setScoring #

import numpy as np
import scipy as sp
import scipy.stats
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



def zscore(gsExpr, x, sigma):
    # assume gsExpr and x have same length
    if len(gsExpr) != len(x):
        print("ERROR: gene set size mismatch")
        return(0.0)
    top = (np.mean(gsExpr) - np.mean(x))
    z = top / sigma
    return(z)


def sampleScoringZV2( inputv ):

    (dir, sample, inputFiles, subgraphfile, genes) = inputv

    # read in the filtered file for sample..
    inputs = open(dir + inputFiles[sample], 'r').read().strip().split("\n")
    sampRes = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)
    sizeMax = 201 # len(sgs[len(sgs)][0])+3

    for i, gs in enumerate(genes):  # for each gene set
        print("std score, gene set "+str(i))
        #levelSet = iciRule(gs, inputs)
        levelSet = range(0,len(inputs))
        m = min(len(gs), 199)
        zs = []

        if m <= sizeMax:

            subgraphs = [sgi for sgi in sgs[(m-1)] if setoverlap(sgi, gs) < int(0.2 * m)]
            for li in levelSet:
                exprMat = inputs[li].strip().split('\t')  # the filtered data
                expr = [float(x) for x in exprMat]        # convert to floads
                gsExpr = np.array([expr[j] for j in gs])  # for this gene set
                subGraphExpr = np.array([[expr[j] for j in gx] for gx in subgraphs])

                gsMean = np.mean(gsExpr)
                subgraphMean =  np.mean( [np.mean(x) for x in subGraphExpr] )
                subgraphSD = np.std( [np.std(x) for x in subGraphExpr] )
                zs.append( (gsMean - subgraphMean) / subgraphSD )

            sampRes.append(np.mean(zs))
        else:
            sampRes.append(0.0)
    return(sampRes)


def setScoringStandardMultiScaleZscoreV2(dir, Nf, filterfiles, subgraphfile, genes, cores):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    inputFiles = open(dir+filterfiles, 'r').read().strip().split('\n')
    sampleList = []

    # make list of inputs ... what's needed for each sample
    inputs = []
    for sample in range(0,len(inputFiles)):
        sampleList.append(sample)
        inputs.append( (dir, sample, inputFiles, subgraphfile, genes)  )  # gather the inputs
    with Pool(cores) as p:
        outputs = p.map(sampleScoringZV2, inputs)

    return( (outputs,sampleList) )

#######################################################################################################
#######################################################################################################
#######################################################################################################


def zscore(gsExpr, x, sigma):
    # assume gsExpr and x have same length
    if len(gsExpr) != len(x):
        print("ERROR: gene set size mismatch")
        return(0.0)
    top = (np.mean(gsExpr) - np.mean(x))
    z = top / sigma
    return(z)


def sampleScoringZNOMS( inputv ):

    (dir, sample, inputs, subgraphfile, genes) = inputv

    # read in the filtered file for sample..
    sampRes = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)
    sizeMax = len(sgs)
    expr = [float(x) for i, x in enumerate(inputs.split('\t')) if i > 0]

    for i, gs in enumerate(genes):  # for each gene set
        levelSet = [0]
        m = len(gs)
        zs = []

        if m <= sizeMax:
            subgraphs = [sgi for sgi in sgs[m] if setoverlap(sgi, gs) < 1]
            for li in levelSet:

                gsExpr = np.array([expr[j] for j in gs])  # for this gene set
                subGraphExpr = np.array([[expr[j] for j in gx] for gx in subgraphs])

                gsMean = np.mean(gsExpr)
                subgraphMean =  np.mean( [np.mean(x) for x in subGraphExpr] )
                subgraphSD = np.std( [np.std(x) for x in subGraphExpr] )
                zs.append( (gsMean - subgraphMean) / subgraphSD )

            sampRes.append(np.mean(zs))
        else:
            sampRes.append(0.0)
    return(sampRes)


def setScoringStandardMultiScaleZscoreNOMS(dir, Nf, exprfile, subgraphfile, genes, cores):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    inputs = open(dir+exprfile,'r').read().strip().split("\n")
    sampleList = []
    inputList = []

    # make list of inputs ... what's needed for each sample
    for sample in range(1,len(inputs)):
        sampleList.append(sample)
        inputList.append( (dir, sample, inputs[(sample)], subgraphfile, genes)  )  # gather the inputs
    with Pool(cores) as p:
        outputs = p.map(sampleScoringZNOMS, inputList)


    return( (outputs,sampleList) )

#######################################################################################################
#######################################################################################################
#######################################################################################################


def setScoringDenovoMultiScale(dir, Nf, exprfile, filterfiles, subgraphfile, genes,levels):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    # rank each of the samples across genes
    #inputs = open(dir + exprfile, 'r').read().strip().split("\n")
    inputFiles = open(dir+filterfiles, 'r').read().strip().split('\n')
    outputs = []
    sampleList = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)
    sizeMax = len(sgs)

    for sample in range(0,len(inputFiles)):
        # read in the filtered file for sample..
        inputs = open(dir + inputFiles[sample], 'r').read().strip().split("\n")
        sampleList.append(sample)
        sampRes = []

        # for each denovo gene set in the trees,
        for i, gs in enumerate(genes):
            # each gene set starts with a rankSum of 0
            # then, in the filtered file, do the rank sum for each level
            levelSet = levels[i]
            m = len(gs)
            rankSum = 0.0

            if m <= sizeMax:
                subGraphSums = np.array([0.0 for gx in sgs[m]])

                for li in levelSet:

                    exprMat = inputs[li].strip().split('\t')

                    # sum up the ranks (r_e).
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    exprRanks = sp.stats.rankdata(expr)
                    rankSum += sum([exprRanks[j] for j in gs])

                    # sum up ranks for sampled subgraphs mean(r_s).
                    # same for scale levels, get subgraph sets, sum ranks
                    subGraphSums += np.array([sum([exprRanks[j] for j in gx]) for gx in sgs[m]]) ######### IN FITLER FILE, yep

                # save r_e / r_s ### ACROSS LEVELS ####
                res0 = sum([1.0 for x in subGraphSums if rankSum > x]) / float(len(subGraphSums))
                sampRes.append(res0)
            else:
                sampRes.append(0.0)
            # end one gene set
        # end one sample
        outputs.append(sampRes)
    return( (outputs,sampleList) )

