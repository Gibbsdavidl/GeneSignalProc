
# setScoring #

import numpy as np
import scipy as sp
import scipy.stats
import extractSubGraphs as sg



def setScoringDenovo(dir, Nf, exprfile, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    # rank each of the samples across gene sp
    inputs = open(dir + exprfile, 'r').read().strip().split("\n")
    outputs = []
    sampleList = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)
    sizeMax = len(sgs)

    for sample in range(1,len(inputs)):

        exprMat = inputs[sample].strip().split('\t')
        sampleList.append(exprMat[0])
        allRes = []

        # for each denovo gene set in the trees,
        for i, gs in enumerate(genes):

            # sum up the ranks (r_e).
            expr = [float(x) for j, x in enumerate(exprMat) if j > 0]
            exprRanks = sp.stats.rankdata(expr)
            rankSum = sum([exprRanks[j] for j in gs])

            # sum up ranks for sampled subgraphs mean(r_s).
            m = len(gs)

            if m <= sizeMax:
                # same for scale levels, get subgraph sets, sum ranks
                subGraphSums = [sum([exprRanks[j] for j in gx]) for gx in sgs[m]]
                # save r_e / r_s
                res0 = sum([1.0 for x in subGraphSums if rankSum > x]) / float(len(subGraphSums))

                allRes.append(res0)
            else:
                allRes.append(np.nan)  # should be a nan, maybe.
        outputs.append(allRes)

    return( (outputs,sampleList) )


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
        levelSet = levels[sample]
        sampleList.append(sample)
        sampRes = []

        # for each denovo gene set in the trees,
        for i, gs in enumerate(genes):
            # each gene set starts with a rankSum of 0
            # then, in the filtered file, do the rank sum for each level
            m = len(gs)
            rankSum = 0.0
            subGraphSums = np.array([0.0 for gx in sgs[m]])

            if m <= sizeMax:
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


def setScoringStandard(dir, Nf, exprfile, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    # rank each of the samples across gene sp
    inputs = open(dir + exprfile, 'r').read().strip().split("\n")
    outputs = []
    sampleList = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)
    sizeMax = len(sgs)

    for sample in range(1,len(inputs)):

        exprMat = inputs[sample].strip().split('\t')
        sampleList.append(exprMat[0])
        allRes = []

        # for each denovo gene set in the trees,
        for i, gs in enumerate(genes):

            # sum up the ranks (r_e).
            expr = [float(x) for j, x in enumerate(exprMat) if j > 0]
            exprRanks = sp.stats.rankdata(expr)
            rankSum = sum([exprRanks[j] for j in gs])

            # sum up ranks for sampled subgraphs mean(r_s).
            m = len(gs)

            if m <= sizeMax:
                # same for scale levels, get subgraph sets, sum ranks
                subGraphSums = [sum([exprRanks[j] for j in gx]) for gx in sgs[m]]
                # save r_e / r_s
                res0 = sum([1.0 for x in subGraphSums if rankSum > x]) / float(len(subGraphSums))

                allRes.append(res0)
            else:
                allRes.append(np.nan)  # should be a nan, maybe.
        outputs.append(allRes)

    return( (outputs,sampleList) )

