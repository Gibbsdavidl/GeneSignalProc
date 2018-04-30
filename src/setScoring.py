

# setScoring #


import numpy as np
import scipy as sp
import extractSubGraphs as sg

def setScoringDenovo(dir, Nf, exprfile, filterfiles, subgraphfile, genes):

    # rank each of the samples across gene sp
    inputs = open(dir + exprfile, 'r').read().strip().split("\n")

    sample = 1  # skip the header
    zzzzz = inputs[sample].strip().split('\t')
    sgs = sg.loadSubGraphs(dir, subgraphfile)
    allRes = []

    # for each denovo gene set in the trees,
    for i, gs in enumerate(genes):
        print(i)

        # sum up the ranks (r_e).
        expr = [float(x) for j, x in enumerate(zzzzz) if j > 0]
        exprRanks = sp.stats.rankdata(expr)
        rankSum = sum([exprRanks[j] for j in gs])

        # sum up ranks for sampled subgraphs mean(r_s).
        m = len(gs)

        # same for scale levels, get subgraph sets, sum ranks
        subGraphSums = [sum([exprRanks[j] for j in gx]) for gx in sgs[m]]
        subGraphMean = np.mean(subGraphSums)

        # save r_e / r_s
        res0 = sum([1.0 for x in subGraphSums if rankSum > x]) / float(len(subGraphSums))
        res1 = rankSum / subGraphMean

        allRes.append( (res0, res1) )

    return(allRes)



def setScoringStandard(exprFile, subraphFile, setFile):

    return(1)
