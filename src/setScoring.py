
# setScoring #

import numpy as np
import scipy as sp
import scipy.stats
import subGraphGenerator as sg



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
    return(-1)


def setScoringStandardMultiScale(dir, Nf, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
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

        # for each gene set,
        for i, gs in enumerate(genes):
            # each gene set starts with a rankSum of 0
            # then, in the filtered file, do the rank sum for
            # each level returned by ICI
            #fout = open(dir+'scorelog_'+str(sample)+'_'+str(i)+'.tsv','w')
            levelSet = iciRule(gs, inputs)
            m = len(gs)
            rankSum = 0.0
            dist = []
            if m <= sizeMax:
                #subGraphSums = np.array([0.0 for gx in sgs[m]])
                #dist = np.array([0.0 for gx in sgs[m]])
                for li in levelSet:

                    exprMat = inputs[li].strip().split('\t')  # not great name ... it's the filtered data
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    #exprRanks = sp.stats.rankdata(expr)

                    gsExpr = np.array([expr[j] for j in gs]) # for this scale-level

                    #rankSum += sum([exprRanks[j] for j in gs])
                    # sum up ranks for sampled subgraphs mean(r_s).
                    # same for scale levels, get subgraph sets, sum ranks
                    subGraphExpr = np.array([ [expr[j] for j in gx] for gx in sgs[m] ])
                    #subGraphSums += np.array([sum([exprRanks[j] for j in gx]) for gx in sgs[m]])  ######### IN FITLER FILE, yep
                    dist += [scipy.stats.ttest_rel(gsExpr, x) for x in subGraphExpr]
                    #dist += [scipy.spatial.distance.cosine(gsExpr, x) for x in subGraphExpr]

                # save r_e / r_s ### ACROSS LEVELS ####
                #fout.write(str(rankSum)+'\t'+'\t'.join([str(z) for z in subGraphSums])+'\n')

                #res0 = sum([1.0 for x in subGraphSums if rankSum > x]) / float(len(subGraphSums))
                #res0 = sum([x[0] for x in dist] ) / float(len(dist))
                res0 = np.mean(dist)
                sampRes.append(res0)
            else:
                sampRes.append(0.0)
            #fout.close()
            # end one gene set
        # end one sample
        #totRes = sum(sampRes)
        #sampRes = [x/totRes for x in sampRes]
        outputs.append(sampRes)
    return( (outputs,sampleList) )
