
# setScoring #

import numpy as np
import scipy as sp
import scipy.stats
import subGraphGenerator as sg
import sys


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


def setScoringStandardMultiScale_median_diffs_t(dir, Nf, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    inputFiles = open(dir + filterfiles, 'r').read().strip().split('\n')
    outputs = []
    sampleList = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)
    sizeMax = len(sgs)

    for sample in range(0, len(inputFiles)):
        # read in the filtered file for sample..
        inputs = open(dir + inputFiles[sample], 'r').read().strip().split("\n")
        sampleList.append(sample)
        sampRes = []

        # for each gene set,
        for i, gs in enumerate(genes):
            levelSet = iciRule(gs, inputs)
            m = len(gs)
            if len(levelSet) > 0 and m <= sizeMax:
                subgraphs = [sgi for sgi in sgs[m] if setoverlap(sgi,gs) < 1]
                dist = np.array([0.0 for gx in subgraphs])
                for li in levelSet:
                    exprMat = inputs[li].strip().split('\t')  # not great name ... it's the filtered data
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    gsExpr = np.array([expr[j] for j in gs])  # for this scale-level

                    subGraphVal = np.array([np.median([expr[j] for j in gx]) for gx in subgraphs])
                    gsVal = np.median(gsExpr)
                    dist += [gsVal - x for x in subGraphVal]

                res0 = np.median(dist) / mad(dist)
                sampRes.append(res0)
                #res0 = scipy.stats.ttest_1samp(dist,0.0)
                #sampRes.append(res0[0])
            else:
                sampRes.append(0.0)
        outputs.append(sampRes)
    return ((outputs, sampleList))


def tstat_twosample_pooled(gsExpr, x):
    # assume gsExpr and x have same length
    if len(gsExpr) != len(x):
        print("ERROR: gene set size mismatch")
        return(0.0)
    n1 = len(gsExpr)-1
    top = np.square(np.mean(gsExpr) - np.mean(x))
    bot = (  (n1*np.var(gsExpr) + n1*np.var(x)) / (2*n1)) * np.sqrt(2.0 / len(gsExpr))
    t = top / bot
    return(t)


def setScoringStandardMultiScaleTwoSampleTPooled(dir, Nf, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    inputFiles = open(dir+filterfiles, 'r').read().strip().split('\n')
    outputs = []
    sampleList = []

    for sample in range(0,len(inputFiles)):
        # read in the filtered file for sample..
        inputs = open(dir + inputFiles[sample], 'r').read().strip().split("\n")
        sampleList.append(sample)
        sampRes = []
        sgs = sg.loadSubGraphs(dir, subgraphfile)
        sizeMax = len(sgs)

        for i, gs in enumerate(genes):  # for each gene set
            levelSet = iciRule(gs, inputs)
            m = len(gs)
            dist = []

            if m <= sizeMax:
                for li in levelSet:
                    exprMat = inputs[li].strip().split('\t')  # not great name ... it's the filtered data
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    gsExpr = np.array([expr[j] for j in gs]) # for this scale-level
                    # the subgraph means and sd depend on *this* sample and *this ICI* level set... hard to precompute that.
                    subGraphExpr = np.array([ [expr[j] for j in gx] for gx in sgs[m] ])
                    dist += [tstat_twosample_pooled(gsExpr,x) for x in subGraphExpr]

                res0 = np.mean(dist)
                sampRes.append(res0)
            else:
                sampRes.append(0.0)

        outputs.append(sampRes)
    return( (outputs,sampleList) )



def setScoringStandardMultiScaleNumpyT(dir, Nf, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    inputFiles = open(dir+filterfiles, 'r').read().strip().split('\n')
    outputs = []
    sampleList = []

    for sample in range(0,len(inputFiles)):
        # read in the filtered file for sample..
        inputs = open(dir + inputFiles[sample], 'r').read().strip().split("\n")
        sampleList.append(sample)
        sampRes = []
        sgs = sg.loadSubGraphs(dir, subgraphfile)
        sizeMax = len(sgs)

        for i, gs in enumerate(genes):  # for each gene set
            levelSet = iciRule(gs, inputs)
            m = len(gs)
            dist = []

            if m <= sizeMax:
                for li in levelSet:
                    exprMat = inputs[li].strip().split('\t')  # not great name ... it's the filtered data
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    gsExpr = np.array([expr[j] for j in gs]) # for this scale-level
                    # the subgraph means and sd depend on *this* sample and *this ICI* level set... hard to precompute that.
                    subGraphExpr = np.array([ [expr[j] for j in gx] for gx in sgs[m] ])
                    dist += [scipy.stats.ttest_rel(gsExpr, x) for x in subGraphExpr]

                res0 = sum([x[0] for x in dist] ) / float(len(dist))
                #res0 = np.median([x[0] for x in dist]) # not as good...
                sampRes.append(res0)
            else:
                sampRes.append(0.0)

        outputs.append(sampRes)
    return( (outputs,sampleList) )





def zscore(gsExpr, x):
    # assume gsExpr and x have same length
    if len(gsExpr) != len(x):
        print("ERROR: gene set size mismatch")
        return(0.0)
    top = (np.mean(gsExpr) - np.mean(x))
    bot = np.std(gsExpr)
    t = top / bot
    return(t)


def setScoringStandardMultiScaleZscore(dir, Nf, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    inputFiles = open(dir+filterfiles, 'r').read().strip().split('\n')
    outputs = []
    sampleList = []

    for sample in range(0,len(inputFiles)):
        # read in the filtered file for sample..
        inputs = open(dir + inputFiles[sample], 'r').read().strip().split("\n")
        sampleList.append(sample)
        sampRes = []
        sgs = sg.loadSubGraphs(dir, subgraphfile)
        sizeMax = len(sgs)

        for i, gs in enumerate(genes):  # for each gene set
            levelSet = iciRule(gs, inputs)
            m = len(gs)
            dist = []

            # compute sd across all sgs

            if m <= sizeMax:
                for li in levelSet:
                    exprMat = inputs[li].strip().split('\t')  # not great name ... it's the filtered data
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    gsExpr = np.array([expr[j] for j in gs]) # for this scale-level
                    # the subgraph means and sd depend on *this* sample and *this ICI* level set... hard to precompute that.
                    subGraphExpr = np.array([ [expr[j] for j in gx] for gx in sgs[m] ])
                    dist += [tstat_twosample_pooled(gsExpr,x) for x in subGraphExpr]

                res0 = np.mean(dist)
                sampRes.append(res0)
            else:
                sampRes.append(0.0)

        outputs.append(sampRes)
    return( (outputs,sampleList) )


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

def setScoringStandardMultiScale_median_diffs(dir, Nf, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
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
            levelSet = iciRule(gs, inputs)
            m = len(gs)
            if len(levelSet) > 0 and m <= sizeMax:

                dist = np.array([0.0 for gx in sgs[m]])
                for li in levelSet:

                    exprMat = inputs[li].strip().split('\t')  # not great name ... it's the filtered data
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    gsExpr = np.array([expr[j] for j in gs]) # for this scale-level

                    subGraphVal = np.array([ np.median([expr[j] for j in gx]) for gx in sgs[m] ])
                    gsVal = np.median(gsExpr)
                    dist += [ gsVal - x for x in subGraphVal ]

                res0 = np.median(dist) #/ mad(dist)
                sampRes.append(res0)
            else:
                sampRes.append(0.0)
        outputs.append(sampRes)
    return( (outputs,sampleList) )



def setScoringStandardMultiScale_median_diffs_prob_above_zero(dir, Nf, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
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
            levelSet = iciRule(gs, inputs)
            m = len(gs)
            if m <= sizeMax:
                dist = np.array([0.0 for gx in sgs[m]])
                for li in levelSet:

                    exprMat = inputs[li].strip().split('\t')  # not great name ... it's the filtered data
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    gsExpr = np.array([expr[j] for j in gs]) # for this scale-level

                    subGraphVal = np.array([ np.median([expr[j] for j in gx]) for gx in sgs[m] ])
                    gsVal = np.median(gsExpr)
                    dist += [ gsVal - x for x in subGraphVal ]

                res0 = sum([1.0 for x in dist if x > 0]) / float(len(dist))
                sampRes.append(res0)
            else:
                sampRes.append(0.5)
        outputs.append(sampRes)
    return( (outputs,sampleList) )




def setScoringStandardMultiScale_mahalanoibis(dir, Nf, filterfiles, subgraphfile, genes):
    # dir: the working directory
    # Nf: number of scales
    # exprfile: the expression file, samples in rows.
    # filterfiles: the list of filtered expression files
    # subgraphfile: the file listing sampled subgraphs
    # genes: the gene sets

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    # rank each of the samples across genes
    inputFiles = open(dir+filterfiles, 'r').read().strip().split('\n')
    outputs = []
    sampleList = []
    sgs = sg.loadSubGraphs(dir, subgraphfile)
    sizeMax = len(sgs)

    for sample in range(0,len(inputFiles)):
        inputs = open(dir + inputFiles[sample], 'r').read().strip().split("\n")
        sampleList.append(sample)
        sampRes = []

        for i, gs in enumerate(genes):
            levelSet = iciRule(gs, inputs)
            m = len(gs)
            if m <= sizeMax:
                dist = []
                for li in levelSet:

                    exprMat = inputs[li].strip().split('\t')  # not great name ... it's the filtered data
                    expr = [float(x) for x in exprMat]  ############## IN FITLER FILE, yep
                    gsExpr = np.array([expr[j] for j in gs]) # for this scale-level

                    subGraphExpr = np.array([ [expr[j] for j in gx] for gx in sgs[m] ])
                    try:
                        u = np.mean(subGraphExpr, axis=0)
                        z = np.vstack(subGraphExpr)
                        c = np.cov(z.T)
                        v = np.linalg.inv(c)  # might barf here
                        d = scipy.spatial.distance.mahalanobis(u, gsExpr, v)
                        dist += [d]
                    except:
                        dist += [0.0]
                res0 = sum(dist)
                sampRes.append(res0)
            else:
                sampRes.append(0.0)
        outputs.append(sampRes)
    return( (outputs,sampleList) )





def TESTsetScoringStandardMultiScale(dir, Nf, filterfiles, subgraphfile, genes):
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
                    #dist += [scipy.stats.ttest_rel(gsExpr, x) for x in subGraphExpr]
                    #dist += [scipy.spatial.distance.cosine(gsExpr, x) for x in subGraphExpr]
                    dist += [(np.mean(gsExpr) - np.mean(x))/np.std(x) for x in subGraphExpr]
                    #try:
                    #    u = np.mean(subGraphExpr, axis=0)
                    #    z = np.vstack(subGraphExpr)
                    #    c = np.cov(z.T)
                    #    v = np.linalg.inv(c)  # might barf here
                    #    d = scipy.spatial.distance.mahalanobis(u, gsExpr, v)
                    #    dist += [d]
                    #except:
                    #    dist += [0.0]
                # save r_e / r_s ### ACROSS LEVELS ####
                #fout.write(str(rankSum)+'\t'+'\t'.join([str(z) for z in subGraphSums])+'\n')

                #res0 = sum([1.0 for x in subGraphSums if rankSum > x]) / float(len(subGraphSums))
                #res0 = sum([x[0] for x in dist] ) / float(len(dist))
                #res0 = scipy.stats.ttest_1samp(dist,0.0)
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


