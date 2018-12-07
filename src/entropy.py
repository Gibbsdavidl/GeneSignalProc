

import gzip
from multiprocessing import Pool
import fit_neg_binom as nbinomfit
from scipy.stats import nbinom
import numpy as np

def formatExprData(dirs,exprfile):
    reformat = []
    samples  = []
    inputs = gzip.open(dirs+exprfile,'rt').read().strip().split("\n")
    header = inputs[0].strip().split('\t')
    genedict = {g:(i-1) for i,g in enumerate(header) if i > 0}
    for i in range(1,len(inputs)):
        bits = inputs[i].strip().split('\t')
        vals = [float(x) for i,x in enumerate(bits) if i > 0]
        samples.append(bits[0])
        reformat.append(vals)
    return( (reformat, samples, genedict) )  # get a matrix, rows==samples, list of samples, and index of genes.




def entropyScoringV1( inputv ):

    (dir, sample, sigs, gdict, genes) = inputv

    # read in the filtered file for sample..
    #sgs = sg.loadSubGraphs(dir, subgraphfile)

    sampRes = [] # list of scores across gene sets

    for i, gs in enumerate(genes):  # for each gene set
        vals = np.array([np.ceil(pow(2.0, x)) for x in sigs if pow(2.0,x) > 2])  # our count data, with init filtering
        f = nbinomfit.fit_nbinom(X=vals)            # fit model to our gene counts data
        rv = nbinom(f['size'], f['prob'])           # make an object to sample from
        gsfiltered = [gi for gi in gs if gi in gdict]
        gidx = [gdict[gi] for gi in gsfiltered]             # get the gene index
        scorevals = [pow(2.0,sigs[i]) for i in gidx]  # then get the values for genes in gs
        ps = rv.cdf(scorevals)                      # look up the probs
        score = sum([ (pi * np.log(pi)) for pi in ps])                   # compute entropy, genes are bins
        sampRes.append(score)
    return(sampRes)




def scoreSets(dirs, geneSets, exprfile, cores):

    # dir: the working directory
    # exprfile: the expression file, samples in rows.
    # genes: the gene sets
    # cores: cores

    # return a matrix of gene set scores (samples X gs)

    print("scoring sets")
    sigs, samps, gdict = formatExprData(dirs, exprfile)

    # make list of inputs ... what's needed for each sample
    inputs = []
    sampleList = []
    for sample in range(0,len(samps)):
        sampleList.append(sample)
        inputs.append( (dirs, sample, sigs[sample], gdict, geneSets) )  # gather the inputs
    with Pool(cores) as p:
        outputs = p.map(entropyScoringV1, inputs)

    return( (outputs,sampleList) )
