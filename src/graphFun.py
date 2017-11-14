# performing segmentation on the scale space

import numpy as np
import igraph as ig
import copy
from scipy import stats

def loadSignal(filename, header=1, column=1):
    signal = []
    fin = open(filename,'r').read().strip().split("\n")
    for line in fin:
        if header == 1:
            header = 0 # skip the header
        else:
            bits = line.split("\t")
            #signal.append(np.log10(float(bits[column])+0.001))
            signal.append(float(bits[column]))
    return(np.array(signal))


def loadGenes(filename, header=1):
    genes = []
    fin = open(filename,'r').read().strip().split("\n")
    for line in fin:
        if header == 1:
            header = 0 # skip header
        else:
            bits = line.split("\t")
            genes.append(bits[0])
    return(np.array(genes))


def loadData(netfile, sigfile):
    sig = np.loadtxt(sigfile)
    mat = np.loadtxt(netfile, delimiter='\t')
    net = ig.Graph.Adjacency(list(mat), mode="undirected")
    return( (sig, mat, net) )


def groupTest(sig, pt, sm_eps, t_eps):
    if len(sig) == 0:
        return(False)
    elif len(sig) < 2:
        return( (np.abs(sig - pt) < sm_eps)[0] )
    else:
        res0 = np.abs(stats.ttest_1samp(sig, pt)[0])
        return(res0 < t_eps)


def segment(net, sm_eps, t_eps, seed, scalespace, expr, level):
    # returns a local segmentation
    # net, an igraph network
    # sm_eps, the difference threshold, used before a t-test can
    # t_eps, t-stat threshold
    # seed, node index
    # scale space
    # the matrix with rows as scales and columns as genes, matches the gene order
    # expr, the signal
    # level, which scale, which is the row of the filtered data.
    #
    q = net.neighbors(seed) # can we get neighbors within a certain distance?
    inset  = [seed]         # the inside of the segment
    willadd = []            # list of nodes that passed the test
    outset = []             # the outside of the segment
    signal = scalespace[level,:]  # the signal at this level in the space
    while (len(q) > 0):  # while we still have nodes in the queue
        #x = q[0] # take the first one out
        testin = copy.deepcopy(inset)
        testin.append(q[0])
        testout_lists = ([net.neighbors(gi) for gi in testin])
        testout = list(set([item for sublist in testout_lists for item in sublist]))
        # try adding one
        # but then the outside needs to be updated
        # what is the difference between the newly formed group (+x) and the new outside
        if groupTest(signal[inset], signal[x], sm_eps, t_eps): # test here. #
            # then we want to keep it
            willadd.append(x)
            q.remove(x)
            nx = [i for i in net.neighbors(x) if i not in inset + outset]
            q += nx
        else:
            # add it to the outset
            q.remove(x)
            outset.append(x)
        if len(willadd) > 0:
            # then we have some nodes to add
            inset += willadd
        willadd = []
        #print("inset size: " + str(len(inset)) + "  mean: " + str(np.mean(signal[inset])) + "  sd: " +str(np.std(signal[inset] )))
    return( (inset, signal[inset], stats.rankdata(signal)[inset], expr[inset] ) )

#segment(gra, 0.01, 37, scalespace, 2)


def segmentSpace(net, sm_eps, t_eps, seed, msr, sig):
    segList = [] # the inset
    sigList = [] # the wavelet coefficients of the set
    ranList = [] # the rank of the gene in the set
    rawList = [] # the expression of the genes
    for i in range(0,len(msr)):
        res0 = segment(net, sm_eps, t_eps, seed, msr, sig, i)
        segList.append(res0[0])
        sigList.append(res0[1])
        ranList.append(res0[2])
        rawList.append(res0[3])
    return( (segList, sigList, ranList, rawList) )


def seedSpace(net, sm_eps, t_eps, seed, msr, sig):
    segList = [] # just the seeds
    sigList = [] # the wavelet coefficient
    ranList = [] # the rank of the expr
    rawList = [] # the actual expression
    for i in range(0,len(msr)):
        res0 = msr[i,seed]
        segList.append([seed])
        sigList.append([res0])
        ranList.append([stats.rankdata(sig)[seed]])
        rawList.append([sig[seed]])
    return( (segList, sigList, ranList, rawList) )


def gini(resList):
    giniList = []
    levs = len(resList[1])
    for li in range(0,levs):
        diffList = []
        n = len(resList[0][li])
        for i in range(0,n):
            for j in range(i,n):
                diffList.append(np.abs(resList[1][li][i] - resList[1][li][j]))
        top = sum(diffList)
        bot = 2*n*sum(resList[0][li])
        giniList.append( (top/bot) )
    return(giniList)


#res1 = segmentSpace(net, 0.5, 10, sig)
#resGini = gini(res1)
