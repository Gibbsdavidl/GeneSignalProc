# performing segmentation on the scale space

import numpy as np
import igraph as ig
#import cairo
from scipy import stats

def loadSignal(filename):
    signal = []
    fin = open(filename,'r').read().strip().split("\n")
    for line in fin:
        bits = line.split("\t")
        signal.append(np.log10(float(bits[0])+0.001))
    return(np.array(signal))


def loadGenes(filename):
    signal = []
    fin = open(filename,'r').read().strip().split("\n")
    for line in fin:
        bits = line.split("\t")
        signal.append(bits[2])
    return(np.array(signal))


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
    q = net.neighbors(seed)
    inset  = [seed]  # the inside of the segment
    willadd = []
    outset = []      # the outside of the segment
    signal = scalespace[level,:]  # the signal at this level in the space
    while (len(q) > 0):  # while we still have nodes in the queue
        x = q[0] # take the first one out
        if groupTest(signal[inset], signal[x], sm_eps, t_eps):
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
