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


def groupTest(inset, outset, eps):
    if len(inset) == 0:
        return(False)
    elif len(sig) < 2:
        return( np.mean(inset) - np.mean(outset) )
    else:
        res0 = stats.ttest_ind(inset, outset, equal_var=True)[0]
        return(res0)


def meanDiff(inset, outset):
    return(np.mean(inset) - np.mean(outset))


def segment(net, seed, scalespace, level, eps):
    # returns a local segmentation, based around the specified seed
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
    meanList = [0.0] # the list of mean differences
    while (len(q) > 0):  # while we still have nodes in the queue
        x = q[0]
        testin = copy.deepcopy(inset) # our test list
        testin.append(x)              # take one out
        testout_lists = ([net.neighbors(gi) for gi in testin])  # the neighbors for each node in the "IN-SET"
        testout = list(set([item for sublist in testout_lists for item in sublist])) # Take union... these are the "OUT-SET"
        thisDiff = meanDiff(signal[testin], signal[testout])
        # what is the difference between the newly formed group (+x) and the new outside
        if thisDiff - meanList[len(meanList)-1] > eps: # if the difference has improved beyond eps #
            # then we want to keep it
            meanList.append(thisDiff)  # keep this difference to compare to later
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
    return( (inset, np.mean(signal[inset]), np.mean(signal[outset]), np.mean(signal[inset]) - np.mean(signal[outset])) )

#segment(gra, 0.01, 37, scalespace, 2)


def segmentSpace(net, eps, msr, n):

    for scale in range(0,len(msr)): # for each scale

        sig1 = copy.deepcopy(msr[scale, :]) # getting out the filtered values for lvl
        sig2 = copy.deepcopy(msr[scale, :]) # getting out the filtered values for lvl
        sig2.sort()
        idxTop = np.where(sig1 >= sig2[len(sig2)-n])[0] # get the top n
        idxBot = np.where(sig1 <= sig2[n])[0]           # get the bottom n

        # get indices for the top and bot genes.
        seed = idxTop[0]
        res0 = segment(net, seed, msr, scale, eps) # get segments for this scale


    return()


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
