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
        x = q[0]          # dequeue the first node in the list
        testin = copy.deepcopy(inset) # our test list
        testin.append(x)              # add this neighbor
        testout_lists = ([net.neighbors(gi) for gi in testin])  # the neighbors for each node in the "IN-SET"
        testout = list(set([item for sublist in testout_lists for item in sublist])) # Take union... these are the "OUT-SET"
        thisDiff = meanDiff(signal[testin], signal[testout]) # abs value difference
        # what is the difference between the newly formed group (+x) and the new outside
        if thisDiff - meanList[len(meanList)-1] > eps: # if the difference has NOT DROPPED beyond eps #
            # then we want to keep it because it's close to the key pt.
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
    return( (level, inset, np.mean(signal[inset]), np.mean(signal[outset]), np.mean(signal[inset]) - np.mean(signal[outset]), meanList ) )


def meanDiff(inset, outset):
    return(abs(np.mean(inset) - np.mean(outset)))


def add1(inset, i):
    x = copy.deepcopy(inset)
    x.append(i)
    return(x)


def getNeighbors(net, inset):
    allNeighbors = ([net.neighbors(gi) for gi in inset])  # the neighbors for each node in the "IN-SET"
    allNeighbors = list(set([item for sublist in allNeighbors for item in sublist]))  # Take union...
    if len(allNeighbors) > 0:
        [allNeighbors.remove(i) for i in inset]
    return(allNeighbors)


def bnb_segment(net, seed, scalespace, level, eps):
    # branch and bound
    # returns a local segmentation, based around the specified seed
    # net, an igraph network
    # eps, 'close enough' threshold
    # seed, node index
    # scale space
    # the matrix with rows as scales and columns as genes, matches the gene order
    # expr, the signal
    # level, which scale, which is the row of the filtered data.
    #
    # start with simply finding the best pair node.
    signal = scalespace[level]  # the signal at this level in the space
    q = net.neighbors(seed)
    best = 100000.0
    keep = -1
    while (len(q) > 0):  # while we still have nodes in the queue
        qi = q.pop()
        x = np.abs(signal[seed] - signal[qi])  # seed is a keypt ... want most similar.
        if x < best : # then let's keep it.
            keep = qi
            best = x

    if keep > -1:
        bestSet = [seed, keep] # now we have a pair.
    else:
        bestSet = [seed]
    # need a queue of possible solutions.
    # this would be a list of 'inset' .. a set of connected nodes in the *in-set*
    # first group of possible solutions would be adding neighbors.
    allNeighbors = getNeighbors(net, bestSet)
    bestScore = meanDiff(signal[bestSet], signal[allNeighbors])
    q = [add1(bestSet, i) for i in allNeighbors]

    ## STUCK HERE ##
    while (len(q) > 0):  # while we still have nodes in the queue
        #print(len(q))
        x = q.pop()          # dequeue the first node in the list
        #if len(x) < (len(signal) * 0.5): # can not have a set larger than half the size of the network!
        y = getNeighbors(net, x) # get the new surrounding neighbors
        thisDiff = meanDiff(signal[x], signal[y]) # abs value difference
        if thisDiff > (bestScore-eps): # if it's 'close enough' to the best then we branch on this solution
            newq = [add1(x, i) for i in y] # try adding these neighbors send off in parallel?
            q += newq
            if thisDiff > bestScore: # but only keep it if it's the best
                bestSet = x
                bestScore = thisDiff
    inMean = np.mean(signal[bestSet])
    outMean = np.mean(signal[getNeighbors(net, bestSet)])
    return( (level, bestSet, inMean, outMean, bestScore) )


def getKeyPts(net, sig):
    keypts = []
    for i in range(0,len(sig)):
        neighborVals = np.array([sig[j] for j in net.neighbors(i)])
        if all(sig[i] < neighborVals) or all(sig[i] > neighborVals) :
            keypts.append(i)
    return(keypts)


def segmentSpace(net, eps, msr):
    setList = []
    for scale in range(0,len(msr)): # for each scale
        sig1 = copy.deepcopy(msr[scale, :]) # getting out the filtered values for lvl
        keypts = getKeyPts(net, sig1)
        for seed in keypts:  # for each keypoint, try to recover a set.
            res0 = bnb_segment(net, seed, msr, scale, eps) # get segments for this scale
            setList.append(res0)
    return(setList)


def connectSets(thisSetList):
    coupled = []
    # for each scale
    for si in range(0,len(thisSetList)):
        for ti in range(si,len(thisSetList)): # get the upper triangle
            s_tup = thisSetList[si]
            t_tup = thisSetList[ti]
            overlap = sum(np.in1d(s_tup[1], t_tup[1]))
            if overlap > 2 and si != ti: # if there's any overlap.
                if (s_tup[2] > 0 and  t_tup[2] > 0) or (s_tup[2] < 0 and  t_tup[2] < 0): # same direction
                    if np.abs(s_tup[0] - t_tup[0]) < 2:                                  # same or adj scale level
                        coupled.append( (overlap, si, ti, s_tup[0], t_tup[0], s_tup[2], t_tup[2]) )
    return(coupled)


def checkGroups(couple, grps):
    # check if one of the sets is already part of a group
    t1 = (couple[3],couple[1])
    t2 = (couple[4],couple[2])
    chk = [ (t1 in x) or (t2 in x) for x in grps]
    if len(chk) > 0:
        return(np.where(chk)[0])
    else:
        return([])

def groupCoupling(coupled, setList):
    # start at the top
    groups = []
    for i in range(0,len(coupled)): # for each coupling
        # if one of these sets is already part of a group
        idx = checkGroups(coupled[i], groups)
        if len(idx) > 0:
            # add it to the group
            groups[idx[0]].add( (coupled[i][3],coupled[i][1]) ) # tuple of scale-level and set ID
            groups[idx[0]].add( (coupled[i][4],coupled[i][2]) )
        # else create a new group.
        else:
            groups.append( set([ (coupled[i][3],coupled[i][1]), (coupled[i][4],coupled[i][2]) ]) )
    return(groups)


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


def compileResults(filenum, thisSetList, setGroups, msr):
    # thisSetList is a list of tuples, (level, gene_idx_list, mean_val, ...)
    # setGroups is a list of lists of tuples (level, set ID)
    res1 = []
    for chainID in range(0,len(setGroups)): # which chain ID
        for thisTuple in setGroups[chainID]: # ti will be a list of tuples
            setMean = thisSetList[thisTuple[1]][2]
            for geneID in thisSetList[thisTuple[1]][1]:
                res0 = [filenum, chainID, thisTuple[1], thisTuple[0], geneID, msr[thisTuple[0]][geneID], setMean]
                res1.append(res0)
    return(res1)


#res1 = segmentSpace(net, 0.5, 10, sig)
#resGini = gini(res1)
