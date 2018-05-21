

# performing graph segmentation on the scale space

import numpy as np
import igraph as ig
import copy
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
import sklearn.preprocessing

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
    return(abs(np.mean(inset) - np.mean(outset)))



def connectSets(thisSetList, overlapSize):
    # here we find set overlaps across the scale-levels
    # and call them 'coupled'
    coupled = []
    for si in range(0, len(thisSetList)):     # for each scale
        ti = si+1                        # and the adjacent scale
        if ti < len(thisSetList):         # don't go past the stack.
            s_tup = thisSetList[si]       # the tuple from segmentation
            t_tup = thisSetList[ti]
            for i in range(0,len(s_tup[1])):          # for each component in the scale-level tuple
                for j in range(0,len(t_tup[1])):
                    overlap = sum(np.in1d(s_tup[2][i], t_tup[2][j]))   # if there's overlap in the nodes
                    if overlap > overlapSize:                           # above some threshold
                        if (s_tup[4][i] > 0 and t_tup[4][j] > 0) or (s_tup[4][i] < 0 and t_tup[4][j] < 0):  # same direction
                            coupled.append((overlap, si, ti, i, j, s_tup[2][i], t_tup[2][j], s_tup[3][i], t_tup[3][j]))
    return (coupled)


def checkGroups(couple, grps):
    # check if one of the sets is already part of a group
    t1 = (couple[1], tuple(couple[5])) # level, set
    t2 = (couple[2], tuple(couple[6])) # level, set
    chk = [ (t1 in x) or (t2 in x) for x in grps] # CHECK HERE
    if len(chk) > 0:
        return(np.where(chk)[0])
    else:
        return([])


def joinSets(coupled, setList):
    # then we create groups of unique nodes
    # that are connected or coupled across levels
    # groups are a list of lists
    # the lists being lists of tuples (level, set)
    groups = []

    for i in range(0,len(coupled)): # for each level

        # if one of these sets is already part of a group
        idx = checkGroups(coupled[i], groups)

        if len(idx) > 0:
            # add it to the group
            groups[idx[0]].add( (coupled[i][1], tuple(coupled[i][5])) ) # tuple of scale-level and set ID
            groups[idx[0]].add( (coupled[i][2], tuple(coupled[i][6])) )

        # else create a new group.
        else:
            groups.append( set([ (coupled[i][1],tuple(coupled[i][5])), (coupled[i][2],tuple(coupled[i][6])) ]) ) # set of tuples

    return(groups)


def compileResults(filenum, thisSetList, setGroups, msr):
    # thisSetList is a list of tuples, (level, gene_idx_list, mean_val, ...)
    # setGroups is a list of lists of tuples (level, set ID)
    res1 = []
    for chainID in range(0,len(setGroups)): # which chain ID
        for thisTuple in setGroups[chainID]: # ti will be a list of tuples
            #setMean = thisSetList[thisTuple[1]][2]
            for geneID in thisTuple[1]:
                res0 = [filenum, chainID, thisTuple[0], geneID, msr[thisTuple[0]][geneID]] #, setMean]
                res1.append(res0)
    return(res1)



################  Segmentation ###############################

def getNeighbors(net, inset):
    allNeighbors = ([net.neighbors(gi) for gi in inset])  # the neighbors for each node in the "IN-SET"
    allNeighbors = list(set([item for sublist in allNeighbors for item in sublist]))  # Take union...
    if len(allNeighbors) > 0:
        [allNeighbors.remove(i) for i in inset]
    return(allNeighbors)


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.nanmedian(arr)
    return np.nanmedian(np.abs(arr - med))

def boolit(x):
    xmed = np.nanmedian(x)
    xmad = mad(x) / 2.0
    if (abs(x) > (xmed + xmad)):
        return(1)
    else:
        return(0)

def connectedComponentLabeling(net, scalespace, level, minsetsize):
    #
    # returns graph components
    #
    # net, an igraph network
    # bins, number of bins to quantize
    # seed, node index
    # scale space
    # the matrix with rows as scales and columns as genes, matches the gene order
    # expr, the signal
    # level, which scale, which is the row of the filtered data.
    #
    # start with simply finding the best pair node.
    signal0 = sklearn.preprocessing.scale(scalespace[level])  # the signal at this level in the space
    signal1 = sklearn.preprocessing.scale(scalespace[level+1])
    signal = signal1 - signal0

    qsig = np.array([boolit(x) for x in signal])
    labels = np.array([-1 for x in range(0,len(qsig))])  # start out with each node having

    # for each node
    currlabel = 1
    for ni in range(0,len(qsig)):
        # if this node is not already labeled.
        if labels[ni] == -1:
            # then give this node the current label
            labels[ni] = currlabel
            #   start a queue, add this node to the queue
            nodequeue = [ni]
            #   while the queue is not empty
            while len(nodequeue) > 0:
                # pop off a node
                p = nodequeue.pop()
                # get the neighborhood of this node.
                nbors = np.array(net.neighbors(p))
                if len(nbors) > 0:
                    # and select the ones that don't have labels
                    nborlabels = labels[nbors]
                    nbors = nbors[nborlabels == -1]
                    # if there are some neighbors without labels
                    if len(nbors) > 0:
                        nborqsig = qsig[nbors]
                        # get list of neighbors in the same quantized level value
                        matched = nbors[nborqsig == qsig[p]]
                        # if there are neighbors in the same quant level
                        if len(matched) > 0:
                            # set the node labels to the smallest label.
                            labels[matched] = currlabel
                            # add those nodes to the queue
                            for mi in matched:
                                nodequeue.append(mi)
        # update the current label
        currlabel +=1


    # then sort out the components that have set sizes greater than minsetsize
    components = set()
    labelSet = set([li for i,li in enumerate(labels) if qsig[i] == 0])
    for li in labelSet:
        if sum(labels == li) > minsetsize:
            components.add(li)

    #   get the mean levels, the set of nodes, etc.
    setidxs = []
    sigvals = []
    meanval = []
    deviati = []
    for ci in list(components):
        idx = np.where(labels == ci)[0]
        setidxs.append(idx)
        sigvals.append(signal[idx])
        meanval.append(np.mean(signal[idx]))
        deviati.append(signal[idx] - np.mean(signal[idx]))

    return( (level, components, setidxs, sigvals, meanval, deviati) )



def segmentSpace(net, msr, minsetsize):
    setList = []
    for scale in range(0,len(msr)-1): # for each scale, and scale +1
        res0 = connectedComponentLabeling(net, msr, scale, minsetsize)  #
        setList.append(res0)
    return(setList)  # so every level will have some sets of nodes.


