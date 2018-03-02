# performing segmentation on the scale space

import numpy as np
import igraph as ig
import copy
from scipy import stats
from sklearn.cluster import KMeans

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


def connectedComponentLabeling(net, scalespace, level, bins, minsetsize):
    #
    # returns graph components, that include the specified seed
    # net, an igraph network
    # bins, number of bins to quantize
    # seed, node index
    # scale space
    # the matrix with rows as scales and columns as genes, matches the gene order
    # expr, the signal
    # level, which scale, which is the row of the filtered data.
    #
    # start with simply finding the best pair node.
    signal = scalespace[level]  # the signal at this level in the space
    #linbins = np.linspace(np.min(signal), np.max(signal), bins) # could get adaptive cuts here !!!!!
    # then we quantize the values.
    #qsig = np.digitize(x=signal, bins=linbins)
    # we have labels for each node.
    km = KMeans(n_clusters=bins, max_iter=500, n_init=21)
    km.fit(signal.reshape(-1, 1))
    qsig = km.labels_

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
    for li in labels:
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



def connectedComponentLabelingGroup(net, scalespaceList, level, bins, minsetsize):
    #
    # returns graph components, that include the specified seed
    # net, an igraph network
    # bins, number of bins to quantize
    # seed, node index
    # scale space
    # the matrix with rows as scales and columns as genes, matches the gene order
    # expr, the signal
    # level, which scale, which is the row of the filtered data.
    #
    # start with simply finding the best pair node.
    signal = scalespace[level]  # the signal at this level in the space

    # for each edge,
    #   compute a t-test
    #
    # !!!!!!!! make qsig !!!!!!!!!!!!

    #linbins = np.linspace(np.min(signal), np.max(signal), bins) # could get adaptive cuts here !!!!!
    # then we quantize the values.
    #qsig = np.digitize(x=signal, bins=linbins)
    # we have labels for each node.
    km = KMeans(n_clusters=bins, max_iter=500, n_init=21)
    km.fit(signal.reshape(-1, 1))
    qsig = km.labels_

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
                        nborqsig = qsig[nbors] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        # get list of neighbors with T statistics above cutoff!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        matched = nbors[nborqsig == qsig[p]]
                        # if there are neighbors that have same mean!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    for li in labels:
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


def segmentSpace(net, bins, msr, minsetsize):
    setList = []
    for scale in range(0,len(msr)): # for each scale
        res0 = connectedComponentLabeling(net, msr, scale, bins, minsetsize)  # get components that have keypts.
        setList.append(res0)
    return(setList)  # so every level will have some sets of nodes.


def segmentSpaceGroup(net, bins, filteredList, minsetsize):
    # going to segment each scale-level... using all the data.
    setList = []
    for scale in range(0,len(msr)): # for each scale
        res0 = connectedComponentLabelingGroup(net, filteredList, scale, bins, minsetsize)  # get components that have keypts.
        setList.append(res0)
    return(setList)  # so every level will have some sets of nodes.


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


def testfun():
    return("asdf")

