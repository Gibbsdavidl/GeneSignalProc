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
    for i in range(0, len(sig)):
        neighborVals = np.array([sig[j] for j in net.neighbors(i)])
        if all(sig[i] < neighborVals) or all(sig[i] > neighborVals):
            keypts.append(i)
    return (keypts)


def segmentSpace(net, eps, msr):
    setList = []
    for scale in range(0, len(msr)):  # for each scale
        sig1 = copy.deepcopy(msr[scale, :])  # getting out the filtered values for lvl
        keypts = getKeyPts(net, sig1)
        for seed in keypts:  # for each keypoint, try to recover a set.
            res0 = bnb_segment(net, seed, msr, scale, eps)  # get segments for this scale
            setList.append(res0)
    return (setList)



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
