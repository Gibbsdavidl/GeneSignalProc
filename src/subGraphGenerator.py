



# this function is to sample subgraphs from the gene-gene functional network.

import igraph as ig
from multiprocessing import Pool
import numpy as np

# first the function that returns subgraphs of a given size

# Determine which nodes fall in sufficiently large connected components

def oneSubgraph( x ):
    (i, G, size) = x
    n = len(G.vs)
    comp = G.components(mode=ig.STRONG)
    compSizes = [sum(i == np.array(comp.membership)) for i in set(comp.membership)]
    valid = [compSizes[i] > size for i in comp.membership]  # only sampling from components that are large enough
    firstNode = np.random.choice(np.where(valid)[0], size=1)[0] # root of the subgraph
    used = [i == firstNode for i in range(0,n)] # Is this node selected?
    neigh = [i in G.neighbors(firstNode) for i in range(0,n)]  #
    z = [ni and (not ui) for ni,ui in zip(neigh, used)]
    for iter in range(1,size): # grow the subgraph to size
        if (sum(z) > 0): # if there's nodes still available
            newNode = np.random.choice(np.where(z)[0], size=1)[0] # then pick a new one
            used[newNode] = True                                  # mark it as used
            for gi in G.neighbors(newNode):                       # and get the new neighbors
                neigh[gi] = True
            z = [ni and (not ui) for ni, ui in zip(neigh, used)]  # what's left?
    return(np.where(used)[0])

def forestFire( x ):
    (i, G, size) = x
    p = 0.35 / (1.0-0.35)
    comp = G.components(mode=ig.STRONG)
    compSizes = [sum(i == np.array(comp.membership)) for i in set(comp.membership)]
    valid = [compSizes[i] > size for i in comp.membership]  # only sampling from components that are large enough
    firstNode = np.random.choice(np.where(valid)[0], size=1)[0] # root of the subgraph
    burning = [firstNode]                                       # node is burning
    queue = [firstNode]                                                  # neighbors
    while len(burning) < size and len(queue) > 0:
        # get node
        v = queue.pop(0)
        # get non-burned edges of node
        ws = [w for w in G.neighbors(v) if w not in burning]
        # number of edges to sample ...
        es = np.random.geometric(p=p, size=1)[0]
        if len(ws) > 0 and es > len(ws):
            sws = ws
        elif len(ws) > 0:
            sws = np.random.choice(a=ws, size=es, replace=False) # selected nodes to burn
        else:
            sws = []
        for si in sws:
            if len(burning) < size:
                burning.append(si)
                queue.append(si)
            else:
                break
    return(burning)


def writeAllSubgraphs(dirs, allSgs):
    fout = open(dirs+'all_subgraphs.txt','w')
    for ais in allSgs:
        if (len(ais) > 0):
            for aij in ais:
                fout.write(','.join([str(x) for x in aij])+'\n')
    fout.close()
    return()

def allSubgraphs(dirs, adjmat, maxSize, numGraphs):
    print("loading network")
    mat = np.loadtxt(dirs+adjmat, delimiter='\t')
    G = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")
    allSgs = [[] for i in range(0,maxSize)] # for each subgraph size
    print("searching for subgraphs")
    for gsize in range(5, maxSize):
        inputs = [(i, G, gsize) for i in range(0,numGraphs)]  # gather the inputs
        with Pool(3) as p:
            #sgs = p.map(oneSubgraph, inputs)             # find the subgraphs
            sgs = p.map(forestFire, inputs)
        allSgs[gsize] = sgs
    writeAllSubgraphs(dirs,allSgs)                            # write them out
    return('all_subgraphs.txt')


def loadSubGraphs(dir, subgraphFile):
    dat = open(dir+subgraphFile,'r').read().strip().split('\n')
    sgdict = dict()

    # get size of each subgraph
    ms = np.array([len(di.split(',')) for di in dat])
    nmax = max(ms)

    for i in range(2,nmax): # for each size class of subgraph
        idx = np.where(ms == i)[0]
        sets = [dat[i] for i in idx]
        splitsets = [ [int(y) for y in x] for x in map(lambda x: x.split(','), sets)]
        sgdict[i] = splitsets

    return(sgdict)
