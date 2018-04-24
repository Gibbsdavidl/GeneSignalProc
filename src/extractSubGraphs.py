



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
    for gsize in range(2, maxSize):
        inputs = [(i, G, gsize) for i in range(0,numGraphs)]  # gather the inputs
        with Pool(3) as p:
            sgs = p.map(oneSubgraph, inputs)             # find the subgraphs
        allSgs[gsize] = sgs
    writeAllSubgraphs(dirs,allSgs)                            # write them out
    return(allSgs)


