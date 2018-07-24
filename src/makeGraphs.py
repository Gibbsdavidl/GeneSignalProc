



# this function is to sample subgraphs from the gene-gene functional network.

import igraph as ig
from multiprocessing import Pool
import numpy as np
import copy
import gzip

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
    (i, G, size, seed, valid) = x
    np.random.seed(seed)
    p = 0.35 / (1.0-0.35)
    burning = []
    while len(burning) < size:
        firstNode = np.random.choice(np.where(valid)[0], size=1)[0] # root of the subgraph
        burning = [firstNode]                                       # node is burning
        queue = [firstNode]             # neighbors
        while len(queue) > 0:
            # get node
            v = queue.pop(0)
            # get non-burned edges of node
            ws = [w for w in G.neighbors(v) if w not in burning]
            # number of edges to sample ...
            es = max(1, np.random.geometric(p=p, size=1)[0])
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
    bret = copy.deepcopy(burning)
    return(bret)


def writeAllSubgraphs(dirs, allSgs, subgraphname):
    fout = gzip.open(dirs+subgraphname,'w')
    for ais in allSgs:
        if (len(ais) > 0):
            for aij in ais:
                fout.write((','.join([str(x) for x in aij])+'\n').encode())
    fout.close()
    return()

# from https://www.peterbe.com/plog/uniqifiers-benchmark
def f5(seq, idfun=None):
   # order preserving
   if idfun is None:
       def idfun(x): return repr(x)
   seen = {}
   result = []
   for item in seq:
       item.sort()
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result


def allSubgraphs(dirs, adjfile, genesetfile, maxSize, numGraphs, cores):
    print("loading network")
    mat = np.loadtxt(dirs+adjfile, delimiter='\t')
    G = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")
    comp = G.components(mode=ig.STRONG)
    allSgs = [[] for i in range(0,maxSize)] # for each subgraph size
    print("searching for subgraphs")
    for gsize in range(5, maxSize):
        compSizes = [sum(i == np.array(comp.membership)) for i in set(comp.membership)]  ## very slow ##
        valid = [compSizes[i] > size for i in comp.membership]  # only sampling from components that are large enough
        print("  working on subgraphs of size " + str(gsize))
        inputs = [(i, G, gsize, np.random.randint(low=1, high=999999999), valid) for i in range(0,numGraphs)]  # gather the inputs
        with Pool(cores) as p:
            sgs = p.map(forestFire, inputs)
        sgsidx = f5(sgs)
        allSgs[gsize] = sgsidx
    subgraphfilename = genesetfile+'_subgraphs.tsv.gz'
    writeAllSubgraphs(dirs, allSgs, subgraphfilename)                            # write them out
    return(subgraphfilename)


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


def procGMT(datadir,GMTfile):
    # return an adjmat given a gmt file
    allgenes = set()
    allsets  = set()
    genesets = dict()  # key is gene set
    setnames = dict()  # key is gene name

    for line in open(datadir+GMTfile, 'r').read().strip().split('\n'):
        bits = line.split('\t')
        genes = [bi for i,bi in enumerate(bits) if i > 1]
        gsname = bits[0]  # the gene set name
        genesets[gsname] = genes # all the genes that belong to this set
        allsets.add(gsname)
        for gi in genes:
            allgenes.add(gi)
            if gi in setnames:
                setnames[gi].append(gsname)
            else:
                setnames[gi] = [gsname]


    allsets = list(allsets)
    allgenes = list(allgenes)

    genebins = []  # list of binary set memberships

    # for each gene
    for gi in allgenes:

        # make a vector of zeros with length equal to setnames
        membership = [0.0 for si in allsets]

        # then get index of sets that a gene *DOES* have, and set those ones.
        for i,si in enumerate(allsets):
            if si in setnames[gi]:  # if gi is in set si
                membership[i] = 1.0

        genebins.append(membership)

    # then can make pairwise scores

    return( (allgenes, allsets, genesets, setnames, genebins) )


def makeAdjMat(allgenes, genebins, setthr):

    scrmat = []
    # for each pair of genes,
    for i,gi in enumerate(allgenes):
        scrvec = [0.0 for i in allgenes]  # all zeros unless
        for j,gj in enumerate(allgenes):
            if i != j:
                a = 0.0; b = 0.0; c=0.0; d=0.0
                # get a=intersection, b in one not in other, c in other not in one, d not in either
                for k in range(0,len(genebins[0])):
                    if genebins[i][k] == 1.0 and genebins[j][k] == 0.0:
                        c += 1.0
                    elif genebins[i][k] == 0.0 and genebins[j][k] == 1.0:
                        b += 1.0
                    elif genebins[i][k] == 1.0 and genebins[j][k] == 1.0:
                        a += 1.0
                    else:
                        d += 1.0

                # compute the score.  0.5*(a/(a+b)+a/(a+c))
                if a > setthr:
                    scr = 0.5*( (a/(a+b)) + (a/(a+c)))
                    scrvec[j] = scr
        scrmat.append(scrvec)
    return(scrmat)


def writeAdjAndAnnot(datadir, filename, adjmat, allgenes):
    fout = gzip.open(datadir+filename+'_adjmat.tsv.gz','wb')
    for ai in adjmat:
        ab = [str(x) for x in ai]
        fout.write(('\t'.join(ab)+'\n').encode())
    fout.close()

    fout = gzip.open(datadir+filename+'_genes.tsv.gz', 'wb')
    for bi in allgenes:
        fout.write((bi+'\n').encode())
    fout.close()
    return(datadir+filename+'_adjmat.tsv.gz')


def makeGraphs(datadir, numgraphs, maxgraphsize, genesetfile, threshold, numCores, adjfile, genefile):
    # build the subgraph sets

    # first have to transform a gmt file to a network, adj mat.
    # write adjmat and gene list for row/col labels

    if adjfile == '' and genefile == '':
        (allgenes, allsets, genesets, setnames, genebins) = procGMT(datadir, genesetfile)
        adjmat = makeAdjMat(allgenes, genebins, threshold)
        adjfile = writeAdjAndAnnot(datadir, genesetfile, adjmat, allgenes)
    #else:
    #    # read gene file to get all genes
    #    allgenes = gzip.open(datadir+genefile).read().strip().split()

    # then with that adjmat, we
    s = allSubgraphs(datadir,adjfile,genesetfile,int(maxgraphsize),int(numgraphs),int(numCores))

    return(1)