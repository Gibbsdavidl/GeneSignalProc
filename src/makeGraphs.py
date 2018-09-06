



# this function is to sample subgraphs from the gene-gene functional network.

import igraph as ig
from multiprocessing import Pool
import numpy as np
import scipy as sp
import copy
import gzip
import pickle
import gc
import monty.itertools as it

# first the function that returns subgraphs of a given size

# Determine which nodes fall in sufficiently large connected components

def forestFire( x ):
    (i, gpickle, size, seed, valid) = x
    G = pickle.load(open(gpickle, "rb"))
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
            if len(burning) >= size:
                break
    bret = copy.deepcopy(burning)
    return(bret)


def writeAllSubgraphs(dirs, allSgs, subgraphname):
    fout = gzip.open(dirs+subgraphname,'w')
    for ais in allSgs:
        if (len(ais) > 0):
            for aij in ais:
                fout.write((','.join([str(x) for x in aij])+'\n').encode('utf-8'))
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


def getMaxSizeList(maxSize):

    #maxSize must be a string
    # and either just an integer
    # or a list of numbers separated by commas
    # or a range defined by two numbers separated by a dash

    if ',' in maxSize:
        x = maxSize.strip().split(',')
        y = [int(m) for m in x]
        return(y)
    elif '-' in maxSize:
        x = maxSize.strip().split('-')
        y = [int(m) for m in x]
        z = range(y[0], y[1])
        w = [int(zi) for zi in z]
        return(z)
    else:
        x = range(1,int(maxSize))
        y = [int(xi) for xi in x]
        return(x)

def allSubgraphs(dirs, edgefile, genesetfile, maxSize, numGraphs, cores):
    print("loading network")
    spmat = sp.sparse.load_npz(dirs+edgefile)
    sources, targets = spmat.nonzero()
    weights = spmat[sources,targets]
    weights = np.array(weights)[0]
    G = ig.Graph.TupleList(edges=zip(sources,targets,weights), directed=False, weights=True)
    print("... graph loaded ... ")
    print("G nodes: " + str(G.vcount()))
    print("G edges: " + str(G.ecount()))

    print("... pruning edges ... ")
    dedges = G.es.select(weight_lt=0.1)
    G.delete_edges(dedges)
    print("... after pruning ... ")
    print("G nodes: " + str(G.vcount()))
    print("G edges: " + str(G.ecount()))

    # clear out memory
    print("... pickling graph ...")
    gpickle = dirs+edgefile+'_pickle.p'
    pickle.dump(G, file=open(gpickle, "wb"))
    sources = []; targets = []; weights = [];

    comp = G.components(mode=ig.STRONG)

    maxSizeList = getMaxSizeList(maxSize)

    allSgs = [[] for i in range(0,max(maxSizeList)+1)] # for each subgraph size
    G = []
    gc.collect()

    print("searching for subgraphs")
    for gsize in maxSizeList:
        compSizes = [sum(i == np.array(comp.membership)) for i in set(comp.membership)]  ## very slow ##
        valid = [compSizes[i] > gsize for i in comp.membership]  # only sampling from components that are large enough
        print("  working on subgraphs of size " + str(gsize))
        for kins in it.chunks(range(0, numGraphs), cores):
            print("... ... generating multi-data ...")
            inputs = [(i, gpickle, gsize, np.random.randint(low=1, high=999999999), valid) for i in kins]  # gather the inputs
            print("... ... jobs at the pool ...")
            with Pool(cores) as p:
                sgs = p.map(forestFire, inputs)
            print("... ... uniqueifying... ")
            allSgs[gsize] += sgs
    subgraphfilename = genesetfile+'_subgraphs.tsv.gz'
    allFilteredSgs = [f5(ai) for ai in allSgs]
    writeAllSubgraphs(dirs, allFilteredSgs, subgraphfilename)                            # write them out
    return(subgraphfilename)


def loadSubGraphs(dir, subgraphFile):

    # now we must check the number of elements in each line

    dat = gzip.open(dir+subgraphFile,'rt').read().strip().split('\n')


    sgdict = dict()

    # get size of each subgraph
    ms = np.array([len(di.split(',')) for di in dat])
    nmax = max(ms)+1

    for i in range(0,nmax): # for each size class of subgraph
        idx = np.where(ms == i)[0]
        sets = [dat[i] for i in idx]
        splitsets = [ [int(y) for y in x] for x in map(lambda x: x.split(','), sets)]
        sgdict[i] = splitsets

    return(sgdict)


def procGMT(datadir,GMTfile, okgenes):
    # return an adjmat given a gmt file
    allgenes = set()   # all gene names
    allsets  = set()   # all set names
    genesets = dict()  # key is gene set
    setnames = dict()  # key is gene name
    if len(okgenes) > 0:
        okgenedict = {v: k for k, v in enumerate(okgenes)}

    for line in open(datadir+GMTfile, 'r').read().strip().split('\n'):
        bits = line.split('\t')
        genes = [bi for i,bi in enumerate(bits) if i > 1]
        if len(okgenes) > 0:
            genes = [gi for gi in genes if gi in okgenedict]  ### FILTER OUT GENES THAT ARE NOT PROTEIN CODING OFFICIAL SYMBOLS
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

    pickle.dump((allgenes, allsets, genesets, setnames), open( datadir+"gmt_proc.p", "wb" ) )

    return( datadir+"gmt_proc.p" )


def computeScoreList (x):
    # here allgenes is a list of genes
    # tuples of (i,j)

    ## LOAD PICKLE

    (idxs, pickleName, setthr) = x
    (allgenes, allsets, genesets, setnames)= pickle.load(open(pickleName, "rb"))

    I = np.array([])
    J = np.array([])
    V = np.array([])
    for (i,j) in idxs:
        gi = allgenes[i]
        gj = allgenes[j]
        a = 0.0;
        b = 0.0;
        c = 0.0;
        # get a=intersection, b in one not in other, c in other not in one, d not in either (not used)
        setgi = set(setnames[gi])
        setgj = set(setnames[gj])
        a = float(len(setgi.intersection(setgj)))
        b = float(len(setgi.difference(setgj)))
        c = float(len(setgj.difference(setgi)))
        # compute the score.  0.5*(a/(a+b)+a/(a+c))
        if a > setthr:
            scr = 0.5 * ((a / (a + b)) + (a / (a + c)))
            if scr > 0.0:
                I = np.append(I, i)
                J = np.append(J, j)
                V = np.append(V, scr)
    return( (I,J,V) )


def makeEdgeList(pickleName, setthr, cores):

    (allgenes, allsets, genesets, setnames) = pickle.load(open(pickleName, "rb"))

    grpSize = 100000
    n = len(allgenes)
    allgenes = []; allsets = []; genesets = []; setnames = [];
    gc.collect()

    I = np.array([])
    J = np.array([])
    V = np.array([])

    idxgrp = []
    for x in it.chunks(it.iuptri(items=range(0, n)), grpSize):
        idxgrp.append(x)
        if len(idxgrp) >= cores:
            print("... ... starting multicore ...")
            inputs = [(xi, pickleName, setthr) for xi in idxgrp]  # gather the inputs
            with Pool(cores) as p:
                srclst = p.map(computeScoreList, inputs)
            for si in srclst:
                (i,j,v) = si
                if len(i) > 0 and len(i) == len(v):
                    I=np.append(I,i)
                    J=np.append(J,j)
                    V=np.append(V,v)
            idxgrp = []
            gc.collect()
    # then need to unpack the I,J,Vs
    sparseEdges = sp.sparse.coo_matrix(( V, (I, J)), shape=(n,n)).tocsr()
    return(sparseEdges)



def writeEdgesAndAnnot(datadir, filename, sparseEdges, pickleName):
    (allgenes, allsets, genesets, setnames) = pickle.load(open(pickleName, "rb"))

    # load in the adjacency matrix
    foutname = (datadir + filename + '_sparseMatrix.npz')
    sp.sparse.save_npz(foutname, sparseEdges, compressed=True)

    # the list of genes that maps the row names of the adjacency matrix
    fout = gzip.open(datadir+filename+'_genes.tsv.gz', 'wb')
    for bi in allgenes:
        fout.write((bi+'\n').encode('utf-8'))
    fout.close()

    return(filename+'_sparseMatrix.npz')


def makeGraphs(datadir, numgraphs, maxgraphsize, genesetfile, threshold, numCores, edgefile, genefile):
    # build the subgraph sets

    # first have to transform a gmt file to a network, adj mat.
    # write adjmat and gene list for row/col labels
    if genefile != '':
        okgenes = gzip.open(datadir+genefile, 'rt').read().strip().split('\n')

    if edgefile == '':
        print('...proc gmt...')
        pickleName = procGMT(datadir, genesetfile, okgenes)
        print('...finding edges...')
        sparseEdges = makeEdgeList(pickleName, float(threshold), int(numCores))
        print('...writing edgelist...')
        edgefile = writeEdgesAndAnnot(datadir, genesetfile, sparseEdges, pickleName)

    # then with that adjmat, we search for subgraphs
    print('...searching for subgraphs...')
    s = allSubgraphs(datadir,edgefile,genesetfile,maxgraphsize,int(numgraphs),int(numCores))

    return(1)



########################

# Archive

def makeAdjMat(allgenes, setnames, setthr):

    scrmat = np.zeros( (len(allgenes),len(allgenes)) )
    n = len(allgenes)
    # for each pair of genes,
    for i in range(0,n):
        if i % 500 == 0:
            print(i)
        for j in range(i,n):
            if i != j:
                gi = allgenes[i]
                gj = allgenes[j]
                a = 0.0; b = 0.0; c=0.0;
                # get a=intersection, b in one not in other, c in other not in one, d not in either (not used)
                setgi = set(setnames[gi])
                setgj = set(setnames[gj])
                a = float(len(setgi.intersection(setgj)))
                b = float(len(setgi.difference(setgj)))
                c = float(len(setgj.difference(setgi)))

                # compute the score.  0.5*(a/(a+b)+a/(a+c))
                if a > setthr:
                    scr = 0.5*( (a/(a+b)) + (a/(a+c)))
                    scrmat[i][j] = scr
                    scrmat[j][i] = scr

    return(scrmat)


def writeAdjAndAnnot(datadir, filename, adjmat, allgenes):
    foutname = (datadir+filename+'_adjmat.tsv.gz')
    np.savetxt(fname=foutname, X=adjmat, delimiter='\t')

    fout = gzip.open(datadir+filename+'_genes.tsv.gz', 'wb')
    for bi in allgenes:
        fout.write((bi+'\n').encode('utf-8'))
    fout.close()

    return(filename+'_adjmat.tsv.gz')


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
