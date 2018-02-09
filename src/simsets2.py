
#simsets.py

# need a set of sets, to find an ordering over all.
# related to contigs, but with repeated sets.

import numpy as np
import string
import copy
import igraph

# function to flatten lists
def f(E):
    if E==[]:
        return []
    elif type(E) != list:
        return [E]
    else:
        a = f(E[0])
        b = f(E[1:])
        a.extend(b)
    return


def binit(x, t):
    if x < t:
        return(1)
    else:
        return(0)


def setmat_to_numeric(setmat):
    cols = len(setmat[0])
    rows = len(setmat)
    m = np.zeros((rows, cols-1))
    for j in range(0,cols-1):
        for i in range(0,rows):
            m[i,j] = float(setmat[i,(j+1)])
    return(m)



def createEvenPartsNetworks(genes, ngenes, nparts):
    # we have gene
    # the number of genes
    # and the number of parts we want
    partSize = int(ngenes/nparts)  # number of genes per part
    netList = []
    npgenes = np.array(genes)
    geneList = []
    for i in range(0,nparts):
        thesegenes = npgenes[(i*partSize):(i*partSize + partSize)] # get these genes in this part
        g = igraph.Graph.Tree(partSize, 3)                         # make a tree network
        g.vs["name"] = thesegenes                                  # name it
        netList.append(g)                                          # add this net-set to the list
        geneList.append(thesegenes)                                #
    return((netList, geneList) )


def setmatrix(genes, sets, setnames):
    m = np.empty( (sets+1, len(genes)+1), dtype='U128')
    partSize = int(len(genes) / nparts)
    curr = 0
    for si in range(0,sets):
        m[si,0] = setnames[si]
        for gi in range(1,len(genes)+1):
            if gi > curr and gi < curr+partSize:
                m[si,gi] = '1'
            else:
                m[si,gi] = '0'
        curr = curr + partSize
    return(m)


def connectnetwork(netList, genes, geneList):
    # build up a new network from the parts.
    g = igraph.Graph()
    idx = 0
    partSize = int(len(genes) / len(geneList)) # have to update vertex names
    for net in netList:                         # for each net representing each set
        allvs = [i + idx for i in net.vs.indices] # update the vertex id
        idx += len(allvs)
        g.add_vertices(allvs)
        for ei in net.es:   # for each edge in the net
            eit = ei.tuple       # this will index into the allvs
            targ = allvs[eit[0]] # target
            sour = allvs[eit[1]] # source
            ewgt = np.random.exponential() ## the edge weight
            g.add_edge(targ, sour, weight = ewgt) # add a new edge into g
    g.vs['name'] = genes
    g.add_edges([(19,20), (39,40), (59,60), (79,80)]) # joining the parts ### AUTOMATE THIS!!
    return(g)


def gen_means_and_sds(sets, idx, jdx, slope):
    # the idx gene set will be following a linear trend.
    # the jdx indexes the time point
    set_means = [np.random.sample() * 0.0 for si in sets] # [np.random.sample() * 1 for si in sets]
    set_means[idx] = jdx * slope
    return(set_means)


def gen_expression(ngenes, sets, set_means, sigma):
    # for each set, generate a mean value and a stddev
    gexpr = np.zeros(ngenes) # expr for each gene
    # need a vector of means with length equal to the number of nodes
    geneMeans = np.concatenate([np.repeat(set_means[i], repeats=len(sets[i])) for i in range(0,len(set_means))])
    gexpr += (np.random.multivariate_normal(mean=geneMeans, cov=sigma))
    return(gexpr)




#####################
# Going to generate:#
# 1. number of sets, and number of genes
# 2. network for each set
# 3. set matrix for each set.
# 4. connect the networks
# 5. get the adjacency matrix.
# 6.
# 7.
############################################################

# first simulate the gene and gene sets
ngenes = 100
nparts = 5
ntime = 10
genes = [''.join([string.ascii_lowercase[i] for i in np.random.randint(low=0, high=26, size=5)]) for j in range(0, ngenes)]
setnames = [''.join([string.ascii_lowercase[i] for i in np.random.randint(low=0, high=26, size=5)]) for j in range(0, nparts)]

# create a list of networks.
(netList, geneList) = createEvenPartsNetworks(genes, ngenes, nparts)

# then generate the set matrix (sema)
sema = setmatrix(genes, nparts, setnames)
# can insert gene names into rows, but have to have binary values as strings
np.savetxt(X=sema, fmt='%s', delimiter='\t', fname="/Users/davidgibbs/Data/SimWaveDat/Go4/setmatrix.tsv")

# then print out the adj matrix
g  = connectnetwork(netList, genes, geneList)
gesc2 = g.get_adjacency(attribute='weight')
np.savetxt(X=gesc2.data, fmt='%s', delimiter='\t', fname="/Users/davidgibbs/Data/SimWaveDat/Go4/scorematrix.tsv")

# need symmetric adjacency matrix with weights
# then 'condition' the matrix
for i in range(0,len(gesc2.data)):
     sumres = sum(filter(None, gesc2.data[i])) + 1.0  # Some NoneTypes in there
     gesc2[i,i] = sumres

# replace the Nones with 0s
for i in range(0,len(gesc2.data)):
    for j in range(0,len(gesc2.data)):
        if gesc2[i,j] is None:
            gesc2[i,j] = 0

# then can generate random expression data from network
# add function to random expression.

# generate the set-based expression levels
expr_file_names = []
for si in range(1,ntime):
    expr_file_names.append('/Users/davidgibbs/Data/SimWaveDat/Go4/exprdat_'+str(si)+'.tsv')
    # get the means for each set
    set_means = gen_means_and_sds(setnames, 1, si, 3)  # set 1 (second set) will follow linear trend with slope 3
    # then simulate the expression
    gexpr = gen_expression(ngenes, geneList, set_means, np.array(gesc2.data))
    # write out the set means for each time point
    np.savetxt(X=np.transpose([set_means]), fmt='%s', delimiter='\t', fname='/Users/davidgibbs/Data/SimWaveDat/Go4/set_means_'+str(si)+'.tsv')
    # and write out the simulated gene expression
    np.savetxt(X=np.transpose([genes, gexpr]),  fmt='%s', delimiter='\t', fname='/Users/davidgibbs/Data/SimWaveDat/Go4/exprdat_'+str(si)+'.tsv')
# keep a list of all the files
np.savetxt(X=expr_file_names, fmt='%s', delimiter='\t', fname="/Users/davidgibbs/Data/SimWaveDat/Go4/filelist.tsv")

# done
print("done")
