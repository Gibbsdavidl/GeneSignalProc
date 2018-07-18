
# create a set of genes, set membership, a weighted network, and
# create two groups of samples with simulated expression data.
#  ... right now the sets are disjoint ...

#simsets.py

# need a set of sets, to find an ordering over all.
# related to contigs, but with repeated sets.

import sys
import numpy as np
import string
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
    # we have
    # the number of genes
    # and the number of parts we want
    partSize = int(ngenes/nparts)  # number of genes per part
    netList = []
    npgenes = np.array(genes)
    geneList = []
    for i in range(0,nparts):
        thesegenes = npgenes[(i*partSize):(i*partSize + partSize)] # get these genes in this part
        # add some extra random genes here too.
        g = igraph.Graph.Tree(partSize, 3)                         # make a tree network
        g.vs["name"] = thesegenes                                  # name it
        netList.append(g)                                          # add this net-set to the list
        geneList.append(thesegenes)                                #
    return((netList, geneList) )


def setmatrix(genes, sets, setnames, nparts):
    m = np.empty( (sets+1, len(genes)+1), dtype='U128')
    partSize = int(len(genes) / nparts)
    curr = 0
    # for each set
    for si in range(0,sets):
        # get the set name
        m[si,0] = setnames[si]
        # then for each gene
        for gi in range(1,len(genes)+1):
            # if it's in the set, mark it with a '1'
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
    lastvs = 0
    for ni, net in enumerate(netList):                         # for each net representing each set
        allvs = [i + idx for i in net.vs.indices] # update the vertex id
        idx += len(allvs)
        g.add_vertices(allvs)
        if ni > 0:  # we're on the next component
            a = np.random.choice(lastvs) # sample one of last vs
            b = np.random.choice(allvs)  # sample one of all vs
            w = np.random.exponential()
            g.add_edge(a,b, weight=w) # make edge.
        lastvs = allvs # make last vs the new vs.
        for ei in net.es:   # for each edge in the net
            eit = ei.tuple       # this will index into the allvs
            targ = allvs[eit[0]] # target
            sour = allvs[eit[1]] # source
            ewgt = np.random.exponential() ## the edge weight
            g.add_edge(targ, sour, weight = ewgt) # add a new edge into g
    g.vs['name'] = genes
    for i in range(0,20):  # number of random edges to add
        a = np.random.choice(allvs)  # sample one of last vs
        b = np.random.choice(allvs)  # sample one of all vs
        w = np.random.exponential()
        if a != b:
            g.add_edge(a, b, weight=w)  # make edge.

    return(g)


def gen_means_and_sds(sets, idx, delta):
    # the idx gene set will be following a linear trend.
    # then there are two set means, one for each phenotype group.
    # set_means = [ [np.random.lognormal(mean=1, sigma=0.25, size=1)[0], np.random.lognormal(mean=1, sigma=0.25, size=1)[0]] for si in sets] # [np.random.sample() * 1 for si in sets]
    # then for
    set_means = [
        [3, 3] for si
        in sets]
    set_means[idx][1] += delta
    return(set_means)


def gen_expression(ngenes, sets, set_means, pheno, sigma, si, samplenames):
    # for each set, generate a mean value and a stddev
    # sets maps gene names to sets.
    gexpr = np.zeros(ngenes) # expr for each gene
    # need a vector of means with length equal to the number of nodes
    geneMeans = np.concatenate([np.repeat(set_means[i][pheno[si]], repeats=len(sets[i])) for i in range(0,len(set_means))])
    gexpr += (np.random.multivariate_normal(mean=geneMeans, cov=sigma))
    #gexpr += (np.random.normal(loc=geneMeans,scale=2 ))
    gexpr = [samplenames[si]] + [str(gi) for gi in gexpr]
    return(gexpr)


def runSim_DisjointSets(homedir, ngenes, nparts, nsamples, deltad):
    # first simulate the gene and gene sets
    print("generating data")

    # first get gene names, set names, and sample names
    genes = [''.join([string.ascii_lowercase[i] for i in np.random.randint(low=0, high=26, size=5)]) for j in range(0, ngenes)]
    setnames = [''.join([string.ascii_lowercase[i] for i in np.random.randint(low=0, high=26, size=5)]) for j in range(0, nparts)]
    samplenames = ['sample_'.join([string.ascii_lowercase[i] for i in np.random.randint(low=0, high=26, size=2)]) for j in range(0, nsamples)]

    # create a list of networks.
    (netList, geneList) = createEvenPartsNetworks(genes, ngenes, nparts)

    # then generate the set matrix (sema)
    sema = setmatrix(genes, nparts, setnames, nparts)
    # can insert gene names into rows, but have to have binary values as strings
    np.savetxt(X=sema, fmt='%s', delimiter='\t', fname=homedir+"setmatrix.tsv")

    # then print out the adj matrix
    g  = connectnetwork(netList, genes, geneList)
    gesc2 = g.get_adjacency(attribute='weight')
    # replace the Nones with 0s
    for i in range(0,len(gesc2.data)):
        for j in range(0,len(gesc2.data)):
            if gesc2[i,j] is None:
                gesc2[i,j] = 0
    np.savetxt(X=gesc2.data, fmt='%s', delimiter='\t', fname=homedir+"scorematrix.tsv")

    # need symmetric adjacency matrix with weights
    # then 'condition' the matrix
    for i in range(0,len(gesc2.data)):
         sumres = sum(filter(None, gesc2.data[i])) + 1.0  # Some NoneTypes in there
         gesc2[i,i] = sumres

    # then can generate random expression data from network
    # generate the set-based expression levels
    expr_profiles = [(["sampleID"]+genes)]
    # generate the phenotype as binary variable.
    pheno = [(np.random.choice([0,1], size=1))[0] for i in range(0,nsamples)]
    # get the means for each set
    set_means = gen_means_and_sds(setnames, 1, deltad)  # set 1 (second set) will follow linear trend with slope x
    for si in range(0,nsamples):
        # then simulate the expression
        expr_profiles.append(gen_expression(ngenes, geneList, set_means, pheno, np.array(gesc2.data), si, samplenames))

    # write out the set means for each time point
    np.savetxt(X=np.transpose(set_means), fmt='%s', delimiter='\t', fname=homedir + 'set_means' + '.tsv')
    # and write out the simulated gene expression
    np.savetxt(X=expr_profiles, fmt='%s', delimiter='\t', fname=homedir + 'exprdat' + '.tsv')
    # and phenotype
    np.savetxt(X=np.transpose(pheno), fmt='%s', delimiter='\t', fname=homedir + 'phenotype' + '.tsv')
    # done
    print("done with data simulation")
    return([homedir,"scorematrix.tsv", 'exprdat.tsv', 'phenotype.tsv', "setmatrix.tsv"])
