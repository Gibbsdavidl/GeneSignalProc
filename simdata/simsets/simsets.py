
#simsets.py

# need a set of sets, to find an ordering over all.
# related to contigs, but with repeated sets.

import numpy as np
import string
import copy

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
        return a

def simsetdata(ngenes, mean_sample_size, sd_sample_size):
    # given a number of genes, generate gene sets, and then thin them out
    # for a smaller number of sets that still contains all genes
    # first we define the set of genes
    genes = [''.join([string.ascii_lowercase[i] for i in np.random.randint(low=0, high=26, size=5)]) for j in range(0,ngenes)]
    # then we define a number of sets that contain all genes
    idx  = 0
    sampsize = round(abs(np.random.normal(mean_sample_size,sd_sample_size)))
    set1 = [genes[i] for i in np.random.randint(low=0, high=len(genes), size=sampsize)]
    sets =[(set1)]
    while(sum([genes[i] in f(sets) for i in range(0,len(genes))]) < len(genes)):
      sampsize = round(abs(np.random.normal(mean_sample_size,sd_sample_size)))
      set1 = [genes[i] for i in np.random.randint(low=0, high=len(genes), size=sampsize)]
      sets.append(set1)
      idx += 1
    # after doing all the set creation...
    # could do some thinning .. drop sets if coverage over genes doesn't drop
    settry = copy.deepcopy(sets)
    for si in sets:
        settry.remove(si)
        # if all genes still represented, then stay dropped
        if sum([genes[i] in f(settry) for i in range(0,len(genes))]) == len(genes):
            print("good drop point, still got: " + str(sum([genes[i] in f(settry) for i in range(0,len(genes))])) + " genes")
        else:
            settry.append(si) # OK, put it back.
    return((genes, settry))

# Next we produce a matrix connecting sets.
def setoverlap(sets):
    mat = [] #np.zeros( (len(sets) * len(sets), 3) )
    idx = 0
    for i in range(0,len(sets)):
        for j in range(i,len(sets)):
            intsize = float(len(set(sets[i]).intersection(set(sets[j]))))
            unisize = float(len(set(sets[i]).union(set(sets[j]))))
            jaccard = intsize/unisize
            if i != j:
                mat.append([i,j,jaccard])
                print(str(i) + "  " + str(j))
                idx += 1
    mat.sort(key=lambda x: x[2])
    mat.reverse()
    return(mat)

# join sets
def setjoin (a,b):
    # returning [----a--][a&b][---b---]
    sa = set(a)
    sb = set(b)
    i = sa.intersection(sb)
    ad = sa.difference(sb)
    bd = sb.difference(sa)
    return(list(ad)+list(i)+list(bd))

# in each round of convalesce
def roundjoin(allsets, scores):
    setshave = [i for i in range(0,len(allsets))]
    setsused = []
    newsets = []
    for i in range(0,len(scores)):
        if (scores[i][0] not in setsused) and (scores[i][1] not in setsused):
            newsets.append(setjoin(allsets[scores[i][0]], allsets[scores[i][1]]))
            print("joined " + str(scores[i][0]) + "  " + str(scores[i][1]))
            setsused.append(scores[i][0])
            setsused.append(scores[i][1])
            setshave.remove(scores[i][0])
            setshave.remove(scores[i][1])
    for j in setshave:
        newsets.append(allsets[scores[i][0]])
    return(newsets)

def fulljoin(sets):
    while(len(sets) > 1):
        mat = setoverlap(sets)
        nextsets = roundjoin(sets,mat)
        mat = setoverlap(nextsets)
        sets = nextsets
    return(sets[0])

def setmatrix(genes, sets):
    m = np.zeros( (len(sets), len(genes)) )
    for gi in range(0,len(genes)):
        for si in range(0,len(sets)):
            if genes[gi] in sets[si]:
                m[si,gi] = 1
    return(m)


def kulczynski2(x, y):
    #    1  0
    # 1  a  b
    # 0  c  d
    a=0.0
    b=0.0
    c=0.0
    d=0.0
    for i in range(0,len(x)):
        if x[i] == 1 and y[i] == 1:
            a += 1.0
        if x[i] == 1 and y[i] == 0:
            b +=1.0
        if x[i] == 0 and y[i] == 1:
            c +=1.0
        if x[i] == 0 and y[i] == 0:
            d +=1.0
    numer = ((a/2.0)*(2*a+b+c))
    denom = ((a+b)*(a+c))
    if denom != 0:
        return(  numer/denom  )
    else:
        return(0)

def setscores(m, genes):
    # take a setmatrix
    scoremat = np.zeros( (len(genes),len(genes)) )
    for gi in range(0,len(genes)):
        for hi in range(0,len(genes)):
            if gi != hi:
                scoremat[gi,hi] = kulczynski2(m[:,gi], m[:,hi])
    return(scoremat)


def binit(x, t):
    if x < t:
        return(1)
    else:
        return(0)


def permscores(scrmat, setmat, ngenes, perms):
    # want to generate simulated vectors
    # that have same number of set memberships,
    # but random set memberships
    #---
    # first get list of set-membership-counts
    num_genes_in_each_set = [sum(x) for x in setmat]
    print(num_genes_in_each_set)
    allscrrs = []
    for pi in range(0,perms):
        p1 = float(np.random.choice(a=num_genes_in_each_set, size=1))/ngenes
        p2 = float(np.random.choice(a=num_genes_in_each_set, size=1))/ngenes
        set1 = [binit(np.random.sample(), p1) for x in range(0, ngenes)]
        set2 = [binit(np.random.sample(), p2) for x in range(0, ngenes)]
        scrr = kulczynski2(set1,set2)
        allscrrs.append(scrr)
    return(np.max(allscrrs))


def apply_threshold(gesc, cutoff):
    gesc2 = copy.deepcopy((gesc))
    gesc2[np.where(gesc2 <= cutoff)] = 0
    return(gesc2)


def gen_expression(gord, sets):
    # for each set, generate a mean value and a stddev
    set_means = [np.random.sample() * 10 for si in sets]
    set_sds   = [np.random.sample() * 5 for si in sets]
    gexpr = np.zeros(len(gord)) # expr for each gene
    nbexpr = np.zeros(len(gord))
    for i,g in enumerate(gord):
        for j,s in enumerate(sets):
            if g in s: # if this gene is in this set
                gexpr[i] += abs(np.random.normal(set_means[j], set_sds[j], size=1))
                nbexpr[i] += np.random.negative_binomial(n=100, p=set_means[j]/10.0, size=1)
    return((gexpr,nbexpr))


#####################
# Going to generate:#
# 1. number of sets
# 2. a binary set-membership matrix
# 3. similarity scores between genes, based on shared set membership
# 4. permutation based score to threshold the similarity scores
# 5. network based on similarity scores
# 6. a gene ordering based on joining sets
# 7. simulated expression values based on shared-set values.
############################################################

# first simulate the gene and gene sets
ngenes = 100
genes, sets = simsetdata(ngenes, 30, 30)

# then generate the set matrix (sema)
sema = setmatrix(genes, sets)
np.savetxt(X=sema, delimiter='\t', fname="setmatrix.tsv")

# then score the gene-gene pairs (gene scores gesc)
gesc = setscores(sema, genes)
# find the permutation based thresholds
cutoff = permscores(gesc, sema, ngenes, 100)
gesc2  = apply_threshold(gesc, cutoff)
np.savetxt(X=gesc2, delimiter='\t', fname="scorematrix.tsv")

# then get the gene ordering (gene ordering gord
gord = fulljoin(sets)
np.savetxt(X=gord, fmt='%s', delimiter='\t', fname="geneorder.tsv")

# generate the set-based expression levels
(gexpr,nbexpr) = gen_expression(gord, sets)
np.savetxt(X=np.transpose([gexpr,nbexpr]), delimiter='\t', fname='exprdat.tsv')

print("done")
