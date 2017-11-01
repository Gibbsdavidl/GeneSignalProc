
#simsets.py

# need a set of sets, to find an ordering over all.
# related to contigs, but with repeated sets.

import numpy as np
import string

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

def simsetdata(ngenes):
    # first we define the set of genes
    genes = [''.join([string.ascii_lowercase[i] for i in np.random.randint(low=0, high=26, size=5)]) for j in range(0,ngenes)]
    # then we define a number of sets that contain all genes
    idx  = 0
    sampsize = round(abs(np.random.normal(80,20)))
    set1 = [genes[i] for i in np.random.randint(low=0, high=len(genes), size=sampsize)]
    sets =[(set1)]
    while(sum([genes[i] in f(sets) for i in range(0,len(genes))]) < len(genes)):
      sampsize = round(abs(np.random.normal(80,20)))
      set1 = [genes[i] for i in np.random.randint(low=0, high=len(genes), size=sampsize)]
      sets.append(set1)
      idx += 1
    return((genes, sets))

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
            scoremat[gi,hi] = kulczynski2(m[:,gi], m[:,hi])
    return(scoremat)


def permscores(scrmat, setmat, perms):
    # want to generate simulated vectors
    # that have same number of set memberships,
    # but random set memberships
    #---
    # first get list of set-membership-counts
    for pi in perms:
        #

# first simulate the gene and gene sets
genes, sets = simsetdata(100)

# then generate the set matrix
sema = setmatrix(genes, sets)
np.savetxt(X=sema, delimiter='\t', fname="setmatrix.tsv")

# then score the gene-gene pairs
gesc = setscores(sema, genes)
np.savetxt(X=gesc, delimiter='\t', fname="scorematrix.tsv")

# then get the gene ordering
gord = fulljoin(sets)
np.savetxt(X=gord, fmt='%s', delimiter='\t', fname="geneorder.tsv")

print("done")
