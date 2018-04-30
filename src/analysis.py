


# code to compare the ground truth to what's recovered from denovoGeneSets

import numpy as np


def jaccard(a,b):
    top = len(np.intersect1d(a,b))
    bot = len(np.union1d(a,b))
    return(float(top) / float(bot))

#
# trees, (in, out) where those are pointers to the denovo_trees file.
def analysis(predacc, genes, trees, means, dirs, setfile, setscores):

    # open the set assignment matrix
    mat = open(dirs + setfile,'r').read().strip().split('\n')
    seti = np.array(mat[1].split('\t'))

    # we are looking at set 1.
    geneidx = [x - 1 for x in np.where(seti == '1')[0]]

    fout = open(dirs+'analyout.tsv','w')
    fout.write("accr\tmean\tngenes\tJI\ttreeidx\tgenes\n")

    # geneList is the target set repeated, for comparison to each tree
    geneList = [list(geneidx) for i in range(0,len(genes))]
    corrJI = [jaccard(list(a),b) for (a,b) in zip(genes, geneList)]

    # put the trees in order of prediction ability
    # ranks = [i for i in range(0,len(predacc))]
    predidx = np.argsort(predacc)
    meanidx = np.argsort(means)
    jiidx = np.argsort(corrJI)

    # trying to find the order to list the results..
    posSum = []
    for i in range(0,len(predacc)):
        posa = (np.where(predidx == i)[0])[0] # where did item i get placed?
        posb = (np.where(jiidx == i)[0])[0] # and item b?
        posSum.append(posa+posb)

    orderedGenes = []

    for i in np.argsort(posSum):
        a = str(predacc[i])
        b = str(means[i])
        c = str(len(genes[i]))
        d = str(trees[i])
        e = str(corrJI[i])
        f = str(genes[i])
        g = str(setscores[i])
        fout.write('\t'.join([a,b,c,e,d,f,g])+'\n')
        print('\t'.join([a,b,c,e,d,g]))

    return(1)

