

import numpy as np


def jaccard(a,b):
    top = len(np.intersect1d(a,b))
    bot = len(np.union1d(a,b))
    return(float(top) / float(bot))


def analysis(predacc, genes, trees, means, dirs, setfile):

    mat = open(dirs + setfile,'r').read().strip().split('\n')
    seti = np.array(mat[1].split('\t'))
    geneidx = [x - 1 for x in np.where(seti == '1')[0]]

    fout = open(dirs+'analyout.tsv','w')
    fout.write("accr\tmean\tngenes\tJI\ttreeidx\n")

    ranks = [i for i in range(0,len(predacc))]
    predidx = np.argsort(predacc)
    meanidx = np.argsort(means)

    # for each entry
    posSum = []
    for i in range(0,len(predacc)):
        posa = (np.where(predidx == i)[0])[0] # where did item i get placed?
        posb = (np.where(meanidx == i)[0])[0]
        posSum.append(posa+posb)

    geneList = [list(geneidx) for i in range(0,len(genes))]
    corrJI = [jaccard(list(a),b) for (a,b) in zip(genes, geneList)]

    for i in np.argsort(posSum):
        a = str(predacc[i])
        b = str(means[i])
        c = str(len(genes[i]))
        d = str(trees[i])
        e = str(corrJI[i])
        fout.write('\t'.join([a,b,c,e,d])+'\n')

    return(0)