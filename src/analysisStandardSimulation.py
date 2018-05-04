


# code to compare the ground truth to what's recovered from denovoGeneSets

import numpy as np



def writeOutputs(dir,sampleList,outputs,idx):

    fout = open(dir+'setscores.txt','w')
    setids = ['set_'+str(i) for i in idx]
    fout.write('\t'.join(['ID'] + setids) +'\n')
    for i in range(0,len(outputs)):
        outstr = []
        for j in idx:
            outstr.append(str(outputs[i][j]))
        outstr = [sampleList[i]] + outstr
        fout.write('\t'.join(outstr)+'\n')
    fout.close()
    return(1)


# trees, (in, out) where those are pointers to the denovo_trees file.
def analysis(predacc, genes, dirs, setfile, setscores, setsamples):
    # predacc - prediction accuracy from random forest
    # genes - list of gene sets
    # dirs - working directory
    # setfile - the matrix of assignments to sets
    # setscores - the file of set scores matrix
    # setsamples - the sample names to write into the set score matrix.


    # open the set assignment matrix
    mat = open(dirs + setfile,'r').read().strip().split('\n')
    seti = np.array(mat[1].split('\t'))

    means = open(dirs + "set_means.tsv",'r').read().strip().split('\n')
    means = (means[1]).split('\t')

    fout = open(dirs+'analyout.tsv','w')
    fout.write("set\taccr\tmean\tngenes\tgenes\n")
    print("set\taccr\tmean\tngenes\n")

    # put the sets in order of prediction ability
    predidx = np.argsort(predacc)
    meanidx = np.argsort(means)

    # trying to find the order to list the results..
    posSum = []
    for i in range(0,len(predacc)):
        posa = (np.where(predidx == i)[0])[0] # where did item i get placed?
        posb = (np.where(meanidx == i)[0])[0] # and item b?
        posSum.append(posa+posb)

    idx = [i for i in np.argsort(posSum)]

    for i in idx:
        seti = str(i)
        a = str(predacc[i])
        b = str(means[i])
        c = str(len(genes[i]))
        f = str(genes[i])
        fout.write('\t'.join([seti,a,b,c,f])+'\n')
        print('\t'.join([seti,a,b,c]))

    writeOutputs(dirs, setsamples, setscores, idx)

    return(1)

