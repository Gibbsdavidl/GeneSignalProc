


# code to compare the ground truth to what's recovered from denovoGeneSets

import numpy as np



#def writeOutputs(dir,sampleList,outputs,idx):
#
#    fout = open(dir+'setscores.tsv','w')
#   setids = ['set_'+str(i) for i in idx]
#    fout.write('\t'.join(['ID'] + setids) +'\n')
#    for i in range(0,len(outputs)):
#        outstr = []
#        for j in idx:
#            outstr.append(str(outputs[i][j]))
#        outstr = [sampleList[i]] + outstr
#        fout.write('\t'.join(outstr)+'\n')
#    fout.close()
#    return(1)

def writeOutputs(dir,sampleList,outputs,idx, outdir):

    fout = open(outdir+'setscores.tsv','w')
    for i in range(0,len(outputs)):
        outstr = []
        for j in idx:
            outstr.append(str(outputs[i][j]))
        outstr = [str(sampleList[i])] + outstr
        fout.write('\t'.join(outstr)+'\n')
    fout.close()
    return(1)

def writeOutputsGSO(dir,sampleList,outputs, filesuffix, outdir):
    fout = open(outdir+filesuffix,'w')
    for i in range(0,len(outputs)):
        outstr = []
        for j in range(0,len(outputs[0])):
            outstr.append(str(outputs[i][j]))
        outstr = [str(sampleList[i])] + outstr
        fout.write('\t'.join(outstr)+'\n')
    fout.close()
    return(1)

# trees, (in, out) where those are pointers to the denovo_trees file.
def analysis(predacc, genes, featureImp, gseaScore, level1Score, outdir):
    # predacc - prediction accuracy from random forest
    # genes - list of gene sets
    # dirs - working directory
    # setscores - the file of set scores matrix
    # setsamples - the sample names to write into the set score matrix.
    # featureImp
    # gseaScore - scores from ssGSEA
    # level1Score - scores from 1 level

    fout = open(outdir+'analyout.tsv','w')
    fout.write("set\tmsgs\t1level\tgsea\tfeatimp\tngenes\tgenes\n")
    print("set\tmsgs\t1level\tgsea\tngenes\n")

    # put the sets in order of prediction ability
    predidx = np.argsort(featureImp)

    for i in predidx:
        seti = str(i)
        a = str(predacc)
        b = str(level1Score)
        d = str(gseaScore)
        #c = str(len(genes[i]))
        #f = str(genes[i])
        g = str(featureImp[i])
        fout.write('\t'.join([seti,a,b,d,g])+'\n')
        print('\t'.join([seti,a,b,d,g]))

    #writeOutputs(dirs, setsamples, setscores, predidx)

    return(1)


# trees, (in, out) where those are pointers to the denovo_trees file.
def analysis2(predacc, genes, featureImp, gseaScore, level1Score, outdir):
    # predacc - prediction accuracy from random forest
    # genes - list of gene sets
    # dirs - working directory
    # setscores - the file of set scores matrix
    # setsamples - the sample names to write into the set score matrix.
    # featureImp
    # gseaScore - scores from ssGSEA
    # level1Score - scores from 1 level

    fout = open(outdir+'analyout.tsv','w')

    #fout.write("set\tmsgs\t1level\tgsea\tfeatimp\tngenes\tgenes\n")
    #print("set\tmsgs\t1level\tgsea\tngenes\n")

    # put the sets in order of prediction ability
    predidx = np.argsort(featureImp)
    i = len(predidx) - 1
    seti = str(i)
    a = str(predacc)
    b = str(level1Score)
    d = str(gseaScore)
    g = str(featureImp[i])
    fout.write('\t'.join([seti,a,b,d,g])+'\n')
    #print('\t'.join([seti,a,b,d,g]))
    fout.close()
    return(1)



# trees, (in, out) where those are pointers to the denovo_trees file.
def analysis3(predacc, genes, featureImp, gseaScore, outdir):
    # predacc - prediction accuracy from random forest
    # genes - list of gene sets
    # dirs - working directory
    # setscores - the file of set scores matrix
    # setsamples - the sample names to write into the set score matrix.
    # featureImp
    # gseaScore - scores from ssGSEA
    # level1Score - scores from 1 level

    fout = open(outdir+'analyout.tsv','w')

    #fout.write("set\tmsgs\tgsea\tfeatimp\tngenes\tgenes\n")
    #print("set\tmsgs\tgsea\tngenes\n")

    # put the sets in order of prediction ability
    predidx = np.argsort(featureImp)

    for i in predidx:
        seti = str(i)
        a = str(predacc)
        d = str(gseaScore)
        c = str(len(genes[i]))
        g = str(featureImp[i])
        fout.write('\t'.join([seti,a,d,c,g])+'\n')
        print('\t'.join([seti,a,d,c,g]))
    fout.close()

    return(1)
