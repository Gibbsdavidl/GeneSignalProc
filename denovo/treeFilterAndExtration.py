

# filter the chains
# extract the trees across samples
# for modeling.

import csv
import numpy as np
import copy

def filterTrees(sampleList, treeIDList, levelList, geneIDList, filteredList, levelThresh, topNTrees):
    # want trees, (defined by sampleID, treeID)
    # for a tree ID, need to have x-levels to levelThresh
    trees = []
    genes = []
    maxabs = []
    minabs = []
    means = []
    meds = []
    levs = []
    nlev = []
    ina = 0
    outa = 1
    for i in range(2,len(sampleList)): # for each row in the tree table.
        if sampleList[outa] == sampleList[i] and treeIDList[outa] ==  treeIDList[i]:
            # keep expanding the tree
            outa = i
        else:
            # examine tree and start a new tree
            levelRange = levelList[ina:(outa+1)]
            if len(set(levelRange)) >= levelThresh:
                # keep tree
                trees.append( (ina, (outa+1)) )
                genes.append( set(geneIDList[ina:outa+1]) )
                maxabs.append( max([abs(x) for x in filteredList[ina:outa+1]]) )
                minabs.append(min([x for x in filteredList[ina:outa + 1]]))
                means.append(np.mean([abs(x) for x in filteredList[ina:outa+1]]) )
                meds.append(np.median([abs(x) for x in filteredList[ina:outa+1]]) )
                levs.append(set(levelList[ina:outa+1]))
                nlev.append(len(set(levelList[ina:outa+1])))
            ina = i
            outa = i+1

    # then take the tree with max abs values greater than the cutoff.
    #maxabs.sort(reverse=True)
    #cutval = maxabs[topNTrees]

    tree2 = []
    gene2 = []
    vals2 = []
    levs2 = []
    nlev2 = []

    med2  = copy.deepcopy(meds)

    med2.sort(reverse=True)
    cutval = med2[int(topNTrees/2)]
    for i in range(0,len(trees)):
        if med2[i] >= cutval:
            tree2.append(trees[i])
            gene2.append(genes[i])
            vals2.append(meds[i])
            levs2.append(levs[i])
            nlev2.append(nlev[i])

    med2.sort(reverse=False)
    cutval = med2[int(topNTrees/2)]
    for i in range(0,len(trees)):
        if meds[i] <= cutval:
            tree2.append(trees[i])
            gene2.append(genes[i])
            vals2.append(meds[i])
            levs2.append(levs[i])
            nlev2.append(nlev[i])

    # want the tree with the most levels.

    return(tree2, gene2, vals2, levs2, nlev2)


def treeFilterAndEx(dirs, filterfiles, treefile, levelThresh, topNTrees):

    #input = open(treefile,'r')

    sampleList = []
    treeIDList = []
    levelList = []
    geneIDList = []
    filteredList = []

    with open(treefile, 'r') as f:
        next(f)  # skip headings
        reader = csv.reader(f, delimiter='\t')
        for sampleID, treeID, level, geneID, filtered in reader:
            sampleList.append(sampleID)
            treeIDList.append(treeID)
            levelList.append(int(level))
            geneIDList.append(int(geneID))
            filteredList.append(float(filtered))

    trees, genes, vals, levs, nlev = filterTrees(sampleList, treeIDList, levelList, geneIDList, filteredList, levelThresh, topNTrees)

    return(trees, genes, vals, levs, nlev)