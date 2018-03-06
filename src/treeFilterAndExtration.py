

# filter the chains
# extract the trees across samples
# for modeling.

import csv


def filterTrees(sampleList, treeIDList, levelList, geneIDList, filteredList, levelThresh, topNTrees):
    # want trees, (defined by sampleID, treeID)
    # for a tree ID, need to have x-levels to levelThresh
    trees = []
    genes = []
    maxabs = []
    ina = 0
    outa = 1
    for i in range(2,len(sampleList)):
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
            ina = i
            outa = i+1

    # then take the tree with max abs values greater than the cutoff.
    maxabs.sort(reverse=True)
    cutval = maxabs[topNTrees]
    tree2 = []
    gene2 = []
    for i in range(0,len(trees)):
        if maxabs[i] >= cutval:
            tree2.append(trees[i])
            gene2.append(genes[i])
    return(tree2, gene2)


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

    trees, genes = filterTrees(sampleList, treeIDList, levelList, geneIDList, filteredList, levelThresh, topNTrees)

    return(trees, genes)