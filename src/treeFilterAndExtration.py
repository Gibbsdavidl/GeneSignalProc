

# filter the chains
# extract the trees across samples
# for modeling.

import csv

def treeFilterAndEx(dirs, filterfiles, treefile, levelThresh):

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
            geneIDList.append(geneID)
            filteredList.append(float(filtered))



    return(1)