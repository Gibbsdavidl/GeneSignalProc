

# filter the chains
# extract the chains across pts
# for modeling.

import csv

def chainFilterAndEx(dirs, filterfiles, chainfile, levelThresh):

    #input = open(chainfile,'r')

    timePtList = []
    chainIDList = []
    levelList = []
    geneIDList = []
    filteredList = []

    with open(chainfile, 'r') as f:
        next(f)  # skip headings
        reader = csv.reader(f, delimiter='\t')
        for timePt, chainID, level, geneID, filtered in reader:
            timePtList.append(timePt)
            chainIDList.append(chainID)
            levelList.append(level)
            geneIDList.append(geneID)
            filteredList.append(filtered)



    return(1)