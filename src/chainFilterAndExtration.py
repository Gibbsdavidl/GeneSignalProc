

# filter the chains
# extract the chains across pts
# for modeling.

import csv

def chainFilterAndEx(dirs, filterfiles, chainfile, levelThresh):

    input = open('dirs'+chainfile)

    timePt = []
    chainID = []
    level = []
    geneID = []
    filtered = []

    with open('data.csv', 'r') as f:
        next(f)  # skip headings
        reader = csv.reader(f, delimiter='\t')
        for timePt, chainID, level, geneID, filtered in reader:
            timePt.append(name)
            chainID.append(age)
            level.append(level)
            geneID.append(geneID)
            filtered.append(filtered)



    return(1)