
# want to find sets that are connected across scales and have
# shared high or low filtered values.
#
# read each filtered file, from filtered file list
# get the top and bottom ranked members from each file and scale
# for each top/bot ranked member
#    fit a 2D haar wavelet, maximize the scale on the graph.
#    built a list of sets for each scale
#    collapse sets if there's overlap
# then across scales, but within files
#    if the discovered sets are found across scales
#       then intersect the sets, to connect them across scales.
# write out the multi-scale sets by file/time point.
# done!
#
# the output should look like, a table:
# File   GeneSetID   NodesInGeneSet   ExistInScales   Values
# 1      27          5,7,19           2,3,4           4.3,7.1,6.5

# File, index to the filelist.txt file
# GeneSetID, the denovo gene set found
# NodesInGeneSet, the list of gene indices found in this gene set, aligns to gene order
# ExistInScales, the list of what scales it was found in.
# Values, mean value in each scale.

import sys, getopt
import numpy as np
import igraph as ig
import graphFun

from datetime import datetime, timedelta

print("denovo starting at:")
started = datetime.now()
print(started)

try:
    opts, args = getopt.getopt(sys.argv[1:],"hd:i:o:g:m:")
except getopt.GetoptError:
    print ('denovoGeneSets.py -d <working dir> -i <filelist> -o <output_prefix> -m <score matrix> -g <gene set name>')
    sys.exit(2)
if len(opts) == 0:
    print ('denovoGeneSets.py -d <working dir> -i <filelist> -o <output_prefix> -m <score matrix> -g <gene set name>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('denovoGeneSets.py -d <working dir> -i <filelist> -o <output_prefix> -m <score matrix> -g <gene set name>')
        sys.exit()
    elif opt in ("-i"):
        filelist = arg
    elif opt in ("-d"):
        dirs = arg
    elif opt in ("-o"):
        outputprefix = arg
    elif opt in ("-g"):
        setname = arg
    elif opt in ("-m"):
        adjmat = arg


# get the input files, and where we will write the output file names
inputs = open(dirs+filelist,'r').read().strip().split("\n")
output = open(dirs+outputprefix+'.txt', 'w')

# build the graph
mat = np.loadtxt(dirs+adjmat, delimiter='\t')
#net = ig.Graph.Weighted_Adjacency(mat.tolist(), mode="undirected")
net = ig.Graph.Adjacency(mat.tolist(), mode="undirected")
net.es['weight'] = mat[mat.nonzero()]

net.write_edgelist("edges.txt")
#net = ig.Read_Adjacency(dirs+adjmat)

# then we create a list of the filtered expression matrices
filteredList = []
setList = []
allResults = []

for i in range(0,len(inputs)):
    print("segmenting file " + str(i))
    msr = np.loadtxt(dirs + inputs[i], delimiter='\t')
    setTupleList = graphFun.segmentSpace(net, 0.01, msr)
    filteredList.append(msr)      # the list of multi-scale-signals
    setList.append(setTupleList)  # each input gets a list of set-tuples

#numScales = len((filteredList[0])[:,0])

# want output to be sets that overlap across scales, for each file.
for i in range(0,len(inputs)):
    print("joining sets into groups " + str(i))
    thisSetList = setList[i]
    # sets are connected if they are on the same level, or adjacent levels
    # and they have at least two nodes in common.
    setConnect = graphFun.connectSets(thisSetList)
    # now sets can be coelesed
    setGroups = graphFun.groupCoupling(setConnect, thisSetList)
    # groups come back as a list of sets of tuples (scale-level, set ID)
    resultsList = graphFun.compileResults(i, thisSetList, setGroups, filteredList[i])
    allResults.append(resultsList)

#     [filenum, chainID, thisTuple[1], thisTuple[0], geneID, msr[ti[0]][gi], setMean]

output.write("TimePt\tChainID\tSetID\tLevel\tGeneID\tFiltered\tMeanValue\n")
for x in allResults:
    for y in x:
        z = map(str, y)
        output.write('\t'.join(z)+'\n')
output.close()

finished = datetime.now()
print(finished)
print("done at")
print(finished)
print("took:")
print(finished-started)
