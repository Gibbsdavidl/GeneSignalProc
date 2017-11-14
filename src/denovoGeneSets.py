
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

print("starting at:")
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
output = open(dirs+'extracted.txt', 'w')

# build the graph
mat = np.loadtxt(dirs+adjmat, delimiter='\t')
gra = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")

# then we create a list of the filtered expression matrices
filteredList = []
for i in range(0,len(inputs)):
    mat = np.loadtxt(dirs + inputs[i], delimiter='\t')
    filteredList.append(mat)


graphFun.segment(gra, 1, 2, 27, filteredList[1], [], 2)
