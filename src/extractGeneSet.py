

# then take gene_sets_matrix, and a gene set name,
# read each filtered file, from filtered file list
# cluster one of them (first?) use that gene ordering.
# subset the filtered values from each file,
# write a new file with that gene ordering.
# NICE to have a heatmap, using this ordering, for each sample.

import sys, getopt
#from sklearn import cluster
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np

from datetime import datetime, timedelta

print("starting at:")
started = datetime.now()
print(started)

try:
    opts, args = getopt.getopt(sys.argv[1:],"hd:i:o:g:m:")
except getopt.GetoptError:
    print ('extractGeneSet.py -d <working dir> -i <filelist> -o <output_prefix> -m <gene set matrix> -g <gene set name>')
    sys.exit(2)
if len(opts) == 0:
    print ('extractGeneSet.py -d <working dir> -i <filelist> -o <output_prefix> -m <gene set matrix> -g <gene set name>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('extractGeneSet.py -d <working dir> -i <filelist> -o <output_prefix> -m <gene set matrix> -g <gene set name>')
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
        setmatname = arg


# get the input files, and where we will write the output file names
inputs = open(dirs+filelist,'r').read().strip().split("\n")
output = open(dirs+'extracted.txt', 'w')

# read in the set membership matrix
setmat = np.loadtxt(dirs+setmatname, delimiter='\t', dtype='U128')

# then we create a list of the filtered expression matrices
filteredList = []
for i in range(0,len(inputs)):
    mat = np.loadtxt(dirs + inputs[i], delimiter='\t')
    filteredList.append(mat)

# then we want the genes in the gene set specified.
# this is the set
seti = (np.where([float(setname == setmat[i][0]) for i in range(0,len(setmat))]))[0][0]
# and these are the indices into gene order.  0 here is 1 in the setmatrix.
gidx = np.where( [ '1' == setmat[3][i] for i in range(1,len(setmat[0]))] )[0]
# gene names


# create a list of sub-matrices.
subMatrixList = []
for fi in filteredList:
    submat = np.zeros( (len(fi), len(gidx) ))
    for ri in range(0,len(fi)):  # for each row in the matrix
        submat[ri, :] = fi[ri, gidx] # copy over the values for this gene set
    subMatrixList.append(submat)

# cluster one of them, and get the ordering.
#ward = cluster.ward_tree(subMatrixList[4])
Z = linkage(np.transpose(subMatrixList[4]), 'ward')
Q = dendrogram(Z, get_leaves=True, distance_sort=True)
geneOrder = Q['leaves']

# reorder each sub-matrix
# make a cool viz of each one.


# write out each sub-matrix
# write it in tidy format.
output.write("TimePt\tScale\tGene\tValue\n")
for i,sm in enumerate(subMatrixList):  # i is the time point
    for ri in range(0, len(sm)):           # ri is the scale
        for ci in range(0, len(sm[ri])):   # ci is the filtered value
            output.write(str(i) +'\t' + str(ri) + '\t' + str(gidx[ci]) + '\t' + str(sm[ri,ci]) + '\n')

print("done")
