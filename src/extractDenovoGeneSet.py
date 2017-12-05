

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

print("extraction starting at:")
started = datetime.now()
print(started)

setid = '-1'
chain = '-1'
timept = '-1'

try:
    opts, args = getopt.getopt(sys.argv[1:],"hd:i:f:o:t:c:s:r:")
except getopt.GetoptError:
    print ('extractDenovoGeneSet.py -d <working dir> -i <filelist> -f <denovo file> -o <output_prefix> -t <time pt> -c <chain> -s <set ID> -r <level filter>')
    sys.exit(2)
if len(opts) == 0:
    print ('extractDenovoGeneSet.py -d <working dir> -i <filelist> -f <denovo file> -o <output_prefix> -t <time pt> -c <chain> -s <set ID')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('extractDenovoGeneSet.py -d <working dir> -i <filelist> -f <denovo file> -o <output_prefix> -t <time pt> -c <chain> -s <set ID')
        sys.exit()
    elif opt in ("-i"):
        filelist = arg
    elif opt in ("-f"):
        denovofile = arg
    elif opt in ("-d"):
        dirs = arg
    elif opt in ("-o"):
        outputprefix = arg
    elif opt in ("-t"):
        timept = arg
    elif opt in ("-c"):
        chain = arg
    elif opt in ("-s"):
        setid = arg
    elif opt in ("-r"):
        lvlfilt = float(arg)


# get the input files, and where we will write the output file names
inputs = open(dirs+filelist,'r').read().strip().split("\n")
denovos = open(dirs+denovofile,'r').read().strip().split('\n')
output = open(dirs+outputprefix+'.txt', 'w')

# then we create a list of the filtered expression matrices
filteredList = []
for i in range(0,len(inputs)):
    mat = np.loadtxt(dirs + inputs[i], delimiter='\t')
    filteredList.append(mat)

# then we want the genes in the gene set specified.
gidx = list()
for line in denovos[1:]:
    bits = line.strip().split('\t')
    if bits[0] == timept and bits[1] == chain:
        if setid == '-1': # did not set the set ID
            gidx.append(bits[4])
        elif setid == bits[2]:
            gidx.append(bits[4])
        # else don't keep this gene

# take only genes that are found across levels
gidx = np.array(gidx)
hidx = [int(x) for x in gidx if sum(x == gidx) > lvlfilt]
hidx = list(set(hidx))

# create a list of sub-matrices.
subMatrixList = []
for fi in filteredList:
    submat = np.zeros( (len(fi), len(hidx) ))
    for ri in range(0,len(fi)):  # for each row in the matrix
        submat[ri, :] = fi[ri, hidx] # copy over the values for this gene set
    subMatrixList.append(submat)

# cluster one of them, and get the ordering.
#ward = cluster.ward_tree(subMatrixList[4])
#Z = linkage(np.transpose(subMatrixList[4]), 'ward')
#Q = dendrogram(Z, get_leaves=True, distance_sort=True)
#geneOrder = Q['leaves']

# write out each sub-matrix
# write it in tidy format.
output.write("TimePt\tScale\tGene\tValue\n")
for i,sm in enumerate(subMatrixList):  # i is the time point
    for ri in range(0, len(sm)):           # ri is the scale
        for ci in range(0, len(sm[ri])):   # ci is the filtered value
            output.write(str(i) +'\t' + str(ri) + '\t' + str(gidx[ci]) + '\t' + str(sm[ri,ci]) + '\n')

print("done")
