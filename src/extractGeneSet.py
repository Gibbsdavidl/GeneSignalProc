

# then take gene_sets_matrix, and a gene set name,
# read each filtered file, from filtered file list
# cluster one of them (first?) use that gene ordering.
# subset the filtered values from each file,
# write a new file with that gene ordering.
# NICE to have a heatmap, using this ordering, for each sample.

import sys, getopt

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
outputlist = open(dirs+'extracted_files_list.txt', 'w')

# read in the set membership matrix
setmat = np.loadtxt(dirs+setmatname, delimiter='\t', dtype='U128')

# then we create a list of the filtered expression matrices
filteredList = []
for i in range(0,len(inputs)):
    mat = np.loadtxt(dirs + inputs[i], delimiter='\t')
    filteredList.append(mat)

# then we want the genes in the gene set specified.

# create a list of sub-matrices.

# cluster one of them, and get the ordering.

# reorder each sub-matrix

# write out each sub-matrix

# make a cool viz of each one.



print("done")