

import sys, getopt
import graphFun
import waveletFun
import numpy as np
import pygsp as gs
import igraph as ig
from datetime import datetime, timedelta

print("starting at:")
started = datetime.now()
print(started)

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:n:m:d:")
except getopt.GetoptError:
    print ('filterSignal.py -m <adj matrix> -d <working dir> -i <filelist> -o <output_prefix> -n <number of scales>')
    sys.exit(2)
if len(opts) == 0:
    print ('filterSignal.py -m <adj matrix> -d <working dir> -i <filelist> -o <output_prefix> -n <number of scales>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('filterSignal.py -m <adj matrix> -d <working dir> -i <filelist> -o <output_prefix> -n <number of scales>')
        sys.exit()
    elif opt in ("-i"):
        filelist = arg
    elif opt in ("-d"):
        dirs = arg
    elif opt in ("-o"):
        outputprefix = arg
    elif opt in ("-n"):
        Nf = int(arg)
    elif opt in ("-m"):
        adjmat = arg

# made the network in pygsp
print("loading network")
mat = np.loadtxt(dirs+adjmat, delimiter='\t')
gra = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")
net = gs.graphs.Graph(W=mat)
net.directed = False

# for each input file
inputs = open(dirs+filelist,'r').read().strip().split("\n")
outputlist = open(dirs+'filtered_files_list.txt', 'w')

for i in range(0,len(inputs)):
    # compute wavelets
    print("working on file " + str(i))
    # process the data
    sig = graphFun.loadSignal(dirs+inputs[i], 0, 1)  # log10 of value + 0.001
    genes = graphFun.loadGenes(dirs+inputs[i], 0)
    msr = waveletFun.waveletFilter(net, sig, Nf)
    # the filtered signal is in shape (Nf, num_nodes)
    np.savetxt(dirs+outputprefix+"_"+str(i)+".txt", msr, delimiter='\t')
    outputlist.write(outputprefix+'_'+str(i)+'.txt'+'\n')
    print("******************************************")

outputlist.close()
print("finished at:")
stopped = datetime.now()
print(stopped)
print(stopped - started)
