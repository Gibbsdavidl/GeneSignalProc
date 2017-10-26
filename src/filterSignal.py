

# the big dive
import graphFun
import waveletFun
import numpy as np
import pygsp as gs
import igraph as ig
from datetime import datetime, timedelta

print("starting at:")
started = datetime.now()
print(started)

# made the network in pygsp
print("loading network")
mat = np.loadtxt("data/pruned_filtered_intersected_K2_matrix.tsv.gz", delimiter='\t')
gra = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")
net = gs.graphs.Graph(W=mat)
net.directed = False

# for each input file
dirs = "data/data_tables/"
inputs = open("fileNames.txt",'w').read().strip().split("\n")

#for i in range(0,len(inputs)):
for i in range(0,1):
    # compute wavelets
    print("working on file " + str(i))
    # process the data
    sig = graphFun.loadSignal(dirs+inputs[i])  # log10 of value + 0.001
    genes = graphFun.loadGenes(dirs+inputs[i])
    msr = waveletFun.waveletFilter(net, sig, Nf)
    # the filtered signal is in shape (Nf, num_nodes)
    np.savetxt(fname="filter_out"+inputs[i]+".txt", msr, delimiter='\t')
    print("******************************************")

print("finished at:")
stopped = datetime.now()
print(stopped)
print(stopped - started)
