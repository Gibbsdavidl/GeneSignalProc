

# make a multi scale representation of the data via filtering

import sys, getopt
import graphFun
import waveletFun
import numpy as np
import pygsp as gs
import igraph as ig
from datetime import datetime, timedelta

def filterData(exprfile, dirs, outputprefix, Nf, adjmat):

    #print ('filterSignal.py -m <adj matrix> -d <working dir> -i <filelist> -o <output_prefix> -n <number of scales>')
    # made the network in pygsp
    print("loading network")
    mat = np.loadtxt(dirs+adjmat, delimiter='\t')
    gra = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")
    net = gs.graphs.Graph(W=mat)
    net.directed = False

    # for each input file
    inputs = open(dirs+exprfile,'r').read().strip().split("\n")  ###########!!!!!!!!!!!!!! NEW FORMAT!!!!!!!!!!!!!!!!!
    outputlist = open(dirs+'filtered_files_list.txt', 'w')
    samplelist = open(dirs+'filtered_sample_list.txt', 'w')

    filteredSignal = []

    for i in range(1,len(inputs)):  # skip header.
        # compute wavelets
        print("working on file " + str(i))
        # process the data
        vals = inputs[i].split('\t')
        sig = np.array([float(x) for x in vals[1:len(vals)]])
        msr = waveletFun.waveletFilter(net, sig, Nf) # list of filtered signal for each sample
        # the filtered signal is in shape (Nf, num_nodes)
        np.savetxt(dirs+outputprefix+str(i)+".txt", msr, delimiter='\t')
        outputlist.write(outputprefix+str(i)+'.txt'+'\n')
        samplelist.write(vals[0] + '\n')
        print("******************************************")

    outputlist.close()
    samplelist.close()
    print("finished filtering data")
    return(['filtered_files_list.txt'])
