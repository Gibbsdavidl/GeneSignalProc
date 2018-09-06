

# make a multi scale representation of the data via filtering

import sys, getopt
import graphFun
import waveletFun
import numpy as np
import pygsp as gs
import igraph as ig
import os
import gzip
import scipy
from datetime import datetime, timedelta


def formatExprData(dirs,exprfile,allgenes):
    reformat = []
    samples  = []
    graphGenes = [gi.decode('utf-8') for gi in allgenes]
    inputs = gzip.open(dirs+exprfile,'rt').read().strip().split("\n")
    header = inputs[0].strip().split('\t')
    genedict = {g:(i-1) for i,g in enumerate(header) if i > 0}
    for i in range(1,len(inputs)):
        thisline = []
        bits = inputs[i].strip().split('\t')
        vals = [float(x) for i,x in enumerate(bits) if i > 0]
        samples.append(bits[0])
        for gi in graphGenes:
            if gi in genedict:
                thisline.append(vals[genedict[gi]])
            else:
                #print('Gene missing in data:' + gi)
                thisline.append(0.0)
        reformat.append(thisline)
    return( (reformat, samples) )




def noFilterData(exprfile, dirs, outputprefix, Nf, adjmat, allgenes, outdir):
    # exprfile - matrix of gene expression
    # dirs - the working directory
    # outputprefix - the prefix put on filter file outputs
    # Nf - number of scale levels
    # adjmat - the file name for the adjacent matrix

    if not os.path.exists(outdir+'filtered_files'):
        os.makedirs(outdir+'filtered_files')

    # for each input file
    sigs, samps = formatExprData(dirs,exprfile, allgenes)     #open(dirs+exprfile,'r').read().strip().split("\n")  ###########!!!!!!!!!!!!!! NEW FORMAT!!!!!!!!!!!!!!!!!
    outputlist = open(outdir+'not_filtered_files_list.txt', 'w')
    samplelist = open(outdir+'not_filtered_sample_list.txt', 'w')

    filteredSignal = []

    # compute wavelets
    # process the data
    print("filtering data")
    for i in range(0,len(sigs)):
        #vals = inputs[i].split('\t')
        sig = np.array([sigs[i]]) # np.array([float(x) for x in vals[1:len(vals)]])
        np.savetxt(outdir+'filtered_files/'+outputprefix+str(i)+".txt", sig, delimiter='\t')
        outputlist.write('filtered_files/'+outputprefix+str(i)+'.txt'+'\n')
        samplelist.write(samps[i] + '\n')

    outputlist.close()
    samplelist.close()
    print("finished not_filtered_files_list data")
    return(['not_filtered_files_list.txt'])



def heatFilterData(exprfile, dirs, outputprefix, Nf, adjmat, allgenes, edgeT, outdir):
    # exprfile - matrix of gene expression
    # dirs - the working directory
    # outputprefix - the prefix put on filter file outputs
    # Nf - number of scale levels
    # adjmat - the file name for the adjacent matrix

    if not os.path.exists(outdir+'filtered_files'):
        os.makedirs(outdir+'filtered_files')

    # made the network in pygsp
    print("loading network")
    mat = scipy.sparse.load_npz(dirs+adjmat)
    mat = mat.multiply(mat >= edgeT)
    #gra = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")
    net = gs.graphs.Graph(W=mat)
    net.directed = False
    net.estimate_lmax()
    #net.compute_fourier_basis()

    # for each input file
    sigs, samps = formatExprData(dirs,exprfile, allgenes)     ### needs to be same order as graph ###
    outputlist = open(outdir+'filtered_files_list.txt', 'w')
    samplelist = open(outdir+'filtered_sample_list.txt', 'w')

    # compute wavelets
    # process the data
    print("filtering data")
    for i in range(0,len(sigs)):
        sig = np.array(sigs[i])
        msr = waveletFun.heatFilter(net, sig, Nf)                                                  # list of filtered signal for each sample
        np.savetxt(outdir+'filtered_files/'+outputprefix+str(i)+".txt", msr, delimiter='\t')         # the filtered signal is in shape (Nf, num_nodes)
        outputlist.write('filtered_files/'+outputprefix+str(i)+'.txt'+'\n')
        samplelist.write(samps[i] + '\n')

    outputlist.close()
    samplelist.close()
    print("finished filtering data")
    return(['filtered_files_list.txt'])


def mexFilterData(exprfile, dirs, outputprefix, Nf, adjmat):
    # exprfile - matrix of gene expression
    # dirs - the working directory
    # outputprefix - the prefix put on filter file outputs
    # Nf - number of scale levels
    # adjmat - the file name for the adjacent matrix

    if not os.path.exists('filtered_files'):
        os.makedirs('filtered_files')

    # made the network in pygsp
    print("loading network")
    mat = np.loadtxt(dirs+adjmat, delimiter='\t')
    gra = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")
    net = gs.graphs.Graph(W=mat)
    net.directed = False
    net.estimate_lmax()
    net.compute_fourier_basis()

    # for each input file
    inputs = open(dirs+exprfile,'r').read().strip().split("\n")  ###########!!!!!!!!!!!!!! NEW FORMAT!!!!!!!!!!!!!!!!!
    outputlist = open(dirs+'filtered_files_list.txt', 'w')
    samplelist = open(dirs+'filtered_sample_list.txt', 'w')

    filteredSignal = []

    # compute wavelets
    # process the data
    print("filtering data")
    for i in range(1,len(inputs)):  # skip header.
        vals = inputs[i].split('\t')
        sig = np.array([float(x) for x in vals[1:len(vals)]])
        msr = waveletFun.waveletFilter(net, sig, Nf) # list of filtered signal for each sample
        # the filtered signal is in shape (Nf, num_nodes)
        np.savetxt(dirs+'filtered_files/'+outputprefix+str(i)+".txt", msr, delimiter='\t')
        outputlist.write('filtered_files/'+outputprefix+str(i)+'.txt'+'\n')
        samplelist.write(vals[0] + '\n')

    outputlist.close()
    samplelist.close()
    print("finished filtering data")
    return(['filtered_files_list.txt'])

