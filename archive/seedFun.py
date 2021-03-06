

# the big dive
import graphFun
import waveletFun
import numpy as np
import pygsp as gs
import igraph as ig

# made the network in pygsp
print("loading network")
mat = np.loadtxt("test/kuldist_selgene_75perc_filter.tsv", delimiter='\t')
gra = ig.Graph.Weighted_Adjacency(list(mat), mode="undirected")
net = gs.graphs.Graph(W=mat)
net.directed = False

# for each input file
dirs = "test/data_tables/"
inputs = [
"GSM2262836_gene_level.tsv",  "GSM2262841_gene_level.tsv",  "GSM2262846_gene_level.tsv",
"GSM2262837_gene_level.tsv",  "GSM2262842_gene_level.tsv",  "GSM2262847_gene_level.tsv",
"GSM2262838_gene_level.tsv",  "GSM2262843_gene_level.tsv",  "GSM2262890_gene_level.tsv",
"GSM2262839_gene_level.tsv",  "GSM2262844_gene_level.tsv",  "GSM2262891_gene_level.tsv",
"GSM2262840_gene_level.tsv",  "GSM2262845_gene_level.tsv"
]

filelabel = ["1h", "d1", "d6", "1h", "d1", "d6", "4h", "d6", "0h", "4h", "d6", "0h", "4h", "d6"]
fileorder = [3,     8,     9,    2,    7,   13,   4,    10,   1,    5,    11,   0,    6,    12]

Nf = 10
sm_eps = 0.05
t_eps = 1.0
#seed = 37  # CD9
#seed = 1188 # TLR1
seed = 702 # TLR4
#seed = 1292 # LAMP1
#seed = 542 # EGR2


ls = [] # labels
os = [] # orders
cc = [] # coefficients
sc = [] # scales
rs = [] # ranks
aw = [] # raw signal

for i in range(0,len(inputs)):
    # compute wavelets
    print("working on file " + str(i))
    sig = graphFun.loadSignal(dirs+inputs[i])  # log10 of value + 0.001
    genes = graphFun.loadGenes(dirs+inputs[i])
    msr = waveletFun.waveletFilter(net, sig, Nf)
    # the filtered signal is in shape (Nf, num_nodes)
    # for given seed, compute segmentations.
    msr = graphFun.seedSpace(gra, sm_eps, t_eps, seed, msr, sig)
    # collect the results for this input
    ls += [filelabel[i] for x in range(0,Nf)]
    os += [fileorder[i] for x in range(0,Nf)]
    cc += [x[0] for x in msr[1]]
    rs += [x[0] for x in msr[2]]
    sc += [i for i in range(0,Nf)]
    aw += [x[0] for x in msr[3]]
    print("******************************************")

fout = open("lamp1_seed_out.txt", 'w')
fout.write('\t'.join(ls) + '\n')
fout.write('\t'.join([str(oi) for oi in os]) + '\n')
fout.write('\t'.join([str(x) for x in cc]) + '\n')
fout.write('\t'.join([str(x) for x in sc]) + '\n')
fout.write('\t'.join([str(x) for x in rs]) + '\n')
fout.write('\t'.join([str(x) for x in aw]) + '\n')
fout.close()
# compare segmentation levels across time points.
