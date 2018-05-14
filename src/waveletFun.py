
import numpy as np
import pygsp as g


#g <- graph.tree(n=40, mode="undirected")
#mat <- as_adj(g, type="both", names=F, sparse=F)
#write.table(mat, file="test2_graph.txt", sep="\t", row.names=F, col.names=F, quote=F)

# made the network in Gephi

def loadNet(netfile):
    mat = np.loadtxt(netfile, delimiter='\t')
    net = g.graphs.Graph(W=mat, directed=False)
    return(net)


def normi(x):
    return(x / max(abs(x)))


def waveletFilter(net, signal, Nf):
    mexhat = g.filters.MexicanHat(net, Nf)
    sighat = mexhat.analyze(signal)
    sighat_transpose = sighat.transpose()
    return(sighat_transpose)


def heatFilter(net, signal, Nf):
    tis = [ti for ti in range(10, Nf*20, 20)]
    hf = g.filters.Heat(net, tau=tis, normalize=False)
    sighat = hf.analyze(signal)
    sighat_transpose = sighat.transpose()
    return(sighat_transpose)

def imageFilteredSig(X):
    def format_coord(x, y):
        col = int(x + 0.5)
        row = int(y + 50.5)
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            z = X[row, col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f' % (x, y)
    fig, ax = plt.subplots()
    ax.imshow(X, interpolation='nearest')
    numrows, numcols = X.shape
    ax.format_coord = format_coord
    plt.show()
