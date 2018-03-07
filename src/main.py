

# the whole enchilada

#################
# DENOVO-PORTION#
#################

# simulate data

   # output expression in time steps, file per time step

# run denovo gene sets

   # outputs a list of trees.

# run denovo extractor

   # output a filtered expr for each tree passing a filter

# run statistical modeller.

   # output model results... predicted output, and summary


#############
# Named Set #
#############

# run denovo extractor

   # output a filtered expr for each tree passing a filter

# run statistical modeller.

   # output model results... predicted output, and summary


import sys, getopt
import simGroups as ss
import filterSignal as fs
#import denovoGroupGeneSets as dg
import denovoGeneSets as dg
import treeFilterAndExtration as cf
import model as mm

def main():

    # defaults
    Nf = 10
    ngenes = 100
    nparts=5
    nsamples=20
    filteredPrefix = "filtered_"
    denovoPrefix = "denovo_"
    levelThresh = 3
    topNTrees = 10

    # parse command line options
    mode = 'x'
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:d:", ["help"])
    except:
        print("for help use --help")
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print("help message here")
            sys.exit(0)
    # process arguments
    for opt, arg in opts:
        if opt == '-h':
            print('main.py ')
            sys.exit()
        elif opt in ("-m"):
            mode = arg
        elif opt in ("-d"):
            datadir = arg
        elif opt in ('-nf'):
            Nf = arg

    print('\nworking in ' + datadir)

    # first simulate the data
    x = ss.runSim(datadir, ngenes=ngenes, nparts=nparts, nsamples=nsamples)

    # filter the data
    y = fs.filterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])

    # recover the trees
    z = dg.denovoGeneSets(filelist=y[0], dirs=x[0], outputprefix=denovoPrefix, adjmat=x[1])

    # filter trees and extract data for modeling
    trees, genes = cf.treeFilterAndEx(dirs=x[1], treefile=z[0], filterfiles=y[0], levelThresh=levelThresh, topNTrees=topNTrees)

    # run models
    m = mm.rfModel(dirs=x[0], exprfile=x[2], pheno=x[3], genes=genes)

    # compare model results to simulation.

    # print out comparison and results.

    return(1)

if __name__ == "__main__":
    sys.exit(main())

