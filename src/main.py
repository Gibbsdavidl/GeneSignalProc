

# the whole enchilada

#################
# DENOVO-PORTION#
#################

# simulate data

   # output expression in time steps, file per time step

# run denovo gene sets

   # outputs a list of chains.

# run denovo extractor

   # output a filtered expr for each chain passing a filter

# run statistical modeller.

   # output model results... predicted output, and summary


#############
# Named Set #
#############

# run denovo extractor

   # output a filtered expr for each chain passing a filter

# run statistical modeller.

   # output model results... predicted output, and summary


import sys, getopt
import simGroups as ss
import filterSignal as fs
import denovoGroupGeneSets as dg
import chainFilterAndExtration as cf

def main():

    # defaults
    Nf = 10
    ngenes = 100
    nparts=5
    ntime=10
    filteredPrefix = "filtered_"
    denovoPrefix = "denovo_"
    levelThresh = 3

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
    x = ss.runSim(datadir, ngenes=ngenes, nparts=nparts, nsamples=ntime)

    # filter the data
    y = fs.filterData(exprfile=x[2], dirs=x[0], outputprefix=filteredPrefix, Nf=Nf, adjmat=x[1])

    # recover the chains
    z = dg.denovoGeneSets(filelist=y[0], dirs=x[0], outputprefix=denovoPrefix, adjmat=x[1])

    # filter chains and extract data for modeling
    w = cf.chainFilterAndEx(dirs=x[1], chainfile=z[0], filterfiles=y[0], levelThresh=levelThresh)


    # run models

    return(1)

if __name__ == "__main__":
    sys.exit(main())

