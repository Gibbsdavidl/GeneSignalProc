


#############
  # MAIN #
#############

import sys, getopt
import makeGraphs as mgs
import standardScoring as std

def argProc(args, opts):
    # process options
    mode = ''    # mode is standard or denovo, etc
    datadir = '' # the working directory
    subgraphs = ''  # the subgraphs file
    genesets = ''
    threshold = 0
    numSubGraphs = 200
    maxSubGraphSize = 201
    Nf = 10      # the number of scale-levels.
    filterType = 'heat'
    numCores = 2
    genefile=''
    subgraphfile = ''
    adjfile = ''

    for o, a in opts:
        if o in ("-h", "--help"):
            print("For help use --help")
            print("Modes available: makegraphs, setscoring")
            print("makegraphs args: ")
            print("    -m Mode")
            print("    -d data dir")
            print("    -t threshold for geneset overlaps")
            print("    -c number of cores")
            print("    -n number of subgraphs")
            print("    -x max size of subgraphs")
            print("    -a adjacency file if available")
            print("    -e gene list file if available")
            print("setscoring args: ")
            print("    -m Mode")
            print("    -d data dir")
            print("    -l number of scale levels")
            print("    -f filter name")
            print("    -s subgraph file")
            print("    -c num cores")
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
        elif opt in ('-n'):
            numSubGraphs = arg
        elif opt in ('-x'):
            maxSubGraphSize = arg
        elif opt in ('-s'):
            subgraphs = arg
        elif opt in ('-f'):
            filterType = arg
        elif opt in ('-c'):
            numCores = arg
        elif opt in ('-g'):
            genesets = arg
        elif opt in ('-t'):
            threshold = arg
        elif opt in ('-l'):
            Nf = arg
        elif opt in ('-e'):
            genefile = arg
        elif opt in ('-a'):
            adjfile = arg

    return(mode,
           datadir,
           numSubGraphs,
           maxSubGraphSize,
           subgraphs,
           filterType,
           numCores,
           genesets,
           threshold,
           Nf,
           genefile,
           subgraphfile,
           adjfile
           )


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:d:n:s:f:c:t:x:a:e:l:", ["help"])
    except:
        print("For help use --help")
        print("Modes available: makegraphs, setscoring")
        print("makegraphs args: ")
        print("    -m Mode")
        print("    -d data dir")
        print("    -t threshold for geneset overlaps")
        print("    -c number of cores")
        print("    -n number of subgraphs")
        print("    -x max size of subgraphs")
        print("    -a adjacency file if available")
        print("    -e gene list file if available")
        print("setscoring args: ")
        print("    -m Mode")
        print("    -d data dir")
        print("    -l number of scale levels")
        print("    -f filter name")
        print("    -s subgraph file")
        print("    -c num cores")
        sys.exit(2)

    (mode,datadir,numSubGraphs,maxSubGraphSize,subgraphs,filterType,numCores,genesets,threshold,Nf,genefile,subgraphfile,adjfile) = argProc(args,opts)

    if mode == 'makegraphs':
        mgs.makeGraphs(datadir, numSubGraphs, maxSubGraphSize, genesets, threshold, numCores, adjfile, genefile)

    elif mode == 'setscoring':
        std.runStandard(datadir, Nf, subgraphs, filterType, numCores, genesets)

    else:
        print("Modes available: makegraphs, setscoring")
        print("-d data dir  -m Mode  -f filter name  -n number of scale-levels -s subgraph file")
    return(1)


if __name__ == "__main__":
    sys.exit(main())

