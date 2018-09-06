


#############
  # MAIN #
#############

import sys, getopt
import makeGraphs as mgs
import standardScoring as std
import os

def argProc(args, opts):
    # process options
    mode = ''    # mode is standard or denovo, etc
    datadir = '' # the working directory
    outputdir = '' # the output directory
    subgraphs = ''  # the subgraphs file
    genesets = ''
    threshold = 0.5
    numSubGraphs = 5
    maxSubGraphSize = 5
    Nf = 10      # the number of scale-levels.
    filterType = 'heat'
    numCores = 2
    genefile=''
    subgraphfile = ''
    adjfile = ''
    exprfile = ''
    phenofile = ''

    for o, a in opts:
        if o in ("-h", "--help"):
            print("For help use --help")
            print("Modes available: makegraphs, setscoring")
            print("makegraphs args: ")
            print("    -m Mode")
            print("    -d data dir")
            print("    -o output dir")
            print("    -t threshold for geneset overlaps")
            print("    -c number of cores")
            print("    -n number of subgraphs")
            print("    -x max size of subgraphs")
            print("    -g genesets file as .gmt")
            print("    -a adjacency file if available")
            print("    -e gene list file if available")
            print("setscoring args: ")
            print("    -m Mode")
            print("    -d data dir")
            print("    -o output dir")
            print("    -r expression data, samples in rows, genes in columns, first col is sample name, first row is gene names")
            print("    -p phenotype file if available, same order as expression data")
            print("    -a adjacency file")
            print("    -s subgraph file")
            print("    -t threshold for geneset overlaps")
            print("    -x threshold for graph edges")
            print("    -e gene list")
            print("    -g genesets file as .gmt")
            print("    -f filter name")
            print("    -l number of scale levels")
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
        elif opt in ('-r'):
            exprfile = arg
        elif opt in ('-p'):
            phenofile = arg
        elif opt in ('-o'):
            outputdir = arg

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
           adjfile,
           exprfile,
           phenofile,
           outputdir
           )


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:d:n:s:f:c:t:x:a:e:l:g:r:p:o:", ["help"])
    except:
        print("For help use --help")
        print("Modes available: makegraphs, setscoring")
        print("makegraphs args: ")
        print("    -m Mode")
        print("    -d data dir")
        print("    -o output dir")
        print("    -g genesets file as .gmt")
        print("    -a adjacency file if available")
        print("    -e gene list file if available")
        print("    -n number of subgraphs")
        print("    -x max size of subgraphs")
        print("    -t threshold for geneset overlaps")
        print("    -c number of cores")
        print("setscoring args: ")
        print("    -m Mode")
        print("    -d data dir")
        print("    -o output dir")
        print("    -r expression data, samples in rows, genes in columns, first col is sample name, first row is gene names")
        print("    -p phenotype file if available, same order as expression data")
        print("    -a adjacency file")
        print("    -s subgraph file")
        print("    -t threshold for geneset overlaps")
        print("    -x threshold for graph edges")
        print("    -e gene list")
        print("    -g genesets file as .gmt")
        print("    -f filter name")
        print("    -l number of scale levels")
        print("    -c num cores")
        sys.exit(2)

    (mode,datadir,numSubGraphs,maxSubGraphSize,subgraphs,filterType,numCores,genesets,threshold,Nf,genefile,subgraphfile,adjfile,exprfile,phenofile,outputdir) = argProc(args,opts)

    if outputdir == '':
        outputdir = datadir

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    if mode == 'makegraphs':
        mgs.makeGraphs(datadir, numSubGraphs, maxSubGraphSize, genesets, threshold, numCores, adjfile, genefile, outputdir)

    elif mode == 'setscoring':
        std.runStandard(datadir, Nf, exprfile, filterType, numCores, subgraphs, genefile, genesets, adjfile, phenofile, threshold, maxSubGraphSize, outputdir)

    else:
        print("Modes available: makegraphs, setscoring")

    print("done")
    return(1)


if __name__ == "__main__":
    sys.exit(main())

