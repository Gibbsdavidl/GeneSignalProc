


#############
  # MAIN #
#############

import sys, getopt
import runDenovoSimulation as rds
import runStandardSimulation as std

def argProc(args, opts):
    # process options
    mode = ''    # mode is standard or denovo, etc
    datadir = '' # the working directory
    subgraphs = ''  # the subgraphs file
    Nf = ''      # the number of scale-levels.

    for o, a in opts:
        if o in ("-h", "--help"):
            print("Modes available: denovo_sim, denovo_sim_rerun")
            print("-m Mode  -d data dir  -nf  number of scale-levels")
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
            Nf = arg
        elif opt in ('-s'):
            subgraphs = arg
    return(mode, datadir, Nf, subgraphs)


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:d:n:s:", ["help"])
    except:
        print("for help use --help")
        print("Modes available: denovo_sim, denovo_sim_reuse_data")
        print("-m Mode  -d data dir  -n number of scale-levels")
        sys.exit(2)

    mode, datadir, Nf, subgraphs = argProc(args,opts)

    if mode == 'denovo_sim':
        rds.runDenovoSim(datadir, Nf, subgraphs)

    elif mode == 'denovo_sim_rerun':
        rds.runDenovoSimRerun(datadir, Nf, subgraphs)

    elif mode == 'standard_sim':
        std.runStandard(datadir, Nf, subgraphs)

    else:
        print("Modes: denovo_sim, denovo_sim_reuse_data ")

    return(1)


if __name__ == "__main__":
    sys.exit(main())

