

# going to use the chains as groups in group lasso with overlaps.


# groups are just the integer that labels the group, with length equal to the number of features.

# too many features??


print("denovo starting at:")
started = datetime.now()
print(started)

try:
    opts, args = getopt.getopt(sys.argv[1:],"hd:i:o:g:m:")
except getopt.GetoptError:
    print ('denovoGroupLasso.py -d <working dir> -i <filename> -o <output_prefix> ')
    sys.exit(2)
if len(opts) == 0:
    print ('denovoGroupLasso.py -d <working dir> -i <filename> -o <output_prefix> ')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('denovoGroupLasso.py -d <working dir> -i <filename> -o <output_prefix> ')
        sys.exit()
    elif opt in ("-i"):
        filein = arg
    elif opt in ("-d"):
        dirs = arg
    elif opt in ("-o"):
        outputprefix = arg

# get the input files, and where we will write the output file names
inputs = open(dirs+filelist,'r').read().strip().split("\n")
output = open(dirs+outputprefix+'.txt', 'w')


