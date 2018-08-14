

import numpy as np

genesetfile = 'h.all.v6.0.symbols.gmt'

genesets = open(genesetfile, 'r').read().strip().split('\n')  # gene sets

datmat = open('gtexSalivary Gland_rpkm.tsv', 'r').read().strip().split('\n')  # gtex data

datgenes = datmat[0].strip().split('\t')  # first row, gene names in data

d = {g: i for i, g in enumerate(datgenes)}  # index to gene names

for gs in genesets:  # for each gene set, we will make three data sets for each level of boost (neg boost too)

    bits = gs.split('\t')

    # the genes start at index=2 #

    boostVals = [-2.0, 1.5, 2.0]

    for repi in [0,1,2]:

        for bi in boostVals:
            # going to make a set of 64 samples.
            pheno = np.random.choice([0.0, 1.0], len(datmat)-1)

            pfout = open('pheno_' + str(repi) + '_' + str(bi) + '_' + bits[0] + '.tsv', 'w')
            pfout.write('\n'.join([str(ppi) for ppi in pheno])+'\n')
            pfout.close()

            gfout = open('data__' + str(bi) + '_' + str(repi) + '_' + bits[0] + '.tsv', 'w')
            gfout.write('\t'.join(datgenes) + '\n')

            # going to write out the whole matrix
            for i,pi in enumerate(pheno):      # This is the number of rows
                print(i)

                datmatvals = datmat[(i+1)].strip().split('\t')

                rowdat = []

                for j,gj in enumerate(datgenes): # for each row / sample in the matrix

                    if j == 0:  # the sample name
                        rowdat.append(datmatvals[j])
                    elif gj in gs and pi == 1.0:  # this gene in is the gene set
                        rowdat.append(np.log2(float(datmatvals[j])+1.0) + np.random.normal(bi, 0.25, 1)[0])
                    else:
                        rowdat.append(np.log2(float(datmatvals[j])+1.0) + np.random.normal(0.0, 0.25, 1)[0])

                # at end of row, write it out
                gfout.write('\t'.join([str(rri) for rri in rowdat])+'\n')
            gfout.close()







