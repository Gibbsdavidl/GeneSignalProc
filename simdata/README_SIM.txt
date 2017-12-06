
# sim

# 1. make some random genes

# 2. create gene sets that have random gene memberships

# 3. remove gene sets if we still have complete coverage over genes

# 4. assign each gene set a parameter set

# 5. generate some value per gene, depending on a sum over distributions, parameterized by set.

# 6. find gene ordering based on overlap of gene sets

# 7. write expression out for gene list, in order.

# 8. filter each file:  GeneSignalProc/src/filterSignal.py -m scorematrix.tsv -d simdata/simsets/ -i filelist.txt -o exprtest1 -n 10
#    >>>filtering, write a list of the filtered files. Each row is a scale.

### need gene set names ###


# then take gene_sets_matrix, and a gene set name,
# read each filtered file, from filtered file list
# cluster one of them (first?) use that gene ordering.
# subset the filtered values from each file,
# write a new file with that gene ordering.
# NICE to have a heatmap, using this ordering, for each sample.


