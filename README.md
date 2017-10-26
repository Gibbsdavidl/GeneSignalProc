# GeneSignalProc
Application of graph signal processing to gene expression data.

1. Filter the gene expression data.
   a.) use a list of file names to process each one.
   b.) each file will have an output, with each scale on a line

2. Search for keypoints and keysets.
   a.) use the file list, again, they will be processed to match the outputs from (1)
   b.) key sets will be written out for each file.

3. Compare the keysets between files,
   a.) based on groups
   b.) based on time series
   c.) based on same set (one keyset defined in one result) across files
   d.) compare overlap of keysets across files

