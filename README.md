# GeneSignalProc
Application of graph signal processing to gene expression data.

1. Filter the gene expression data.
   a.) use a list of file names to process each one.
   b.) each file will have an output, with each scale on a line

2. Search for trees that connect segments in scale-levels.
   a.) use connected component labeling in each scale-level (sets)
   b.) connect sets across scale-levels to build trees
   c.) filter trees that lack depth

3. Compare the trees
   a.) based on predictive capability using RandomForests
   b.) rank the trees based on cross validation

To Run:

python3.6 GeneSignalProc/src/main.py -m sim -d /Users/davidgibbs/Data/SimWaveDat/Go6/
