

# example of the annotation file

#SAMPID	SMATSSCR	SMCENTER	SMPTHNTS	SMRIN	SMTS
#GTEX-1117F-0003-SM-58Q7G		B1			Blood
#GTEX-1117F-0003-SM-5DWSB		B1			Blood
#GTEX-1117F-0226-SM-5GZZ7	0	B1	"2 pieces, ~15% vessel stroma, rep delineated"	6.8	Adipose Tissue
#GTEX-1117F-0426-SM-5EGHI	0	B1	"2 pieces, !5% fibrous connective tissue, delineated (rep)"	7.1	Muscle
#GTEX-1117F-0526-SM-5EGHJ	0	B1	"2 pieces, clean, Monckebeg medial sclerosis, rep delineated"	8	Blood Vessel
#GTEX-1117F-0626-SM-5N9CS	1	B1	"2 pieces, up to 4mm aderent fat/nerve/vessel, delineated"	6.9	Blood Vessel
#GTEX-1117F-0726-SM-5GIEN	1	B1	"2 pieces, no abnormalities"	6.3	Heart
#GTEX-1117F-1326-SM-5EGHH	1	B1	"2 pieces, diffuse mesothelial hyperplasia; ~10% vessel/fibrous tissue (delineated)"	5.9	Adipose Tissue
#GTEX-1117F-2226-SM-5N9CH	1	B1	"1 piece vascular tissue, probably ovarian hilum.  Not GTEx target tissue"	6.6	Ovary



#1.2
#56238	8555
#Name	Description	GTEX-111CU-1826-SM-5GZYN	GTEX-111FC-0226-SM-5N9B8	GTEX-111VG-2326-SM-5N9BK	GTEX-111YS-2426-SM-5GZZQ
#ENSG00000223972.4	DDX11L1	0	0	0	0
#ENSG00000227232.4	WASH7P	6.50895977020264	10.7456922531128	6.67049884796143	6.3844690322876
#ENSG00000243485.2	MIR1302-11	0	0	0	0
#ENSG00000237613.2	FAM138A	0	0	0	0
#ENSG00000268020.2	OR4G4P	0	0	0	0
#ENSG00000240361.1	OR4G11P	0	0	0	0
#ENSG00000186092.4	OR4F5	0	0	0	0


# two imput files

# read annotation file.
# make a dict
# get col 0, and col 5

# then read expr file:
#  skip first two lines
#  then get line 2, list of sample IDs

# for each type of tissue:
    #  new file needs first row SampID, GeneIDs...
    #  each row is a sample, then values.



annotfile = '/Users/davidgibbs/Data/gtex/GTEx_Data_V6_Annotations_SampleAttributesDS.txt'
#exprfile = '/Users/davidgibbs/Data/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct'
exprfile = '/Users/davidgibbs/Data/gtex/gtex_1000.tsv'
genefile = '/Users/davidgibbs/Data/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.genes.txt'

annots = open(annotfile,'r').read().strip().split('\n')
sampDict = dict()
tissueDict = dict()
tissueSet = set()
for i in range(2,len(annots)):
    bits = annots[i].strip().split('\t')
    sampDict[bits[0]] = bits[5]
    tissueSet.add(bits[5])
    if bits[5] not in tissueDict:
        tissueDict[bits[5]] = [bits[0]]
    else: # already have this tissue type
        tissueDict[bits[5]].append(bits[0])


print(len(sampDict))
print(len(tissueDict))
print(len(tissueSet))

# get the gene list
genes = open(genefile, 'r').read().strip().split('\n')
genes = genes[3:]

print(len(genes))
print(genes[0:10])

proc = 3

for tissue in list(tissueSet):
    if proc == 0:
        exit(1)
    else:
        proc -= 1
    print(tissue)
    if tissue != '':

# example of the annotation file

#SAMPID	SMATSSCR	SMCENTER	SMPTHNTS	SMRIN	SMTS
#GTEX-1117F-0003-SM-58Q7G		B1			Blood
#GTEX-1117F-0003-SM-5DWSB		B1			Blood
#GTEX-1117F-0226-SM-5GZZ7	0	B1	"2 pieces, ~15% vessel stroma, rep delineated"	6.8	Adipose Tissue
#GTEX-1117F-0426-SM-5EGHI	0	B1	"2 pieces, !5% fibrous connective tissue, delineated (rep)"	7.1	Muscle
#GTEX-1117F-0526-SM-5EGHJ	0	B1	"2 pieces, clean, Monckebeg medial sclerosis, rep delineated"	8	Blood Vessel
#GTEX-1117F-0626-SM-5N9CS	1	B1	"2 pieces, up to 4mm aderent fat/nerve/vessel, delineated"	6.9	Blood Vessel
#GTEX-1117F-0726-SM-5GIEN	1	B1	"2 pieces, no abnormalities"	6.3	Heart
#GTEX-1117F-1326-SM-5EGHH	1	B1	"2 pieces, diffuse mesothelial hyperplasia; ~10% vessel/fibrous tissue (delineated)"	5.9	Adipose Tissue
#GTEX-1117F-2226-SM-5N9CH	1	B1	"1 piece vascular tissue, probably ovarian hilum.  Not GTEx target tissue"	6.6	Ovary



#1.2
#56238	8555
#Name	Description	GTEX-111CU-1826-SM-5GZYN	GTEX-111FC-0226-SM-5N9B8	GTEX-111VG-2326-SM-5N9BK	GTEX-111YS-2426-SM-5GZZQ
#ENSG00000223972.4	DDX11L1	0	0	0	0
#ENSG00000227232.4	WASH7P	6.50895977020264	10.7456922531128	6.67049884796143	6.3844690322876
#ENSG00000243485.2	MIR1302-11	0	0	0	0
#ENSG00000237613.2	FAM138A	0	0	0	0
#ENSG00000268020.2	OR4G4P	0	0	0	0
#ENSG00000240361.1	OR4G11P	0	0	0	0
#ENSG00000186092.4	OR4F5	0	0	0	0


# two imput files

# read annotation file.
# make a dict
# get col 0, and col 5

# then read expr file:
#  skip first two lines
#  then get line 2, list of sample IDs

# for each type of tissue:
    #  new file needs first row SampID, GeneIDs...
    #  each row is a sample, then values.



annotfile = '/Users/davidgibbs/Data/gtex/GTEx_Data_V6_Annotations_SampleAttributesDS.txt'
#exprfile = '/Users/davidgibbs/Data/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct'
exprfile = '/Users/davidgibbs/Data/gtex/gtex_1000.tsv'
genefile = '/Users/davidgibbs/Data/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.genes.txt'

annots = open(annotfile,'r').read().strip().split('\n')
sampDict = dict()
tissueDict = dict()
tissueSet = set()
for i in range(2,len(annots)):
    bits = annots[i].strip().split('\t')
    sampDict[bits[0]] = bits[5]
    tissueSet.add(bits[5])
    if bits[5] not in tissueDict:
        tissueDict[bits[5]] = [bits[0]]
    else: # already have this tissue type
        tissueDict[bits[5]].append(bits[0])


print(len(sampDict))
print(len(tissueDict))
print(len(tissueSet))

# get the gene list
genes = open(genefile, 'r').read().strip().split('\n')
genes = genes[3:]

print(len(genes))
print(genes[0:10])

        sampList = tissueDict[tissue] # get list of samples
        exprs = open(exprfile, 'r')   # open the expression file
        linect = 0                    # count lines as we move through it
        exprMat = ''
        sampVec = ''
        idxVec = ''
        print('    ' + str(len(sampList)))
        skip1 = exprs.readline()
        skip2 = exprs.readline()
        line  = exprs.readline()   # has sample lines
        sampLine = line.strip().split('\t')  # split up the line
        sampVec = [si for si in sampLine if si in sampList]
        idxVec = [i for i, si in enumerate(sampLine) if si in sampList]

        for line in exprs:
            bits = line.strip().split('\t') # split up the line
            exprVec = [float(bits[i]) for i in idxVec] # expression for one gene
            if len(sampVec) != len(exprVec):
                print("error in vector lengths!!")
                exit()
            if len(sampVec) == 0:
                continue
            if exprMat == '':
                exprMat = [[] for x in sampVec]  # if we haven't constructed the matrix yet, then build it
                sampFound = len(sampVec)
            for i, si in enumerate(exprVec):     # for each sample
                exprMat[i].append(exprVec[i])    # append *this* gene's expression to each sample's vector.

        fout = open('gtex' + tissue + '_rpkm.tsv', 'w')
        fout.write('\t'.join(['SampID'] + genes + ['\n']))
        for i,x in enumerate(exprMat):  # each row indexes a sample ID
            fout.write('\t'.join([sampVec[i]] + [str(xi) for xi in x] + ['\n']))







