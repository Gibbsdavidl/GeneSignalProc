'''
This file contains script for emulating the ssGSEA algorithm from Barbie et al. 2009

Original code by:  AndrewZhaoLuo

'''

import numpy as np
#from simulation import *

def calculate_enrichment_score(gene_set, expressions, omega):
    """
    Given a gene set, a map of gene names to expression levels, and a weight omega, returns the ssGSEA
    enrichment score for the gene set as described by *D. Barbie et al 2009*
    :requires: every member of gene_set is a key in expressions
    :param gene_set: a set of gene_names in the set
    :type gene_set: set
    :param expressions: a dictionary mapping gene names to their expression values
    :type expressions: dict
    :param omega: the weighted exponent on the :math:`P^W_G` term.
    :type omega: float
    :returns: an array representing the intermediate Enrichment Scores for each step along the sorted gene list.
              To find the total enrichment score, take the sum of all values in the array.
    """

    #first sort by absolute expression value, starting with the highest expressed genes first
    keys_sorted = sorted(expressions, key=expressions.get, reverse=True)

    #values representing the ECDF of genes in the geneset
    P_GW_numerator = 0
    P_GW_denominator = 0

    #determining denominator value
    i = 1 #current rank stepping through listing of sorted genes
    for gene in keys_sorted:
        if gene in gene_set:
            P_GW_denominator += i ** omega
        i += 1

    P_GW = lambda : P_GW_numerator / P_GW_denominator

    #values representing the ECDF of genes not in the geneset
    P_NG_numerator = 0
    P_NG_denominator = len(expressions) - len(gene_set)
    P_NG = lambda : P_NG_numerator / P_NG_denominator

    #integrate different in P_GW and P_NG
    i = 1 #current index in the traversal of sorted genes
    scores = []
    for gene in keys_sorted:
        if gene in gene_set:
            P_GW_numerator += i ** omega
        else:
            P_NG_numerator += 1

        scores.append(P_GW() - P_NG())
        i += 1

    return sum(scores)



def buildExprDict(exprdat, geneNames):
    edict = dict()
    bits = exprdat.split('\t')
    for ai, bi in enumerate(bits[1:]):
        edict[ai] = float(bi)
    return(edict)


def printSSGSEAResults(res0, dirs):
    fout = open(dirs+'ssgsea_scores.tsv','w')
    for lx in res0:
        drddddd = [str(x) for x in lx]
        fout.write('\t'.join(drddddd)+'\n')
    return(0)


def scoreSets(dirs, geneSets, exprfile, omega):
    inputs = open(dirs + exprfile, 'r').read().strip().split("\n")
    geneNames = inputs[0].split('\t')
    allResults = []
    for i in inputs[1:]:
        thisResult = []
        for g in geneSets:
            di = buildExprDict(i, geneNames)
            scr = calculate_enrichment_score(g, di, omega)
            thisResult.append(scr)
        allResults.append(thisResult)
    printSSGSEAResults(allResults, dirs)
    return(allResults)
