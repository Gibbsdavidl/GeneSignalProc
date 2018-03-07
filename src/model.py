
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn import datasets

import numpy as np

def rfModelTest():
    iris = datasets.load_iris()
    clf = RandomForestClassifier(max_depth=2, random_state=0)
    scores = cross_val_score(clf, iris.data, iris.target, cv=10)
    return(scores)



def rfModel(dirs, exprfile, pheno, genes):
    ys = [int(yi) for yi in open(dirs + pheno,'r').read().strip().split('\n')]
    inputs = open(dirs + exprfile, 'r').read().strip().split("\n")
    inputHeader = inputs.pop(0)
    scoreList = []
    # for each gene set
    for i, gs in enumerate(genes):
        print(i)
        # subset the data into xs
        xs = []
        for j, rowj in enumerate(inputs):
            rowj =  rowj.strip().split('\t')
            rowj.pop(0)
            numj = [float(jj) for ji,jj in enumerate(rowj) if ji in gs]
            xs.append(numj)
        #    run cross validation on this gene set
        clf = RandomForestClassifier(max_depth=2)
        scores = cross_val_score(clf, xs, ys, cv=10)
        scoreList.append(np.mean(scores))
    #    take the mean for CV score
    print(scoreList)
    # return the results as a table, tree, gene-set, score

    return(scoreList)