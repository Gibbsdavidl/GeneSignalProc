
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



def rfModel(dirs, exprfile, genes):
    inputs = open(dirs + exprfile, 'r').read().strip().split("\n")

    # for each gene set
    #    subset the data
    #    run cross validation
    #    take the mean for CV score

    # return the results as a table, tree, gene-set, score

    return(1)