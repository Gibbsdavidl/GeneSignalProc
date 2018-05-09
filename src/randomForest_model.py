
# discover if a gene set can predict a given phenotypes using randomForest

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn import datasets
from sklearn.model_selection import LeaveOneOut

import numpy as np

def rfModelTest():
    iris = datasets.load_iris()
    clf = RandomForestClassifier(max_depth=2, random_state=0)
    scores = cross_val_score(clf, iris.data, iris.target, cv=10)
    return(scores)



def rfModel(dirs, exprfile, pheno, genes, cvs):
    print("random forest")
    ys = [int(yi) for yi in open(dirs + pheno,'r').read().strip().split('\n')]
    inputs = open(dirs + exprfile, 'r').read().strip().split("\n")
    inputHeader = inputs.pop(0)
    scoreList = []
    # for each gene set
    for i, gs in enumerate(genes):
        # subset the data into xs
        xs = []
        for j, rowj in enumerate(inputs):
            rowj = rowj.strip().split('\t')
            rowj.pop(0)
            numj = [float(jj) for ji,jj in enumerate(rowj) if ji in gs]
            xs.append(numj)
        #    run cross validation on this gene set
        clf = RandomForestClassifier(max_depth=5)
        scores = cross_val_score(clf, xs, ys, cv=cvs)
        scoreList.append(np.mean(scores))
    #    take the mean for CV score
    # return the results as a table, tree, gene-set, score
    return(scoreList)



def rfModelSetScores(dirs, inputs, pheno, genes, cvs):
    print("random forest")
    ys = [int(yi) for yi in open(dirs + pheno,'r').read().strip().split('\n')]
    #inputs = open(dirs + exprfile, 'r').read().strip().split("\n")
    #inputHeader = inputs.pop(0)
    # use all reported gene sets for prediction
    xs = np.array([ np.array(x) for x in inputs ])
    clf = RandomForestClassifier(max_depth=5, n_estimators=100)
    clf.fit(xs, ys)
    cvscores = cross_val_score(clf, xs, ys, cv=cvs, n_jobs=4)
    featImp = clf.feature_importances_
    scoreMean = mean(cvscores)
    #    take the mean for CV score
    # return the results as a table, tree, gene-set, score
    return(scoreMean, cvscores, clf, featImp)