
import json
import numpy as np
from sklearn.cross_validation import train_test_split
from lasagne import layers
import lasagne
from lasagne.updates import sgd, nesterov_momentum
from nolearn.lasagne import NeuralNet
import numpy as np
import theano
from collections import Counter
from model import model

floatX = theano.config.floatX


def make_xy_per_class(samples, X, Y, features):
    labels = Counter([i.split("|")[2] for i in samples])
    feature_dict = {i.split("|")[0]: True for i in features}

    x_train = []
    y_train = []
    x_test = []
    y_test = []
    x_test_names = []

    for label in list(set(Y)):
        index = Y == label
        x_sel = X[index]
        y_sel = Y[index]
        s_sel = samples[index]

        train_size = 1
        # if len(y_sel) < 5: train_size = 0.9999;

        ixtrain, ixtest, ytrain, ytest = train_test_split(
            range(x_sel.shape[0]), y_sel, random_state=0, train_size=train_size)

        x_train += x_sel[ixtrain].tolist()
        y_train += ytrain.tolist()

        x_test_p = []
        y_test_p = []
        x_train_p = []
        y_train_p = []

        for i in ixtest:
            try:
                a = feature_dict[s_sel[i].split("|")[0]]
                x_train_p.append(x_sel[i])
                y_train_p.append(y_sel[i])
            except:
                x_test_p.append(x_sel[i])
                y_test_p.append(y_sel[i])
                x_test_names.append(s_sel[i])

        x_test += x_test_p
        y_test += y_test_p
        x_train += x_train_p
        y_train += y_train_p

    return [np.array(x_train, dtype="float32"), np.array(y_train, dtype="int32"), np.array(x_test, dtype="float32"), np.array(y_test, dtype="int32"), x_test_names]


def main(data):

    X = data['X_red']
    Y = data['Y']
    sample_names = np.array(data['samples'])
    features = data['features']
    Y_rev = data['Y_rev']

    # x_train, y_train, x_test, y_test, x_test_names = make_xy_per_class(sample_names, X, Y, features)

    x_train = X
    y_train = Y

    print "dataset: ", X.shape, Y.shape
    print "training set: ", x_train.shape, y_train.shape
    # print "testing set: ", x_test.shape, y_test.shape

    clf = NeuralNet(
        layers=model(x_train.shape[1], len(set(Y))),
        update=nesterov_momentum,
        update_learning_rate=0.01,
        update_momentum=0.9,
        regression=False,
        max_epochs=100,
        verbose=2,
    )

    clf.fit(x_train, y_train)

    return {
        'clf': clf,
        'parameters': {
            "features": features,
            "Y_rev": Y_rev,
            "input_nodes": x_train.shape[1],
            "output_nodes": len(set(Y)),
        }
    }
