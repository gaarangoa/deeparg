import json

# train a classifier
from sklearn.feature_extraction import DictVectorizer
from sklearn import preprocessing
import numpy as np


from sklearn.metrics import classification_report, accuracy_score
from lasagne import layers
from lasagne import init

import lasagne

from lasagne.updates import sgd,nesterov_momentum
from nolearn.lasagne import NeuralNet

import numpy as np

from sklearn.utils import shuffle
import theano

floatX = theano.config.floatX


min_prob = {'aminoglycoside':  0.17, 
            'macrolide-lincosamide-streptogramin': 0.15, 
            'polymyxin': 0.18, 
            'nitrofuratoin': 0.5, 
            'mupirocin': 0.5, 
            'bacitracin': 0.3, 
            'quinolone': 0.16, 
            'chloramphenicol': 0.46, 
            'streptothricin': 0.5, 
            'noARGs': 0.01, 
            'puromycin': 0.5, 
            'thiostrepton': 0.5, 
            'unknown': 0.5, 
            'sulfonamide': 0.5, 
            'aminocoumarin': 0.5, 
            'beta_lactam': 0.2, 
            'kasugamycin': 0.5, 
            'peptide': 0.17, 
            'tetracenomycin': 0.5, 
            'multidrug': 0.16, 
            'tetracycline': 0.13, 
            'rifampin': 0.5, 
            'fosmidomycin': 0.5, 
            'triclosan': 0.2, 
            'fosfomycin': 0.5, 
            'na_antimicrobials': 0.5, 
            'qa_compound': 0.5, 
            'pyrimethamine': 0.5, 
            'elfamycin': 0.5, 
            'glycopeptide': 0.15, 
            'tunicamycin': 0.01, 
            'fusidic_acid': 0.3}



def make_xy(alignments, Features):
    SF = {i:True for i in Features}
    #   aignments is a dict of dicts. Each key is the read_id and the value is 
    #   a dict where the keys are the features(genes) and the values are the BitScores
    samples = alignments.keys()
    X = alignments.values()
    #   make sure that all the features are included in the alignments. To do so, we just 
    #   need to take one of the samples and add all the featuers. 
    for i in SF:
        try:
            a = X[0][i]
        except:
            X[0].update({i:0})
    
    Y=[i.split("|")[2] for i in samples]
    
    return [samples, X ,Y]

def best_hit(item):
    for j in align[i]:
        if float(align[i][j])>BH:
            BH = align[i][j]
            sel=j
    if BH>=50:
        return sel.split('|')[3]
    else:
        return 'unknown'

def main(deepL, alignments):
    S, X, Y = make_xy(alignments, deepL['features'])
    h = DictVectorizer(sparse=False)

    X = h.fit_transform(X)
    X = X.astype(floatX)

    # normalize feature matrix
    min_max_scaler = preprocessing.MinMaxScaler()
    X=min_max_scaler.fit_transform(X)

    clf = deepL['clf']

    proba = clf.predict_proba(X)
    #preds = clf.predict(X)
    preds = []
    deepL['Y_rev'].update({10000:"unclassified"})
    labels = []
    for ipx,p in enumerate(proba):
        px = np.argmax(p)
        if p[px] > min_prob[deepL['Y_rev'][px]]:
            preds.append(px)
            labels.append(Y[ipx])
        else:
            labels.append(Y[ipx])
            preds.append(10000)

    report = classification_report(labels, [deepL['Y_rev'][i] for i in preds]);
    print(report);

    return {'preds':[deepL['Y_rev'][i] for i in preds], 
            'labels': labels,
            'Y': Y,
            'proba':proba, 
            'report':report,
            'samples':S}
