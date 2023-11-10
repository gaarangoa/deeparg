import json
from sklearn.feature_extraction import DictVectorizer
from sklearn import preprocessing
import numpy as np

from lasagne import layers
from lasagne import init

import lasagne
from lasagne.updates import sgd,nesterov_momentum
from nolearn.lasagne import NeuralNet

import theano

import sys
import process_blast
import cPickle
import options as opt
from model import model

floatX = theano.config.floatX

min_prob = {'aminoglycoside': 0.17, 
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
    
    return [samples, X]

def main(deepL, alignments):
    S, X = make_xy(alignments, deepL['features'])
    h = DictVectorizer(sparse=False)

    X = h.fit_transform(X)
    X = X.astype(floatX)

    # normalize feature matrix
    min_max_scaler = preprocessing.MinMaxScaler()
    X=min_max_scaler.fit_transform(X)
    
    clf = NeuralNet(
        layers=model(deepL['input_nodes'], deepL['output_nodes']),
        update=nesterov_momentum,
        update_learning_rate=0.01,
        update_momentum=0.9,
        regression=False,
        max_epochs=100,
        verbose=2,
    );

    clf.load_params_from(opt.path+"/model/model.pkl");

    proba = clf.predict_proba(X);

    preds = []
    probs = []
    deepL['Y_rev'].update({10000:"unclassified"})

    for ipx,p in enumerate(proba):
        px = np.argmax(p);
        if p[px] > min_prob[deepL['Y_rev'][px]]:
            preds.append(px);
            probs.append(p[px]);
        else:
            preds.append(10000);
            probs.append(0);

    return [[S[ix], deepL['Y_rev'][i], probs[ix]] for ix,i in enumerate(preds)]


fi = sys.argv[1];
fo = open(sys.argv[2], 'w');
# fi = opt.path+"/db/argsdb.reads.align.test.tsv";

print "Loading deep learning model";
deepL = cPickle.load(open(opt.path+"/model/metadata.pkl"));

print "Loading sample to analyze";
align = process_blast.make_alignments_json(fi, iden=0, eval=1e-5, len=5, BitScore=True);

print "Predicting AR-like reads";
predict = main(deepL, align);

print "Computing relative abundance";

for i in predict:
    fo.write("\t".join([str(j) for j in i])+"\n");

fo.close()