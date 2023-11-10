import numpy as np 
from sklearn.feature_extraction import DictVectorizer
from sklearn import preprocessing
# from sklearn.ensemble import ExtraTreesClassifier
# from sklearn.feature_selection import SelectFromModel
import theano
from collections import Counter
floatX = theano.config.floatX

# This scritp is designed to get the data from the processed blast file and generate a XY dataset
# of features. This script is made for testing purposes.



def make_xy2(alignments):
    #   alignments is a dict of dicts. Each key is the read_id and the value is 
    #   a dict where the keys are the features(genes) and the values are the BitScores
    print "getting sample names and matrix X"
    samples = alignments.keys()
    X = alignments.values()
    
    # XV = []
    # for sample in alignments:
    #     sfeat = alignments[sample]
    #     for feat in sfeat:
    #         XV.append([feat, sfeat[feat]])
    
    #   make sure that all the features are included in the alignments. To do so, we just 
    #   need to take one of the samples and add all the featuers. 

    #   first_hit = M[0]
    print "getting labels"
    Y=[i.split("|")[2] for i in samples]
    

    Y_dict = {i:ix for ix,i in enumerate(set(Y))}
    Y_rev = {y:x for x,y in Y_dict.iteritems()}
    
    Y = np.array([Y_dict[i] for i in Y], dtype='int32')

    print "Formating X to matrix"
    h = DictVectorizer(sparse=False)

      

    X = h.fit_transform(X)
    X = X.astype(floatX)

     

    feature_names=h.get_feature_names()
    print "Normalizing the matrix X"
    min_max_scaler = preprocessing.MinMaxScaler()
    NX=min_max_scaler.fit_transform(X)
    

    print "Feature selection"
    # Feature selection
    # Umbalanced dataset, give weigths to the small classes:

    # weights = Counter(Y)
    # weights = {i:1.2-float(weights[i])/(max(weights.values())) for i in weights}
    
    # clf = ExtraTreesClassifier( n_estimators=2000, 
    #                             bootstrap = True, 
    #                             n_jobs=5, 
    #                             # class_weight=weights, 
    #                             max_features="auto")
    # clf = clf.fit(NX, Y)
    # model = SelectFromModel(clf, prefit=True)
    
    # NX = model.transform(NX)
    
    # feature_dict={ feature_names[ix]:True for ix,i in enumerate(model.get_support()) if i==True}
    # feature_names = feature_dict.keys()

    return {"samples":samples, 
            "X": X, 
            "Y": Y, 
            "X_red":NX, 
            "features":feature_names,
            "Y_dict":Y_dict,
            "Y_rev":Y_rev}



