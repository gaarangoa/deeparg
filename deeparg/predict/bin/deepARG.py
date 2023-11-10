# import modules outside the path
import os
import sys
# cwd = "/".join(os.getcwd().split("/")[:-2])
# sys.path.append(cwd)

import json
from sklearn.feature_extraction import DictVectorizer
from sklearn import preprocessing
import numpy as np

from lasagne import layers
from lasagne import init

import lasagne
from lasagne.updates import sgd, nesterov_momentum
from nolearn.lasagne import NeuralNet

import theano

from deeparg.predict.bin import process_blast
import cPickle
# import options as opt
from deeparg.predict.bin.model import model

import math

from itertools import islice

from tqdm import tqdm
import logging 

logger = logging.getLogger()

floatX = theano.config.floatX

def chunks(data, size=10000):
    it = iter(data)
    for i in xrange(0, len(data), size):
        yield {k: data[k] for k in islice(it, size)}


def make_xy(alignments, Features):
    SF = {i: True for i in Features}
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
            X[0].update({i: 0})

    return [samples, X]


def main(deepL, clf, alignments, version):
    S, X = make_xy(alignments, deepL['features'])
    h = DictVectorizer(sparse=False)

    # print("Liberating ram memory")
    del alignments

    X = h.fit_transform(X)
    X = X.astype(floatX)

    # normalize feature matrix
    min_max_scaler = preprocessing.MinMaxScaler()
    X = min_max_scaler.fit_transform(X)

    proba = clf.predict_proba(X)

    preds = []
    probs = []
    deepL['Y_rev'].update({10000: "unclassified"})

    SP = []

    for ix, p in enumerate(proba):
        npx = np.argsort(p)
        px = npx[-1]

        if 0.0 < p[px] < 0.9:
            preds.append(npx[-2])
            probs.append(p[npx[-2]])
            SP.append(S[ix])
            preds.append(px)
            probs.append(p[px])
            SP.append(S[ix])
        elif p[px] >= 0.9:
            preds.append(px)
            probs.append(p[px])
            SP.append(S[ix])

    return [[SP[ix], deepL['Y_rev'][i], probs[ix]] for ix, i in enumerate(preds)]


import operator


def process(fin, fon, iden, version, evalue, prob, minCoverage, pipeline, version_m, args):

    fi = fin  # first parameter is the input file

    # fi = opt.path+"/db/argsdb.reads.align.test.tsv";
    logger.info("Loading deep learning model ...")
    deepL = cPickle.load(
        open(args.data_path+"/model/"+version_m+"/metadata"+version+".pkl"))
    clf = NeuralNet(
        layers=model(deepL['input_nodes'], deepL['output_nodes']),
        update=nesterov_momentum,
        update_learning_rate=0.01,
        update_momentum=0.9,
        regression=False,
        max_epochs=100,
        verbose=2,
    )

    clf.load_params_from(args.data_path+"/model/"+version_m+"/model"+version+".pkl")

    # print deepL['features']

    logger.info("loading gene lengths")
    glen = {i.split()[0]: float(i.split()[1]) for i in open(
        args.data_path+"/database/"+version_m+"/features.gene.length")}

    logger.info("Loading sample to analyze")
    [align, BH] = process_blast.make_alignments_json(
        fi, iden=iden, eval=evalue, coverage=minCoverage, BitScore=True, Features=deepL['features'], glen=glen, pipeline=pipeline)

    # print align
    logger.info("Predicting ARG-like reads: Running deepARG" +
          version+" model version "+version_m)
    logger.info("input dataset is splitted into chunks of 10000 reads")

    chunks_input = chunks(align, size=10000)
    predict = []
    for _chunk in tqdm(chunks_input, total=int(len(align)/10000), unit="chunks"):
        predict += main(deepL, clf, _chunk, version)

    # print predict
    logger.info("Predicting ARGs")

    fo = open(fon+'.ARG', 'w')  # second parameter is the output file
    fo2 = open(fon+'.potential.ARG', 'w')
    fo.write("#ARG\tquery-start\tquery-end\tread_id\tpredicted_ARG-class\tbest-hit\tprobability\tidentity\talignment-length\talignment-bitscore\talignment-evalue\tcounts\n")
    fo2.write("#ARG\tquery-start\tquery-end\tread_id\tpredicted_ARG-class\tbest-hit\tprobability\tidentity\talignment-length\talignment-bitscore\talignment-evalue\tcounts\n")
    for i in tqdm(predict):

        # 1 Get the alignments for that sample
        x_align = align[i[0]]
        # 2 Get only the alignments with the predicted label
        x_align = {o: x_align[o] for o in x_align.keys() if "|"+i[1]+"|" in o}
        # 3 Compute the best Hit

        if x_align:
            x_bh = max(x_align.iteritems(), key=operator.itemgetter(1))[0]
            bs_bh = x_align[x_bh]
            # print(x_bh, align[i[0]][x_bh])
            # print(BH[i[0]])
            if i[2] >= prob:
                fo.write("\t".join([
                    # gene where read is from (subtype)
                    x_bh.split("|")[-1].upper(),
                    BH[i[0]][2][8],  # alignment gene start
                    BH[i[0]][2][9],  # alignment gene end
                    i[0],  # read-id
                    i[1],  # predicted type
                    x_bh,  # best hit
                    str(i[2]),  # probability
                    BH[i[0]][2][2],  # identity
                    BH[i[0]][2][3],  # alignment length
                    BH[i[0]][2][-1],  # bitscore
                    BH[i[0]][2][-2],  # evalue
                    '1'  # count
                ])+"\n"
                )
            else:
                x_bh = BH[i[0]][0]
                bs_bh = BH[i[0]][1]
                fo2.write("\t".join([
                    # gene where read is from (subtype)
                    x_bh.split("|")[-1].upper(),
                    BH[i[0]][2][8],  # alignment gene start
                    BH[i[0]][2][9],  # alignment gene end
                    i[0],  # read-id
                    i[1],  # predicted type
                    x_bh,  # best hit
                    str(i[2]),  # probability
                    BH[i[0]][2][2],  # identity
                    BH[i[0]][2][3],  # alignment length
                    BH[i[0]][2][-1],  # bitscore
                    BH[i[0]][2][-2],  # evalue
                    '1'  # count
                ])+"\n"
                )
        else:
            x_bh = BH[i[0]][0]
            bs_bh = BH[i[0]][1]
            fo2.write("\t".join([
                # gene where read is from (subtype)
                x_bh.split("|")[-1].upper(),
                BH[i[0]][2][8],  # alignment gene start
                BH[i[0]][2][9],  # alignment gene end
                i[0],  # read-id
                i[1],  # predicted type
                "undefined",  # best hit
                str(i[2]),  # probability
                BH[i[0]][2][2],  # identity
                BH[i[0]][2][3],  # alignment length
                BH[i[0]][2][-1],  # bitscore
                BH[i[0]][2][-2],  # evalue
                # x_bh, # problematic-classification: when the prediction class has a different best feature (the arg cannot be defined - probably because errors in the database )
                '1'
            ])+"\n"
            )

        # remove entries with low coverage
        # if float(BH[i[0]][2][3])/glen[x_bh] < minlenper: continue

    fo.close()
    fo2.close()
