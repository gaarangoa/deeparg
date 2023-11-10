
import json
import utils

path = '/Users/gustavoarango/Documents/projects/ARGdb/db/'

metadata = json.load(open(path+'metadata_with_multilabels.json'))
clusters = json.load(open(path+'CLUSTERS.csv.json'))
alignments = json.load(open(path+'alignment.tsv.json'))

iden = 90

DB={}

for cluster in clusters:
    cluster_stats = utils.consensus_type(clusters[cluster], metadata, alignments)
    DB.update(utils.cluster_annotation(cluster_stats, 90, 1e-10, metadata))


vari = 'AnFactor'
from collections import Counter
Counter([DB[i][vari] for i in DB])


json.dump(DB, open(path+"ANNOTATION_FULL2.json",'w'))



