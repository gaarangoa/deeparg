import json

# this script takes a output from BLAST or diamond tab delimited format and 
# creates a dictionary (key values) used latter for the representation of the
# matches in the database. 

def make_alignments_json(fname, iden=0, eval=1e-5, len=25, BitScore=True):
    alignments = {}
    print "traversing file ..."

    if BitScore==True:
        measure = 11
    else:
        measure = 2

    for i in open(fname):
        i = i.strip().split("\t")
        l1 = max([float(i[7]), float(i[6])])
        l2 = max([float(i[9]), float(i[8])])

        coverage = min([l1,l2])/max([l1,l2])

        # if i[0].split("|")[1]=="UNIPROT" and float(i[2])==100 and coverage==1: continue # check if there is duplicate, only enabled in training, remove it for testing
        if float(i[2])<iden or float(i[10])>eval or int(i[3])-int(i[4])<len: continue
        try:
            alignments[i[0]].update({
                                    i[1]:float(i[measure])
                                    })
        except:
            alignments[i[0]] = {
                                    i[1]:float(i[measure])
                                    }
    print 'storing file ...'
    #json.dump(alignments, open(fname+".BitScoreMatrix.json",'w'))
    return alignments

