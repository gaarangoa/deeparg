import json
from tqdm import tqdm
# this script takes a output from BLAST or diamond tab delimited format and
# creates a dictionary (key values) used latter for the representation of the
# matches in the database.


def make_alignments_json(fname, iden=50, eval=1e-5, coverage=0.8, BitScore=True, Features={}, glen={}, pipeline='reads'):
    alignments = {}
    BHit = {}
    SF = {i: True for i in Features}

    # print SF
    print([iden, eval, coverage, BitScore, pipeline])
    print("traversing input file {} ...".format(fname))
    if BitScore == True:
        measure = 11
    else:
        measure = 2

    for i in tqdm(open(fname), unit="reads"):
        i = i.strip().split("\t")
        l1 = [float(i[7])-float(i[6])]
        l2 = [float(i[9])-float(i[8])]

        # coverage = l1/l2

        # if i[0].split("|")[1]=="UNIPROT" and float(i[2])==100 and coverage==1: continue # check if there is a duplicate, only enabled in training, remove it for testing
        if float(i[2]) < iden:
            continue  # if the alignment has an identity below the threshold
        if float(i[10]) > eval:
            continue  # if the alignment has an evalue greater than the threshold

        # This is done for full genes
        if pipeline == 'genes':
            if float(int(i[3]))/glen[i[1]] < coverage:
                continue  # if the length of the alignment is less than the minimum coverage of the read

        # for reads the minimum length is 25 aminoacids i[3] is length in aminoacids
        # i[6] and i[7] are nucleotides
        if pipeline == 'reads':
            if float(int(i[3])) <= coverage:
                continue

        try:
            if SF[i[1]]:
                alignments[i[0]].update({
                                        i[1]: float(i[measure])
                                        })
        except:
            try:
                if SF[i[1]]:
                    alignments[i[0]] = {
                        i[1]: float(i[measure])
                    }
            except Exception as inst:
                print(inst)
                pass

        # compute the best hit for each entry
        try:
            if SF[i[1]]:
                try:
                    if BHit[i[0]][1] < float(i[measure]):
                        BHit[i[0]] = [i[1], float(i[measure]), i]
                except Exception as e:
                    BHit[i[0]] = [i[1], float(i[measure]), i]
        except:
            pass

    #json.dump(alignments, open(fname+".BitScoreMatrix.json",'w'))
    print(len(alignments), " reads passed the filters and ready for prediction")
    return [alignments, BHit]
