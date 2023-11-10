import sys
from Bio import SeqIO


fi=sys.argv[1] # fastafile from besthit

for record in SeqIO.parse(open(fi), "fasta"):
    idx = record.id
    print idx+"\t"+str(len(record.seq))