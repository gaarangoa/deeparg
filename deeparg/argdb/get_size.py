import sys
from Bio import SeqIO

fi = sys.argv[1]

for record in SeqIO.parse(fi, "fasta"):
    gid = record.id;
    leng = len(str(record.seq))
    print(gid+"\t"+str(leng)) 