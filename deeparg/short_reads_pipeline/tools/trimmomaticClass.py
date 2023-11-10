import os
import sys


def pairedEnd(R1, R2):
    print(sys.path)
    try:
        os.popen(" ".join(
            ['trimmomatic PE -phred33 ',
             R1,
             R2,
             R1+'.paired',
             R1+'.unpaired',
             R2+'.paired',
             R2+'.unpaired',
             'LEADING:3',
             'TRAILING:3',
             'SLIDINGWINDOW:4:15',
             'MINLEN:36'
             ])).read()
        return True
    except:
        return False
