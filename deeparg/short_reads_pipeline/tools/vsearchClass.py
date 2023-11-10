import os
import sys


def merge(R1, R2, sample):
    print(sys.path)
    try:
        x = os.popen(" ".join(
            ['vsearch --fastq_mergepairs ',
             R1,
             '--reverse ', R2,
             '--fastaout ', R1+'.merged',
                             '--fastaout_notmerged_fwd', R1+'.unmerged',
                             '--fastaout_notmerged_rev', R2+'.unmerged'
             ])).read()

        x = os.popen(' '.join([
            'cat ',
            R1+'.merged',
            R1+'.unmerged',
            R2+'.unmerged',
            '>',
            sample+'.clean'
        ])).read()
        return True
    except:
        return False
