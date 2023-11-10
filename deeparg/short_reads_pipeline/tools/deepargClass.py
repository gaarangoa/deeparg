import os
import sys


def run(R, data):
    # print sys.path
    try:
        cmd = " ".join(
            ['deeparg predict',
             '--type nucl',
             '--model SS',
             '-d', data['deep_arg_parameters']['data_path'],
             '-i', R,
             '-o', R+'.deeparg',
             '--arg-alignment-identity', str(data['deep_arg_parameters']['identity']),
             '--min-prob', str(data['deep_arg_parameters']['probability']),
             '--arg-alignment-evalue', str(data['deep_arg_parameters']['evalue']),
             ])
        print(cmd)
        x = os.popen(cmd).read()
        return True
    except Exception as inst:
        print(str(inst))
        return False


def dsize(path_to_data):
    return {i.split()[0].split("|")[-1].upper(): i.split() for i in open('{}/database/v2/features.gene.length'.format(path_to_data))}
