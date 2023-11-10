import argparse
from deeparg.short_reads_pipeline.pipeline import pairedEndPipelineClass
import sys
import os
import time
import base64
import json


def main(args):

    deep_arg_parameters = dict(
        identity=args.deeparg_identity,
        probability=args.deeparg_probability,
        evalue=args.deeparg_evalue,
        data_path=args.deeparg_data_path
    )

    parameters = dict(
        coverage=args.gene_coverage,
        identity_16s_alignment=args.bowtie_16s_identity
    )

    data = dict(
        pairedR1File=args.forward_pe_file,
        pairedR2File=args.reverse_pe_file,
        deep_arg_parameters=deep_arg_parameters,
        sample_output_file=args.output_file,
        parameters=parameters,
    )

    pipe = pairedEndPipelineClass.PairedEnd(data)
    pipe.run()
