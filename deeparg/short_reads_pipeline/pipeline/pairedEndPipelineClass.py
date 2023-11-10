import deeparg.short_reads_pipeline.tools.trimmomaticClass as trim
import deeparg.short_reads_pipeline.tools.vsearchClass as vsearch
import deeparg.short_reads_pipeline.tools.deepargClass as deeparg
import deeparg.short_reads_pipeline.quantification.quantificationClass as quant
import deeparg.short_reads_pipeline.pipeline.d16spipelineClass as D16sPipe
import deeparg.short_reads_pipeline.quantification.normalizationClass as norm

import os

d16sPipe = D16sPipe.d16sPipe()


class PairedEnd():
    def __init__(self, data):
        self.info = ''
        self.pairedR1File = data['pairedR1File']
        self.pairedR2File = data['pairedR2File']
        self.sample_name = data['sample_output_file']
        self.data = data

    def run(self):

        print('Step 1: Trimming and QC using Trimmomatic')
        if not trim.pairedEnd(self.pairedR1File, self.pairedR2File):
            return 0

        print('\n\n\nStep 2: Merging paired end reads using Vsearch')
        if not vsearch.merge(self.pairedR1File+'.paired', self.pairedR2File+'.paired', self.sample_name):
            return 0

        print('\n\n\nStep 3: Run DeepARG-SS to identify ARG-like reads')
        if not deeparg.run(self.sample_name+'.clean', self.data):
            return 0

        print('\n\n\nStep 4: Quantification of ARG-like counts')
        if not quant.merge(self.sample_name+'.clean.deeparg.mapping.ARG', self.data['deep_arg_parameters']['data_path']):
            return 0

        print('\n\n\nStep 5: Normalize to 16S rRNAs - this may take a while ...')
        if not d16sPipe.run(self.sample_name+'.clean', ggdata="{}/data/gg13/dataset".format(self.data['deep_arg_parameters']['data_path'], )):
            return 0

        norm.normalize(
            self.sample_name + '.clean.sorted.bam.merged.quant',
            self.sample_name + '.clean.deeparg.mapping.ARG.merged.quant',
            float(self.data['parameters']['coverage'])/100,
            self.data['parameters']
        )
