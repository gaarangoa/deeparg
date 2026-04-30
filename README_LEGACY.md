# DeepARG Legacy Runtime

This document preserves the original DeepARG Python 2.7 / Theano / Lasagne
usage notes. The legacy implementation is deprecated and no longer part of the
active package code in this branch. New usage should use the modern
PyTorch/Transformers runtime in `README.md`.

## Legacy Requirements

The historical DeepARG runtime was designed for:

    * Python 2.7.15 - 2.7.18
    * Theano 0.8.2
    * Lasagne 0.1
    * nolearn 0.6
    * scikit-learn 0.19.1
    * DIAMOND 0.9.24

The modern runtime no longer requires these packages.

## Legacy Installation With conda

Create a Python 2.7 environment:

    conda create -n deeparg_legacy python=2.7.18
    conda activate deeparg_legacy

Install DIAMOND:

    conda install -c bioconda diamond==0.9.24

Optional tools used by the historical short reads pipeline:

    conda install -c bioconda trimmomatic
    conda install -c bioconda vsearch
    conda install -c bioconda bedtools==2.29.2
    conda install -c bioconda bowtie2==2.3.5.1
    conda install -c bioconda samtools

Install DeepARG:

    pip install git+https://github.com/gaarangoa/deeparg.git

Download the original data:

    deeparg download_data -o /path/to/deeparg-data/

The original data is hosted at:

    https://zenodo.org/records/8280582

## Legacy Prediction

Classify ORFs from assembled contigs:

    deeparg predict \
        --model LS \
        -i ./test/ORFs.fa \
        -o ./test/X \
        -d /path/to/deeparg-data/ \
        --type nucl \
        --min-prob 0.8 \
        --arg-alignment-identity 30 \
        --arg-alignment-evalue 1e-10 \
        --arg-num-alignments-per-entry 1000

## Legacy Model Choices

Use `LS` for long sequences or gene-like sequences:

    deeparg predict \
        --model LS \
        --type nucl \
        --input-file /path/file.fasta \
        --output-file /path/out

Use `SS` for short sequence reads:

    deeparg predict \
        --model SS \
        --type nucl \
        --input-file /path/file.fasta \
        --output-file /path/out

Protein input is also supported:

    deeparg predict \
        --model LS \
        --type prot \
        --input-file /path/file.faa \
        --output-file /path/out

## Legacy Output

DeepARG writes two files:

    * `*.ARG`
    * `*.potential.ARG`

`*.ARG` contains sequences with prediction probability greater than or equal to
the configured probability threshold. `*.potential.ARG` contains lower
probability ARG-like predictions.

Columns:

    * ARG_NAME
    * QUERY_START
    * QUERY_END
    * QUERY_ID
    * PREDICTED_ARG_CLASS
    * BEST_HIT_FROM_DATABASE
    * PREDICTION_PROBABILITY
    * ALIGNMENT_BESTHIT_IDENTITY (%)
    * ALIGNMENT_BESTHIT_LENGTH
    * ALIGNMENT_BESTHIT_BITSCORE
    * ALIGNMENT_BESTHIT_EVALUE
    * COUNTS

## Legacy Short Reads Pipeline

Run the historical short reads pipeline:

    deeparg short_reads_pipeline \
        --forward_pe_file ./test/F.fq.gz \
        --reverse_pe_file ./test/R.fq.gz \
        --output_file ./test/reads \
        -d /path/to/deeparg-data \
        --bowtie_16s_identity 100

This path may require Trimmomatic, vsearch, bedtools, bowtie2, and samtools.

## Notes

The legacy notes are kept for reproducibility and reference. The old runtime
depends on Python 2 and scientific Python packages that are difficult to install
on modern systems. For new runs, use the PyTorch/Transformers runtime documented
in `README.md`.
