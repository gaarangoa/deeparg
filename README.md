
# DeepARG

DeepARG predicts Antibiotic Resistance Genes (ARGs) from metagenomic sequences.
This repository now includes a modern PyTorch/Transformers runtime for the
original deepARG-SS and deepARG-LS models.

The default runtime uses the Hugging Face model bundle:

    gaarangoa/deeparg

That bundle contains:

    * deepARG-SS weights and config
    * deepARG-LS weights and config
    * DeepARG v2 DIAMOND databases
    * DIAMOND 2.1.24 binaries for Linux x86_64 and Linux aarch64
    * Greengenes gg13 Bowtie2 index for optional 16S normalization
    * source files for the modern pipeline

## Quick Start

Install `uv` if needed:

    curl -LsSf https://astral.sh/uv/install.sh | sh

Create and sync the environment from this repository:

    uv sync
    source .venv/bin/activate

Run DeepARG with minimal parameters:

    deeparg predict \
        --model LS \
        -i ./test/ORFs.fa \
        -o ./test/X

For short reads use `SS`:

    deeparg predict \
        --model SS \
        -i ./test/ORFs.fa \
        -o ./test/X

Run the modern short reads pipeline on an already-clean FASTA:

    deeparg short_reads_pipeline \
        --input-fasta ./test/ORFs.fa \
        --output-file ./test/reads

The first run downloads the model, database, and matching DIAMOND binary from
Hugging Face into the local Hugging Face cache.

## Python API

    import deeparg

    outputs = deeparg.predict(
        input_file="./test/ORFs.fa",
        output_file="./test/X",
        model="LS",
    )

    print(outputs["arg"])
    print(outputs["potential_arg"])

The returned paths are:

    * `<output>.mapping.ARG`
    * `<output>.mapping.potential.ARG`
    * `<output>.align.daa.tsv`

Short reads pipeline API:

    result = deeparg.short_reads_pipeline(
        input_fasta="./test/ORFs.fa",
        output_file="./test/reads",
    )

    print(result.arg_quant_file)

## What Gets Downloaded

By default, the modern pipeline loads everything from Hugging Face:

    * model repo: `gaarangoa/deeparg`
    * model subfolder: `LS` or `SS`, selected from `--model`
    * database: `<LS|SS>/database/features.dmnd`
    * gene lengths: `<LS|SS>/database/features.gene.length` when needed
    * 16S normalization index: `gg13/dataset.*`
    * DIAMOND binary:
        * `bin/linux-x86_64/diamond`
        * `bin/linux-aarch64/diamond`

DIAMOND resolution order:

    1. `--diamond-path`
    2. `diamond` on `PATH`
    3. `<data-path>/bin/diamond`
    4. bundled DIAMOND from `gaarangoa/deeparg`

If you already have a native DIAMOND installed, DeepARG will use it. Otherwise,
on Linux x86_64 or Linux aarch64, DeepARG downloads and uses the bundled binary.

## Common Commands

Use nucleotide input:

    deeparg predict \
        --model LS \
        --type nucl \
        -i input.fasta \
        -o output/sample

Use protein input:

    deeparg predict \
        --model LS \
        --type prot \
        -i proteins.faa \
        -o output/sample

Use a specific DIAMOND binary:

    deeparg predict \
        --model LS \
        -i input.fasta \
        -o output/sample \
        --diamond-path /path/to/diamond

Use a local DeepARG data directory instead of the HF database:

    deeparg predict \
        --model LS \
        -i input.fasta \
        -o output/sample \
        --data-path /path/to/deeparg-data

Use another Hugging Face repo or local `save_pretrained` folder:

    deeparg predict \
        --model LS \
        -i input.fasta \
        -o output/sample \
        --hf-model-path USER_OR_ORG/deeparg

Force CPU:

    deeparg predict \
        --model LS \
        -i input.fasta \
        -o output/sample \
        --hf-device cpu

## Short Reads Pipeline

The modern short reads pipeline keeps the legacy workflow but runs DeepARG
through the PyTorch/Transformers API:

    1. optional paired-end trimming with Trimmomatic
    2. optional paired-end merge with vsearch
    3. DeepARG-SS prediction
    4. ARG interval merge and quantification in Python
    5. optional 16S normalization with Bowtie2 and samtools

Minimal run from an already-clean FASTA:

    deeparg short_reads_pipeline \
        --input-fasta input.clean.fasta \
        --output-file output/sample

Run from paired FASTQ files:

    deeparg short_reads_pipeline \
        --forward-pe-file sample_R1.fastq.gz \
        --reverse-pe-file sample_R2.fastq.gz \
        --output-file output/sample

Enable 16S normalization. The `gg13` Bowtie2 index is downloaded from
`gaarangoa/deeparg` automatically:

    deeparg short_reads_pipeline \
        --input-fasta input.clean.fasta \
        --output-file output/sample \
        --normalize-16s

Raw paired FASTQ mode requires `trimmomatic` and `vsearch` on `PATH`.
`--normalize-16s` additionally requires `bowtie2` and `samtools` on `PATH`.
The model, ARG database, DIAMOND binaries, and `gg13` index are all in the HF
bundle.

If Python HTTP clients hang while `curl` can reach Hugging Face, force IPv4:

    HF_FORCE_IPV4=1 deeparg predict \
        --model LS \
        -i input.fasta \
        -o output/sample

## Output

DeepARG writes two main files:

    * `*.ARG`
    * `*.potential.ARG`

`*.ARG` contains predictions with probability greater than or equal to
`--min-prob`, which defaults to `0.8`. `*.potential.ARG` contains lower
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

## Development With uv

Create or refresh the environment:

    uv sync

Update dependency resolution:

    uv lock --upgrade
    uv sync

Run the CLI from the environment:

    .venv/bin/deeparg predict --help

Run a quick local test:

    HF_FORCE_IPV4=1 .venv/bin/deeparg predict \
        --model LS \
        -i ./test/ORFs.fa \
        -o ./data/test_hf_real \
        --hf-device cpu

The project is configured by:

    * `pyproject.toml`
    * `uv.lock`

## Training and Conversion

Train a new Transformers model from a DIAMOND/BLAST TSV:

    python -m deeparg.modern.train_transformers \
        --train-alignments /path/to/train_reads.tsv \
        --output-dir /path/to/deeparg-ss-transformers \
        --pipeline reads \
        --identity 30 \
        --evalue 1 \
        --coverage 25

Convert existing DeepARG Lasagne weights:

    python -m deeparg.modern.export_legacy \
        --metadata-pkl /path/to/model/v2/metadata_SS.pkl \
        --model-pkl /path/to/model/v2/model_SS.pkl \
        --database-dir /path/to/database/v2 \
        --feature-lengths /path/to/database/v2/features.gene.length \
        --output-dir /path/to/deeparg-ss-transformers

Upload a bundle to Hugging Face:

    HF_FORCE_IPV4=1 python -m deeparg.modern.upload_hf \
        --repo-id USER_OR_ORG/deeparg \
        --folder /path/to/deeparg-transformers

## Legacy Runtime

The original DeepARG runtime was based on Python 2.7, Theano, Lasagne, and
nolearn. It is deprecated in this repository; the active codebase now contains
only the PyTorch/Transformers runtime described above.

Legacy installation and usage notes were moved to:

    README_LEGACY.md

## Citation

If you use DeepARG in published research, please cite:

Arango-Argoty GA, Garner E, Pruden A, Heath LS, Vikesland P, Zhang L.
DeepARG: A deep learning approach for predicting antibiotic resistance genes
from metagenomic data. Microbiome 2018 6:23.

https://doi.org/10.1186/s40168-018-0401-z

## License

DeepARG is under the MIT license. Please also review the commercial restrictions
of the databases used during the mining process, including CARD, ARDB, and
UniProt.

## Contact

gustavo1@vt.edu
