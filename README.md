
# DeepARG
A deep learning based approach to predict Antibiotic Resistance Genes (ARGs) from metagenomes. It provides two models,deepARG-SS and deepARG-LS.

<!-- https://zenodo.org/records/8280582/files/deeparg.zip?download=1 -->

## Latest Release 
* updated on Nov 10 - 2023
* deeparg 1.0.2: Added to pip
* Fastq input - short reads pipeline fixed

<!-- ## Web Service
We have released a web service to process raw sequences (paired end) using deepARG. You will get absolute and relative abundances of each submitted sample. You can find the website at http://bench.cs.vt.edu/deeparg -->

#### DeepARG output
DeepARG generates two files: *.ARG that contains the sequences with a probability >= --prob (0.8 default) and *.potential.ARG with sequences containing a probability < --prob (0.8 default). The *.potential.ARG file can still contain ARG-like sequences, howevere, it is necessary inspect its sequences.

The output format for both files consists of the following fields:

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

# Installation
DeepARG is under Python 2.7, therefore, it is recommended to run it via virtual environment or via docker.  

## Instal via miniconda
Install miniconda https://docs.conda.io/en/latest/miniconda.html

<!-- ## Dependencies
DeepARG requires the next software installed:

    1 diamond version 0.9.24 (http://www.diamondsearch.org/index.php?pages/installation/)
    2 Python version 2.7.14 - 2.7.18 -->

<!-- ## Installation -->
<!-- ### Install using virtualenv
(optional) if not pip is installed in yor machine:

    1. wget https://bootstrap.pypa.io/get-pip.py
    2. python get-pip.py

Create a virutal environment with virtualenv

    pip3 install virtualenv
    python3 -m virtualenv env
    source env/bin/activate

Install deeparg with pip and download the data required by deeparg

    pip install deeparg==1.0.2
    Download data from zenodo manually following this link: https://zenodo.org/records/8280582/files/deeparg.zip?download=1 and save it to your local directory
    3. Run deeparg initiate_environment -i /path/to/local/directory/deeparg/

To re-activate the virtual environment:

    source env/bin/activate -->
    
### Use conda environment
Create a virtual environment with conda:

    conda create -n deeparg_env python=2.7.18
    source activate deeparg_env

Install diamond with conda (inside virtual environment): 

    conda install -c bioconda diamond==0.9.24

Optional (used for short reads pipeline): 

    conda install -c bioconda trimmomatic
    conda install -c bioconda vsearch
    conda install -c bioconda bedtools==2.29.2
    conda install -c bioconda bowtie2==2.3.5.1
    conda install -c bioconda samtools

Install deeparg with pip and download the data required by deeparg

    pip install git+https://github.com/gaarangoa/deeparg.git
    deeparg download_data -o /path/to/local/directory/

Activate virtual environment

    conda activate deeparg_env

Deactivate the virtual environment:

    conda deactivate

<!-- ### Known instalation issues
See issues on repository

### Docker

    docker pull gaarangoa/deeparg:latest

To run deeparg using docker just type:

    docker run --rm -v $PWD:/data/ gaarangoa/deeparg:latest deeparg --help

Note that input parameters are under /data/ directory.  -->

### Example:
In this example, we will classify a set of ORFs from a set of assembled contigs. The fasta file contains gene sequences (nucleotides).

    deeparg predict \
        --model LS \
        -i ./test/ORFs.fa \
        -o ./test/X \
        -d /path/to/data/ \
        --type nucl \
        --min-prob 0.8 \
        --arg-alignment-identity 30 \
        --arg-alignment-evalue 1e-10 \
        --arg-num-alignments-per-entry 1000

### Usage

    usage: deeparg predict 
        -h, --help            show this help message and exit
        --model MODEL         Select model to use (short sequences for reads | long
                                sequences for genes) SS|LS [No default]
        -i INPUT_FILE, --input-file INPUT_FILE
                                Input file (Fasta input file)
        -o OUTPUT_FILE, --output-file OUTPUT_FILE
                                Output file where to store results
        -d DATA_PATH, --data-path DATA_PATH
                                Path where data was downloaded [see deeparg download-
                                data --help for details]
        --type TYPE           Molecular data type prot/nucl [Default: nucl]
        --min-prob MIN_PROB   Minimum probability cutoff [Default: 0.8]
        --arg-alignment-identity ARG_ALIGNMENT_IDENTITY
                                Identity cutoff for sequence alignment [Default: 50]
        --arg-alignment-evalue ARG_ALIGNMENT_EVALUE
                                Evalue cutoff [Default: 1e-10]
        --arg-alignment-overlap ARG_ALIGNMENT_OVERLAP
                                Alignment read overlap [Default: 0.8]
        --arg-num-alignments-per-entry ARG_NUM_ALIGNMENTS_PER_ENTRY
                                Diamond, minimum number of alignments per entry
                                [Default: 1000]
        --model-version MODEL_VERSION
                                Model deepARG version [Default: v2]

### Usage examples:

Go to the deeparg-ss directory and run any of the following commands:

Input is a FASTA file:

    1) Annotate gene-like sequences when the input is a nucleotide FASTA file:
        deeparg predict --model LS --type nucl --input /path/file.fasta --out /path/to/out/file.out

    2) Annotate gene-like sequences when the input is an amino acid FASTA file:
        deeparg predict --model LS --type prot --input /path/file.fasta --out /path/to/out/file.out

    3) Annotate short sequence reads when the input is a nucleotide FASTA file:
        deeparg predict --model SS --type nucl --input /path/file.fasta --out /path/to/out/file.out

    3) Annotate short sequence reads when the input is a protein FASTA file (unusual case):
        deeparg predict --model SS --type prot --input /path/file.fasta --out /path/to/out/file.out

# Short reads pipeline

## Usage

    deeparg short_reads_pipeline [-h] --forward_pe_file FORWARD_PE_FILE
                                        --reverse_pe_file REVERSE_PE_FILE
                                        --output_file OUTPUT_FILE
                                        [-d DEEPARG_DATA_PATH]
                                        [--deeparg_identity DEEPARG_IDENTITY]
                                        [--deeparg_probability DEEPARG_PROBABILITY]
                                        [--deeparg_evalue DEEPARG_EVALUE]
                                        [--gene_coverage GENE_COVERAGE]
                                        [--bowtie_16s_identity BOWTIE_16S_IDENTITY]

    optional arguments:
    -h, --help            show this help message and exit
    --forward_pe_file FORWARD_PE_FILE
                            forward mate from paired end library
    --reverse_pe_file REVERSE_PE_FILE
                            reverse mate from paired end library
    --output_file OUTPUT_FILE
                            save results to this file prefix
    -d DEEPARG_DATA_PATH, --deeparg_data_path DEEPARG_DATA_PATH
                            Path where data was downloaded [see deeparg download-
                            data --help for details]
    --deeparg_identity DEEPARG_IDENTITY
                            minimum identity for ARG alignments [default 80]
    --deeparg_probability DEEPARG_PROBABILITY
                            minimum probability for considering a reads as ARG-
                            like [default 0.8]
    --deeparg_evalue DEEPARG_EVALUE
                            minimum e-value for ARG alignments [default 1e-10]
    --gene_coverage GENE_COVERAGE
                            minimum coverage required for considering a full gene
                            in percentage. This parameter looks at the full gene
                            and all hits that align to the gene. If the overlap of
                            all hits is below the threshold the gene is discarded.
                            Use with caution [default 1]

### Example

    deeparg short_reads_pipeline \
        --forward_pe_file ./test/F.fq.gz \
        --reverse_pe_file ./test/R.fq.gz \
        --output_file ./test/reads \
        -d ~/Desktop/darg \
        --bowtie_16s_identity 100


## About
If you use deepARG in published research, please cite:

Arango-Argoty GA, Garner E, Pruden A, Heath LS, Vikesland P, Zhang L. DeepARG: A deep learning approach for predicting antibiotic resistance genes from metagenomic data. Microbiome20186:23
https://doi.org/10.1186/s40168-018-0401-z.

## Database
Database is hosted in Zenodo: https://zenodo.org/records/8280582

## License
deepARG is under the MIT licence. However, please take a look at te comercial restrictions of the databases used during the mining process (CARD, ARDB, and UniProt).

## Contact
If need any asistance please contact: gustavo1@vt.edu

<!-- ## Update pip
python setup.py sdist bdist_wheel
twine upload dist/* -->