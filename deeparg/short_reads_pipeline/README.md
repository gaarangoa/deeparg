## Short reads pipeline analysis

This pipeline has been designed to work under the docker contanier, it will process the reads as the same way as the online version of deepARG. 

        optional arguments:
        -h, --help            show this help message and exit
        --forward_pe_file FORWARD_PE_FILE
                                forward mate from paired end library
        --reverse_pe_file REVERSE_PE_FILE
        root@54addf3cb7b3:/data/short_reads_pipeline# python short_reads_pipeline.py --help
        usage: short_reads_pipeline.py [-h] --forward_pe_file FORWARD_PE_FILE
                                    --reverse_pe_file REVERSE_PE_FILE --output_file
                                    OUTPUT_FILE
                                    [--deeparg_identity DEEPARG_IDENTITY]
                                    [--deeparg_probability DEEPARG_PROBABILITY]
                                    [--deeparg_evalue DEEPARG_EVALUE]
                                    [--gene_coverage GENE_COVERAGE]
                                    [--bowtie_16s_identity BOWTIE_16S_IDENTITY]
                                    [--path_to_executables PATH_TO_EXECUTABLES]

        optional arguments:
        -h, --help            show this help message and exit
        --forward_pe_file FORWARD_PE_FILE
                                forward mate from paired end library
        --reverse_pe_file REVERSE_PE_FILE
                                reverse mate from paired end library
        --output_file OUTPUT_FILE
                                save results to this file prefix
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
        --bowtie_16s_identity BOWTIE_16S_IDENTITY
                                minimum identity a read as a 16s rRNA gene [default
                                0.8]
        --path_to_executables PATH_TO_EXECUTABLES
                                path where the deepARG program is installed