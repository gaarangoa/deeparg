FROM conda/miniconda2

# install dependencies
RUN apt update
RUN apt-get --yes install default-jre
RUN apt-get --yes install bedtools
RUN apt-get --yes install bowtie2
RUN apt-get --yes install samtools
RUN apt-get --yes install wget 

# Install DeepARG 1.0.2
RUN pip install git+https://github.com/gaarangoa/deeparg.git

# Download data
RUN deeparg download_data -o ~/deeparg/

RUN yes | conda install -c bioconda diamond==0.9.24
RUN yes | conda install -c bioconda trimmomatic
RUN yes | conda install -c bioconda vsearch
RUN yes | conda install -c bioconda bedtools
RUN yes | conda install -c bioconda bowtie2
RUN yes | conda install -c bioconda samtools

