FROM ubuntu:18.04

# install dependencies from pip3

RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python-biopython \
                   python3-pip \
                   python3-pysam \
                   wget \
		   libz-dev \
		   bedtools \
                   unzip && \
    pip3 install biopython \
                 ete3 \
		 cython \
		 pysam \
		 alfa