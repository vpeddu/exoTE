FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

# install dependencies from pip3

RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python-biopython \
                   python3-pip \
                   python3-pysam \
                   wget \
		           libz-dev \
		           bedtools \
                   samtools \
                   r-base \
                   unzip && \
    pip3 install biopython \
                 ete3 \
		 cython \
		 pysam \
		 alfa 
RUN apt-get -y update && \
    apt-get install -y default-jdk \
    r-cran-rjava \
    r-cran-nloptr \
    libssh2-1-dev \
    libcurl4-gnutls-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev

RUN wget https://github.com/mhammell-laboratory/TEtranscripts/archive/refs/heads/master.zip && \
    unzip master.zip && \
    cd TEtranscripts-master &&\
    python3 setup.py install


RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('ape',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('viridis')"
RUN R -e "BiocManager::install('reshape2')"
RUN R -e "BiocManager::install('tidyverse')"
RUN R -e "BiocManager::install('Rsubread')"
RUN R -e "BiocManager::install('Rsamtools')"
RUN R -e "BiocManager::install('biomaRt')"
RUN R -e "install.packages('vdiffr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggstance',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplotify',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('chromoMap',dependencies=TRUE, repos='http://cran.rstudio.com/')"