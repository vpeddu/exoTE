FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

# install dependencies from pip3

RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python-biopython \
                   python3-pip \
                   python3-pysam \
                   wget \
                   curl \
                   libxml2-dev \
                   libssl-dev \
                   libcurl4-openssl-dev \
		           libz-dev \
		           bedtools \
                   samtools \
                   r-base \
                   software-properties-common \
                   unzip && \
    pip3 install biopython \
                 ete3 \
		 cython \
		 pysam \
		 alfa 

# updating R 
RUN apt install -y --no-install-recommends software-properties-common dirmngr && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt -y install r-base



RUN apt-get -y update && \
    apt-get install -y default-jdk
    #r-cran-rjava
    #r-cran-nloptr 
    #libssh2-1-dev \
    #libcurl4-gnutls-dev \
    # libgit2-dev \
    # libfontconfig1-dev \
    # libcairo2-dev \
    # libxml2-dev \
    # libcurl4-openssl-dev \
    # libssl-dev

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
RUN R -e "BiocManager::install('rtracklayer')"
RUN R -e "install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('vdiffr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggstance',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplotify',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('chromoMap',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggpmisc',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('ensembldb')"



