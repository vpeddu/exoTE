# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-ubuntu:20.04

RUN apt-get -y update && apt-get install -y default-jdk r-cran-rjava r-cran-nloptr libssh2-1-dev libgit2-dev libfontconfig1-dev libcairo2-dev
RUN apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev

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
RUN R -e "BiocManager::install('rtracklayer')"
RUN apt-get install -y haskell-platform python3-pip 
RUN apt-get install -y pandoc
#RUN wget https://www.github.com/jgm/pandoc/releases/download/2.14.0.2/pandoc-2.14.0.2-linux-amd64.tar.gz && \
#    tar xvzf pandoc-2.14.0.2-linux-amd64.tar.gz --strip-components 1 -C /bin/