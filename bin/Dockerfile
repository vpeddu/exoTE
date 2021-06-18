FROM continuumio/miniconda3

WORKDIR /app

# TODO: Once the repo gets pushed update the dockerfile to pull from the repo
# wget exists in continuumio/miniconda3 so it should work 
RUN wget https://raw.githubusercontent.com/vpeddu/exoTE/main/bin/env.yml

#COPY env.yml .
RUN conda env create -f env.yml


SHELL ["conda", "run", "-n", "TE_Flow", "/bin/bash", "-c"]
#ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "TE_Flow"]

EXPOSE 5003
# The code to run when container is started:
#ENTRYPOINT ["conda", "activate", "TE_Flow"]
