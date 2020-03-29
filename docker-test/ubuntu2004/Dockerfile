FROM ubuntu:20.04

ARG repo=https://github.com/idaes/idaes-pse.git
ARG branch=master

RUN apt-get -y update; \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get -y upgrade; \
    apt-get -y install libgfortran4 libgomp1 liblapack3 libblas3 git wget

##
# Pull down the repo
##
RUN mkdir repo

COPY ./extras/ /repo/

RUN cd repo; \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; \
    bash ~/miniconda.sh -b -p $HOME/miniconda; \
    eval "$(/root/miniconda/bin/conda shell.bash hook)"; \
    git clone ${repo}; \
    cd idaes-pse; \
    git checkout ${branch}; \
    echo "cloned repo"; \
    conda create -n idaes python=3.7 pip psutil; \
    conda activate idaes; \
    python setup.py develop; \
    cd ..; \
    # Once we update get-extensions to get specific linux builds need to modify a bit
    mv idaes-lib-ubuntu2004-64.tar.gz idaes-lib-linux-64.tar.gz; \
    mv idaes-solvers-ubuntu2004-64.tar.gz idaes-solvers-linux-64.tar.gz; \
    idaes get-extensions --url file:////repo/

RUN cd repo; \
    eval "$(/root/miniconda/bin/conda shell.bash hook)"; \
    conda activate idaes; \
    cd idaes-pse/idaes/generic_models; \
    pytest; \
    cd ../core; \
    pytest; \
    cd ../power_generation; \
    pytest
