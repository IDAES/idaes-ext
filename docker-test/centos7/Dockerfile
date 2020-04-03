FROM centos:centos7
MAINTAINER The CentOS Project <cloud-ops@centos.org>

ARG repo=https://github.com/idaes/idaes-pse.git
ARG branch=master

RUN yum -y update; yum clean all
RUN yum -y install lapack blas libgfortran libgomp git wget

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
    conda create -n idaes python=3.7.4 pip psutil; \
    conda activate idaes; \
    python setup.py develop; \
    cd ..; \
    idaes get-extensions --url file:////repo/ --platform centos7

RUN cd repo; \
    eval "$(/root/miniconda/bin/conda shell.bash hook)"; \
    conda activate idaes; \
    cd idaes-pse/idaes/generic_models; \
    pytest; \
    cd ../core; \
    pytest; \
    cd ../power_generation; \
    pytest
