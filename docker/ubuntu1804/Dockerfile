FROM ubuntu:18.04

ARG repo=https://github.com/idaes/idaes-ext.git
ARG branch=master

RUN apt-get -y update; \
    apt-get -y upgrade; \
    apt-get -y install ubuntu-dev-tools liblapack-dev libblas-dev unzip

RUN apt-get -y install gfortran libboost-all-dev cmake

##
# Pull down the repo
##
RUN mkdir repo

COPY ./extras/ /repo/

RUN cd repo; \
    git clone ${repo} ; \
    cd idaes-ext ; \
    git checkout ${branch}; \
    echo "cloned repo"

##
# Build all the things
##
RUN pwd; cd repo/idaes-ext; \
    bash scripts/compile_libs.sh ubuntu1804; \
    bash scripts/compile_solvers.sh ubuntu1804
