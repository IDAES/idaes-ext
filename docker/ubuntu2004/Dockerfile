FROM ubuntu:20.04

ENV repo https://github.com/eslickj/idaes-ext.git
ENV branch dock

RUN apt-get -y update; \
    apt-get -y upgrade; \
    apt-get -y install ubuntu-dev-tools liblapack-dev libblas-dev unzip

RUN apt-get -y install gfortran libboost-all-dev

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
    bash scripts/compile_libs.sh ubuntu2004; \
    bash scripts/compile_solvers.sh ubuntu2004
