FROM centos:centos7
MAINTAINER The CentOS Project <cloud-ops@centos.org>

ARG repo=https://github.com/idaes/idaes-ext.git
ARG branch=master

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran boost-devel
RUN yum -y install wget unzip patch which file lapack-devel blas-devel
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
    bash scripts/compile_libs.sh centos7; \
    bash scripts/compile_solvers.sh centos7
