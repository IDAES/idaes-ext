FROM centos:centos8
MAINTAINER The CentOS Project <cloud-ops@centos.org>

ARG repo=https://github.com/idaes/idaes-ext.git
ARG branch=master

RUN yum -y update; yum clean all
RUN yum -y install git make gcc gcc-c++ gcc-gfortran boost-devel cmake
RUN yum -y install wget unzip patch which file
RUN yum -y install dnf-plugins-core
RUN yum config-manager --set-enabled PowerTools
RUN dnf --enablerepo=PowerTools -y install lapack-devel
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
    bash scripts/compile_libs.sh centos8; \
    bash scripts/compile_solvers.sh centos8
