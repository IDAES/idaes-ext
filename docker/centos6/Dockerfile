FROM centos:centos6
MAINTAINER The CentOS Project <cloud-ops@centos.org>

ENV repo https://github.com/eslickj/idaes-ext.git
ENV branch dock

RUN yum -y update; yum clean all
RUN yum -y install make gcc gcc-c++ gcc-gfortran
RUN yum -y install wget unzip patch which file lapack-devel blas-devel
RUN yum -y install http://opensource.wandisco.com/centos/6/git/x86_64/wandisco-git-release-6-1.noarch.rpm
RUN yum -y install git


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

RUN cd repo/idaes-ext; \
    wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.bz2;\
    tar --bzip2 -xf /boost_1_61_0.tar.bz2

##
# Build all the things
##
RUN pwd; cd repo/idaes-ext; \
    export BOOST_HEADER=/repo/boost_1_61_0; \
    bash scripts/compile_libs.sh centos6; \
    bash scripts/compile_solvers.sh centos6
