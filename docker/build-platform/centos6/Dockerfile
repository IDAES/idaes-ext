FROM centos:centos6
MAINTAINER The CentOS Project <cloud-ops@centos.org>

RUN yum -y update; yum clean all
RUN yum -y install make gcc gcc-c++ gcc-gfortran
RUN yum -y install wget unzip patch which file lapack-devel blas-devel
RUN yum -y install http://opensource.wandisco.com/centos/6/git/x86_64/wandisco-git-release-6-1.noarch.rpm
RUN yum -y install git
RUN yum -y install epel-release
RUN yum -y install cmake3
RUN ln -s /usr/bin/cmake3 /usr/bin/cmake

##
# Pull down the repo
##
RUN mkdir repo

RUN cd /repo; \
    wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.bz2;\
    tar --bzip2 -xf boost_1_61_0.tar.bz2
