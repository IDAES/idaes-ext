FROM centos:centos7
MAINTAINER The CentOS Project <cloud-ops@centos.org>

RUN yum -y update; yum clean all
RUN yum -y install lapack blas libgfortran libgomp git wget

RUN mkdir repo
