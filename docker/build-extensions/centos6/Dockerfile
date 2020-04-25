FROM centos6_idaes_build:latest

ARG repo=https://github.com/idaes/idaes-ext.git
ARG branch=master

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
    export BOOST_HEADER=/repo/boost_1_61_0; \
    bash scripts/compile_libs.sh centos6; \
    bash scripts/compile_solvers.sh centos6
