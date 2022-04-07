FROM eslickj/idaes-ext-ubuntu-2204-build:latest

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
    bash scripts/compile_solvers.sh ubuntu2204; \
    bash scripts/compile_libs.sh ubuntu2204
