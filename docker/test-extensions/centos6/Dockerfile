FROM centos6_idaes_test:latest

ARG repo=https://github.com/idaes/idaes-pse.git
ARG branch=master

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
    pip install -e .; \
    cd ..; \
    idaes get-extensions --url file:////repo/ --platform centos6

RUN cd repo; \
    eval "$(/root/miniconda/bin/conda shell.bash hook)"; \
    conda activate idaes; \
    cd idaes-pse/idaes; \
    pytest tests/test_have_pynumero.py; \
    pytest tests/test_ipopt.py;\
    cd ./generic_models; \
    pytest; \
    cd ../core; \
    pytest; \
    cd ../power_generation; \
    pytest
