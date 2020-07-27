FROM eslickj/idaes-ext-windows-test:latest

ARG repo=https://github.com/idaes/idaes-pse.git
ARG branch=master

SHELL ["powershell", "-command"]

WORKDIR c:/windows/temp

COPY ./extras/ /repo/
RUN Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
RUN cmd /c %UserProfile%\Miniconda3\condabin\conda.bat init
RUN cmd /c conda create -n idaes python=3.7.4 pip psutil git pytest

RUN cd \repo; \
    conda activate idaes; \
    git clone $Env:repo;\
    cd idaes-pse; \
    git checkout $env:branch; \
    pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org -r requirements.txt; \
    pip install -e .; \
    idaes get-extensions --url file:///c:/repo/; \
    cd idaes; \
    pytest tests\test_have_pynumero.py; \
    pytest tests\test_ipopt.py;\
    cd generic_models; \
    pytest; \
    cd ../core; \
    pytest; \
    cd ../power_generation; \
    pytest
