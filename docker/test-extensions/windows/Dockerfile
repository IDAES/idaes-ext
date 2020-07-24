FROM mcr.microsoft.com/windows:1909

ARG repo=https://github.com/idaes/idaes-pse.git
ARG branch=master

SHELL ["powershell", "-command"]

WORKDIR c:/windows/temp

RUN Invoke-WebRequest -UserAgent 'DockerCI' -outfile Miniconda3-latest-Windows-x86_64.exe https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
RUN Start-Process Miniconda3-latest-Windows-x86_64 -ArgumentList '/InstallationType=JustMe /RegisterPython=0 /S /D=%UserProfile%\Miniconda3' -Wait
RUN mkdir C:\repo
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
    cd idaes-pse\idaes; \
    pytest tests\test_have_pynumero.py; \
    pytest tests\test_ipopt.py;\
    cd generic_models; \
    pytest; \
    cd ../core; \
    pytest; \
    cd ../power_generation; \
    pytest
