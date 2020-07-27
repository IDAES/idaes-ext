FROM eslickj/idaes-ext-windows-build:latest

ARG repo=https://github.com/idaes/idaes-ext.git
ARG branch=master

SHELL ["C:\\msys64\\usr\\bin\\bash.exe", "-l", "-c"]
ENV MSYSTEM MINGW64

COPY ./extras/ c:/repo/
RUN 'export PATH=/c/msys64/mingw64/bin:$PATH && cd /c/repo/ && git clone $repo && cd idaes-ext && git checkout $branch && sh scripts/compile_libs.sh windows';
RUN 'export PATH=/c/msys64/mingw64/bin:$PATH && cd /c/repo/idaes-ext && sh scripts/compile_solvers.sh windows';
