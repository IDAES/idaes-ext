FROM mcr.microsoft.com/windows:20H2

ENV MSYSTEM MINGW64
ENV CMAKE_GENERATOR "MSYS Makefiles"

WORKDIR c:/windows/temp

SHELL ["powershell", "-Command", "$ErrorActionPreference = 'Stop'; $ProgressPreference = 'SilentlyContinue';"]

RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
  Invoke-WebRequest -UseBasicParsing -uri "https://github.com/msys2/msys2-installer/releases/download/nightly-x86_64/msys2-base-x86_64-latest.sfx.exe" -OutFile msys2.exe; \
  .\msys2.exe -y -oC:\; \
  Remove-Item msys2.exe ; \
  function msys() { C:\msys64\usr\bin\bash.exe @('-lc') + @Args; } \
  msys ' '; \
  msys 'pacman --noconfirm -Syuu'; \
  msys 'pacman --noconfirm -Syuu'; \
  msys 'pacman --noconfirm -Scc';

RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
  function msys() { C:\msys64\usr\bin\bash.exe @('-lc') + @Args; } \
  msys 'pacman -S --needed --noconfirm --noprogressbar mingw-w64-x86_64-toolchain diffutils mingw-w64-x86_64-cmake mingw-w64-x86_64-make'; \
  msys 'pacman -S --needed --noconfirm --noprogressbar mingw-w64-x86_64-boost unzip patch make git python mingw-w64-x86_64-lapack';

RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
  function msys() { C:\msys64\usr\bin\bash.exe @('-lc') + @Args; } \
  msys 'mkdir /c/repo'; \
  msys 'mkdir /c/repo/petsc-dist';

WORKDIR c:/repo

# The PETSc build really doesn't like MinGW on docker, so need to put extra effort in
RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
  function msys() { C:\msys64\usr\bin\bash.exe @('-lc') + @Args; } \
  msys 'cd /c/repo && wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.0.tar.gz'; \
  msys 'cd /c/repo && tar -zxvf petsc-3.17.0.tar.gz';
SHELL ["C:\\msys64\\usr\\bin\\bash.exe", "-l", "-c"]
RUN 'cd /c/repo/petsc-3.17.0 && export CPATH=/c/repo/petsc-dist/include && /usr/bin/python configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 --download-metis-cmake-arguments=-G"MSYS Makefiles" --download-metis --download-mumps --with-mumps-serial=1 --with-make-exec=mingw32-make --prefix=/c/repo/petsc-dist'

# Why have I done this weird thing?  It will allow the command to fail, but keep going
SHELL ["powershell", "-Command", "$ErrorActionPreference = 'Stop'; $ProgressPreference = 'SilentlyContinue';"]
RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
  function msys() { C:\msys64\usr\bin\bash.exe @('-lc') + @Args; } \
  msys 'cd /c/repo/petsc-3.17.0 && make PETSC_DIR=/c/repo/petsc-3.17.0 PETSC_ARCH=arch-mswin-c-opt all'
SHELL ["C:\\msys64\\usr\\bin\\bash.exe", "-l", "-c"]
RUN 'cd /c/repo/petsc-3.17.0 && mingw32-make PETSC_DIR=/c/repo/petsc-3.17.0 PETSC_ARCH=arch-mswin-c-opt all'
RUN 'cd /c/repo/petsc-3.17.0 && /usr/bin/python config/install.py'

RUN cmd /c 'taskkill /F /FI "MODULES eq msys-2.0.dll"'
