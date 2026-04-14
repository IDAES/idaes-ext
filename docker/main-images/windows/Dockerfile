FROM mcr.microsoft.com/windows/servercore:ltsc2022

ENV MSYSTEM=MINGW64 \
    CHERE_INVOKING=1 \
    CMAKE_GENERATOR="MSYS Makefiles" \
    MSYS2_PATH_TYPE=inherit

WORKDIR C:/Windows/Temp

SHELL ["powershell", "-NoLogo", "-ExecutionPolicy", "Bypass", "-Command", "$ErrorActionPreference='Stop'; $ProgressPreference='SilentlyContinue';"]

# MSYS2 installation
RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    Invoke-WebRequest -UseBasicParsing -Uri "https://github.com/msys2/msys2-installer/releases/download/2025-08-30/msys2-base-x86_64-20250830.sfx.exe" -OutFile msys2.exe; \
    Start-Process -NoNewWindow -Wait -FilePath msys2.exe -ArgumentList '-y', '-oC:\'; \
    Remove-Item msys2.exe -Force

WORKDIR C:/repo

# Bring in HSL tarball
# COPY ./coinhsl.tar.gz C:/repo/coinhsl.tar.gz

# Download PETsc
RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    Invoke-WebRequest -UseBasicParsing -Uri "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.20.6.tar.gz" -OutFile petsc-3.20.6.tar.gz; \
    & 'C:\msys64\usr\bin\bash.exe' -lc 'cd /c/repo && tar -zxvf petsc-3.20.6.tar.gz && mv petsc-3.20.6 petsc'

# Install packages
RUN ["C:\\msys64\\usr\\bin\\bash.exe", "-lc", "\
set -euo pipefail; \
pacman -Sy --noconfirm; \
pacman -S --needed --noconfirm --noprogressbar \
  mingw-w64-x86_64-toolchain \
  mingw-w64-x86_64-cmake \
  mingw-w64-x86_64-make \
  mingw-w64-x86_64-lapack \
  mingw-w64-x86_64-boost \
  mingw-w64-x86_64-gcc-fortran \
  mingw-w64-x86_64-meson \
  mingw-w64-x86_64-ninja; \
pacman -S --needed --noconfirm --noprogressbar \
  base-devel git wget unzip patch which file diffutils findutils python; \
pacman --noconfirm -Scc; \
mkdir -p /c/repo/petsc-dist /c/repo/coinhsl-dist; \
rm -rf /c/\\$Recycle.Bin/* \
"]

# PETSc build
RUN ["C:\\msys64\\usr\\bin\\bash.exe", "-lc", "\
    cd /c/repo/petsc && \
    export CPATH=/c/repo/petsc-dist/include && \
    /usr/bin/python configure \
      --with-debug=0 \
      --with-shared=0 \
      --with-mpi=0 \
      --with-fortran-bindings=0 \
      --download-metis \
      --download-mumps \
      --with-64-bit-indices \
      --with-mumps-serial=1 \
      --download-metis-cmake-arguments='-DCMAKE_POLICY_VERSION_MINIMUM=3.5' \
      --with-make-exec=mingw32-make \
      --prefix=/c/repo/petsc-dist \
"]

RUN ["C:\\msys64\\usr\\bin\\bash.exe", "-lc", "\
    cd /c/repo/petsc && \
    make PETSC_ARCH=arch-mswin-c-opt all \
"]

RUN ["C:\\msys64\\usr\\bin\\bash.exe", "-lc", "\
    cd /c/repo/petsc && /usr/bin/python config/install.py \
"]

# Set entrypoint
CMD ["C:\\msys64\\usr\\bin\\bash.exe", "-l"]
