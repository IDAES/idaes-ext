FROM mcr.microsoft.com/windows:20H2

SHELL ["powershell", "-command"]

WORKDIR c:/windows/temp
ENV MSYSTEM MINGW64
ENV CMAKE_GENERATOR "MSYS Makefiles"

SHELL ["powershell", "-Command", "$ErrorActionPreference = 'Stop'; $ProgressPreference = 'SilentlyContinue';"]

RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
  Invoke-WebRequest -UseBasicParsing -uri "https://github.com/msys2/msys2-installer/releases/download/2022-01-28/msys2-base-x86_64-20220128.sfx.exe" -OutFile msys2.exe; \
  .\msys2.exe -y -oC:\; \
  Remove-Item msys2.exe ; \
  function msys() { C:\msys64\usr\bin\bash.exe @('-lc') + @Args; } \
  msys ' '; \
  msys 'pacman --noconfirm -Syuu'; \
  msys 'pacman --noconfirm -Syuu'; \
  msys 'pacman --noconfirm -Scc'; \
  msys 'pacman -S --needed --noconfirm --noprogressbar mingw-w64-x86_64-toolchain diffutils mingw-w64-x86_64-cmake'; \
  msys 'pacman -S --needed --noconfirm --noprogressbar mingw-w64-x86_64-boost unzip patch make git python mingw-w64-x86_64-lapack mingw-w64-x86_64-msmpi';
