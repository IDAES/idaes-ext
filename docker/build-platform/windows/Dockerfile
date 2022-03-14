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
  msys 'pacman -S --needed --noconfirm --noprogressbar mingw-w64-x86_64-toolchain diffutils mingw-w64-x86_64-cmake'; \
  msys 'pacman -S --needed --noconfirm --noprogressbar mingw-w64-x86_64-boost unzip patch make git python mingw-w64-x86_64-lapack mingw-w64-x86_64-msmpi';

RUN [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
  function msys() { C:\msys64\usr\bin\bash.exe @('-lc') + @Args; } \
  msys 'wget https://raw.githubusercontent.com/eslickj/idaes-ext/main/scripts/compile_petsc_windows.sh'; \
  msys 'sh ./compile_petsc_windows.sh';


RUN taskkill /F /FI "MODULES eq msys-2.0.dll"
