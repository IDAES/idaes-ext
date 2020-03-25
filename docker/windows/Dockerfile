FROM mcr.microsoft.com/windows:1909

ENV repo https://github.com/eslickj/idaes-ext.git
ENV branch dock

SHELL ["powershell", "-command"]

WORKDIR c:/windows/temp

RUN Invoke-WebRequest -UserAgent 'DockerCI' -outfile msys2-x86_64.tar.xz http://repo.msys2.org/distrib/msys2-x86_64-latest.tar.xz
RUN Invoke-WebRequest -UserAgent 'DockerCI' -outfile 7zsetup.exe http://www.7-zip.org/a/7z1900-x64.exe

RUN Start-Process .\7zsetup -ArgumentList '/S' -Wait
RUN Start-Process -FilePath "C:\Program` Files\7-Zip\7z.exe" -ArgumentList e, "C:\Windows\Temp\msys2-x86_64.tar.xz", `-oC:\Windows\Temp\ -NoNewWindow -PassThru -Wait
RUN Start-Process -FilePath "C:\Program` Files\7-Zip\7z.exe" -ArgumentList x, "C:\Windows\Temp\msys2-x86_64.tar", `-oC:\ -NoNewWindow -PassThru -Wait

SHELL ["C:\\msys64\\usr\\bin\\bash.exe", "-l", "-c"]

RUN 'exit'; \
    'pacman -Syuu --needed --noconfirm --noprogressbar'; \
    'pacman -Syu  --needed --noconfirm --noprogressbar'; \
    'pacman -Su   --needed --noconfirm --noprogressbar'; \
    'pacman -Scc --noconfirm';

RUN 'pacman -S --needed --noconfirm --noprogressbar mingw-w64-x86_64-toolchain diffutils'
RUN 'pacman -S --needed --noconfirm --noprogressbar mingw-w64-x86_64-boost unzip patch make git mingw-w64-x86_64-lapack'

ENV MSYSTEM MINGW64

RUN 'mkdir /c/repo/'
COPY ./extras/ c:/repo/
RUN 'export PATH=/c/msys64/mingw64/bin:$PATH   && cd /c/repo/ && git clone $repo && cd idaes-ext && git checkout $branch && sh scripts/compile_libs.sh windows';
RUN 'export PATH=/c/msys64/mingw64/bin:$PATH   && cd /c/repo/idaes-ext && sh scripts/compile_solvers.sh windows';


WORKDIR c:/repo/
