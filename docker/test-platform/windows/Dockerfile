FROM mcr.microsoft.com/windows:1909

SHELL ["powershell", "-command"]

WORKDIR c:/windows/temp

RUN Invoke-WebRequest -UserAgent 'DockerCI' -outfile Miniconda3-latest-Windows-x86_64.exe https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
RUN Start-Process Miniconda3-latest-Windows-x86_64 -ArgumentList '/InstallationType=JustMe /RegisterPython=0 /S /D=%UserProfile%\Miniconda3' -Wait
RUN mkdir C:\repo
RUN Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
RUN cmd /c %UserProfile%\Miniconda3\condabin\conda.bat init
RUN cmd /c conda create -n idaes python=3.7.4 pip psutil git pytest
RUN Invoke-WebRequest -UserAgent 'DockerCI' -outfile Miniconda3-latest-Windows-x86_64.exe https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
RUN Start-Process Miniconda3-latest-Windows-x86_64 -ArgumentList '/InstallationType=JustMe /RegisterPython=0 /S /D=%UserProfile%\Miniconda3' -Wait
