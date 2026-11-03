# Build Windows tarballs (written to tarballs\) in the Windows build image.
# Put coinhsl.tar.gz in extras\ to build with HSL. Requires docker in Windows-container mode.
# Example: .\build.ps1 -Branch main -SolverArgs "--without-hsl"
param(
  [string]$Repo = "https://github.com/IDAES/idaes-ext.git",
  [string]$Branch = "main",
  [string]$Image = "idaes-ext-windows:x86_64",
  [string]$SolverArgs = "",
  [string]$Memory = "8g"
)

$ErrorActionPreference = "Stop"

$flavor = "windows"
$mname = "x86_64"
$container = "idaes-ext-$flavor-$mname-build-tmp"
$bash = "C:\msys64\usr\bin\bash.exe"
$here = $PSScriptRoot

function Invoke-Docker {
  & docker @args
  if ($LASTEXITCODE -ne 0) { throw "docker $args failed (exit $LASTEXITCODE)" }
}

function Remove-Container {
  $ErrorActionPreference = "Continue"
  & docker rm -f $container 2>&1 | Out-Null
}

function Invoke-InContainer([string]$cmd) {
  Invoke-Docker exec $container $bash -lc $cmd
}

Write-Host "build.ps1 arguments:"
Write-Host "    Repo: $Repo"
Write-Host "    Branch: $Branch"
Write-Host "    Image: $Image"
Write-Host "    Memory: $Memory"
Write-Host "    Args to compile_solvers.sh: $SolverArgs"

Remove-Container
try {
  # Windows containers can't always be copied into while running, so copy before starting
  Invoke-Docker create --name $container --memory $Memory -t $Image
  $hsl = Join-Path $here "extras\coinhsl.tar.gz"
  if (Test-Path $hsl) {
    Invoke-Docker cp $hsl "${container}:c:/repo/coinhsl.tar.gz"
  }
  Invoke-Docker start $container

  Invoke-InContainer "cd /c/repo && git clone $Repo && cd idaes-ext && git checkout $Branch"
  Invoke-InContainer "cd /c/repo/idaes-ext && bash scripts/compile_solvers.sh $flavor $SolverArgs"
  Invoke-InContainer "cd /c/repo/idaes-ext && bash scripts/compile_libs.sh $flavor"
  Invoke-Docker stop $container

  $out = Join-Path $here "tarballs"
  New-Item -ItemType Directory -Force -Path $out | Out-Null
  Invoke-Docker cp "${container}:c:/repo/idaes-ext/dist-functions/idaes-functions-$flavor-$mname.tar.gz" $out
  Invoke-Docker cp "${container}:c:/repo/idaes-ext/dist/idaes-solvers-$flavor-$mname.tar.gz" $out
  Write-Host "Tarballs written to $out"
}
finally {
  Remove-Container
}
