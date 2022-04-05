$flavor = $args[0]
$buildarg_1 = $args[1]  # just use this to pass in --no-cache or some such

$repo = "https://github.com/idaes/idaes-ext.git"
$branch = "main"
$mname = "x86_64"

IF ($flavor -eq "windows"){
  $wdir = "c:/repo"
}
ELSEIF ($flavor -eq "centos7"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "centos8"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "rockylinux8"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "ubuntu1804"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "ubuntu2004"){
  $wdir = "/repo"
}
ELSE{
  echo "Please specify a flavor in {windows, el7, el8, ubuntu1804, ubuntu2004}."
  exit 1
}

xcopy /y /E extras ${flavor}\extras\

cd ${flavor}
docker build --rm ${buildarg_1} --build-arg repo=${repo} --build-arg branch=${branch} -t ${flavor}_build_itmp .
Remove-Item extras -Recurse -Force -Confirm:$false
docker run --name ${flavor}_build_tmp -dt ${flavor}_build_itmp:latest
docker stop ${flavor}_build_tmp
docker cp ${flavor}_build_tmp:${wdir}/idaes-ext/dist-lib/idaes-lib-${flavor}-${mname}.tar.gz .
try{
  docker cp ${flavor}_build_tmp:${wdir}/idaes-ext/dist-solvers/idaes-solvers-${flavor}-${mname}.tar.gz .
}
catch{
  echo "Solvers were not built."
}
try{
  docker cp ${flavor}_build_tmp:${wdir}/idaes-ext/dist-petsc/idaes-petsc-${flavor}-${mname}.tar.gz .
}
catch{
  echo "PETSc was not built."
}
docker rm ${flavor}_build_tmp
docker rmi ${flavor}_build_itmp
cp *.tar.gz ../tarballs/
Remove-Item *.tar.gz
cd ..
