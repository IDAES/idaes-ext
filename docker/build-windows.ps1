$flavor = $args[0]
$buildarg_1 = $args[1]  # just use this to pass in --no-cache or some such

$repo = "https://github.com/eslickj/idaes-ext.git"
$branch = "master"

IF ($flavor -eq "windows"){
  $wdir = "c:/repo"
}
ELSEIF ($flavor -eq "centos6"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "centos7"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "centos8"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "ubuntu1804"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "ubuntu1910"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "ubuntu2004"){
  $wdir = "/repo"
}
ELSE{
  echo "Please specify a flavor in {windows, centos6, centos7, centos8,"
  echo "                            ubuntu1804, ubuntu1910, ubuntu2004}."
  exit 1
}

xcopy /y /E extras ${flavor}\extras\

cd ${flavor}
docker build --rm ${buildarg_1} --build-arg repo=${repo} --build-arg branch=${branch} -t ${flavor}_build .
Remove-Item extras -Recurse -Force -Confirm:$false
docker run --name ${flavor}_build_tmp -dt ${flavor}_build:latest
docker stop ${flavor}_build_tmp
docker cp ${flavor}_build_tmp:${wdir}/idaes-ext/dist-lib/idaes-lib-${flavor}-64.tar.gz .
try{
  docker cp ${flavor}_build_tmp:${wdir}/idaes-ext/dist-solvers/idaes-solvers-${flavor}-64.tar.gz .
}
catch{
  echo "Solvers were not built."
}
docker rm ${flavor}_build_tmp
cd ..
