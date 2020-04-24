$flavor = $args[0]
$buildarg_1 = $args[1]  # just use this to pass in --no-cache or some such

$repo = "https://github.com/idaes/idaes-pse.git"
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
docker build --rm ${buildarg_1} --build-arg repo=${repo} --build-arg branch=${branch} -t ${flavor}_idaes_test .
Remove-Item extras -Recurse -Force -Confirm:$false
docker run --name ${flavor}_test_tmp -dt ${flavor}_idaes_test:latest
docker stop ${flavor}_test_tmp
docker rm ${flavor}_test_tmp
cd ..
