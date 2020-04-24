$flavor = $args[0]
$buildarg_1 = $args[1]  # just use this to pass in --no-cache or some such

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

cd ${flavor}
docker build --rm ${buildarg_1} -t ${flavor}_idaes_test .
cd ..
