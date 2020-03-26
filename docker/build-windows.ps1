$flavor = $args[0]

IF ($flavor -eq "windows"){
  $repo = "c:/repo"
} ELSEIF ($flavor -eq "centos7"){
  $repo = "/repo"
}
ELSE{
  echo "Please specify a flavor in {windows, centos7}."
  exit 1
}

cd ${flavor}
docker build --rm -t ${flavor}_build .
docker run --name ${flavor}_build_tmp -dt ${flavor}_build:latest
docker stop ${flavor}_build_tmp
docker cp ${flavor}_build_tmp:${repo}/idaes-ext/dist-lib/idaes-lib-${flavor}-64.tar.gz .
try{
  docker cp ${flavor}_build_tmp:${repo}/idaes-ext/dist-solvers/idaes-solvers-${flavor}-64.tar.gz .
}
catch{
  echo "Solvers were not built."
}
docker rm ${flavor}_build_tmp
cd ..
