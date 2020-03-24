$flavor = $args[0]

IF ($flavor -eq "windows"){
  $repo = "c:/repo/idaes-ext/dist-lib"
} ELSEIF ($flavor -eq "centos7"){
  $repo = "/repo/idaes-ext/dist-lib"
}
ELSE{
  echo "Please specify a flavor in {windows, centos7}."
  exit 1
}

cd ${flavor}
docker build --rm -t ${flavor}_build .
docker run --name ${flavor}_build_tmp -dt ${flavor}_build:latest
docker stop ${flavor}_build_tmp
docker cp ${flavor}_build_tmp:${repo}/idaes-lib-${flavor}-64.tar.gz .
docker cp ${flavor}_build_tmp:${repo}/idaes-solvers-${flavor}-64.tar.gz .
docker rm ${flavor}_build_tmp
cd ..
