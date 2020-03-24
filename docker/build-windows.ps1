$flavor = "windows"

IF ($flavor -eq "windows"){
  $repo = "c:/repo/idaes-ext/dist-lib"
} ELSE{
  $repo = "/repo/idaes-ext/dist-lib"
}

cd ${flavor}
docker build --rm -t ${flavor}_build .
docker run --name ${flavor}_build_tmp -dt ${flavor}_build:latest
docker stop ${flavor}_build_tmp
docker cp ${flavor}_build_tmp:${repo}/idaes-lib-${flavor}-64.tar.gz .
docker rm ${flavor}_build_tmp
cd ..
