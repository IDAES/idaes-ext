$flavor = $args[0]
$buildarg_1 = $args[1]  # just use this to pass in --no-cache or some such

# The 3rd and 4th arguments provided will be interpreted as repo and branch.
# (If you don't want to use buildarg_1, just pass in an empty string.)
# The 5th argument is for the tag of the Docker image. It will default
# to `latest` if not set
$repo = $args[2]
$branch = $args[3]
$tag = $args[4]

# If repo and branch are not provided, use default values
IF ($repo -eq $null){
  $repo = "https://github.com/idaes/idaes-ext.git"
}
IF ($branch -eq $null){
  $branch = "main"
}
IF ($tag -eq $null){
  $tag = "latest"
}

$mname = "x86_64"

IF ($flavor -eq "windows"){
  $wdir = "c:/repo"
}
ELSEIF ($flavor -eq "el9"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "el8"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "ubuntu2404"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "ubuntu2004"){
  $wdir = "/repo"
}
ELSEIF ($flavor -eq "ubuntu2204"){
  $wdir = "/repo"
}
ELSE{
  echo "Specify flavor in {windows, el8, el9, ubuntu2004, ubuntu2204, ubuntu2404}."
  exit 1
}

xcopy /y /E extras ${flavor}\extras\

cd ${flavor}
docker build --rm ${buildarg_1} --build-arg repo=${repo} --build-arg branch=${branch} -t ${flavor}_build_itmp .
Remove-Item extras -Recurse -Force -Confirm:$false
docker run --name ${flavor}_build_tmp -dt ${flavor}_build_itmp:latest
docker stop ${flavor}_build_tmp
docker cp ${flavor}_build_tmp:${wdir}/idaes-ext/dist-functions/idaes-functions-${flavor}-${mname}.tar.gz .
try{
  docker cp ${flavor}_build_tmp:${wdir}/idaes-ext/dist/idaes-solvers-${flavor}-${mname}.tar.gz .
}
catch{
  echo "Solvers were not built."
}
docker rm ${flavor}_build_tmp
docker rmi ${flavor}_build_itmp
cp *.tar.gz ../tarballs/
Remove-Item *.tar.gz
cd ..
