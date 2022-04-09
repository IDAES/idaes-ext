flavor=$1
mname=$2
buildarg_1=$3  # just use this to pass in --no-cache or some such

repo="https://github.com/idaes/idaes-ext.git"
branch="main"

if [ "$flavor" = "windows" ]; then
  echo "Cannot build Windows binaries on Linux"
  exit 1
elif [ "$flavor" = "el7" ]; then
  wdir="/repo"
elif [ "$flavor" = "el8" ]; then
  wdir="/repo"
elif [ "$flavor" = "ubuntu1804" ]; then
  wdir="/repo"
elif [ "$flavor" = "ubuntu2004" ]; then
  wdir="/repo"
elif [ "$flavor" = "ubuntu2204" ]; then
  wdir="/repo"
else
  echo "Specify flavor in {el7, el8, ubuntu1804, ubuntu2004, ubuntu2204}."
  exit 1
fi

cp -r extras $flavor/extras/
cd $flavor

docker build --rm $buildarg_1 --build-arg repo=$repo --build-arg branch=$branch -t "$flavor"_build_itmp .
rm -rf extras
docker run --name "$flavor"_build_tmp -dt "$flavor"_build_itmp:latest
docker stop "$flavor"_build_tmp

docker cp "$flavor"_build_tmp:"$wdir"/idaes-ext/dist-lib/idaes-lib-"$flavor"-"$mname".tar.gz .
docker cp "$flavor"_build_tmp:"$wdir"/idaes-ext/dist-solvers/idaes-solvers-"$flavor"-"$mname".tar.gz .
docker cp "$flavor"_build_tmp:"$wdir"/idaes-ext/dist-petsc/idaes-petsc-"$flavor"-"$mname".tar.gz .

docker rm "$flavor"_build_tmp
docker rmi "$flavor"_build_itmp

mv *.tar.gz ../tarballs/
cd ..
