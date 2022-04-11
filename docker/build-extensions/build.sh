flavor=$1

repo="https://github.com/eslickj/idaes-ext.git"
branch="main"

if [ "$flavor" = "windows" ]; then
  echo "Cannot build Windows binaries on Linux"
  exit 1
elif [ "$flavor" = "el7" ]; then
  image="eslickj/idaes-ext-ubuntu-2004-build:latest"
  mname="x86_64"
  wdir="/repo"
elif [ "$flavor" = "el8" ]; then
  image="eslickj/idaes-ext-ubuntu-2004-build:latest"
  mname="x86_64"
  wdir="/repo"
elif [ "$flavor" = "ubuntu1804" ]; then
  image="eslickj/idaes-ext-ubuntu1804-build:latest"
  mname="x86_64"
  wdir="/repo"
elif [ "$flavor" = "ubuntu2004" ]; then
  image="eslickj/idaes-ext-ubuntu-2004-build:latest"
  mname="x86_64"
  wdir="/repo"
elif [ "$flavor" = "ubuntu2004_aarch64" ]; then
  image="ubuntu2004_build_aarch64:latest"
  flavor="ubuntu2004"
  mname="aarch64"
  wdir="/repo"
elif [ "$flavor" = "ubuntu2204" ]; then
  image="eslickj/idaes-ext-ubuntu-2204-build:latest"
  mname="x86_64"
  wdir="/repo"
else
  echo "Specify flavor in {el7, el8, ubuntu1804, ubuntu2004, ubuntu2204}."
  exit 1
fi


docker run --name "$flavor"_"$mname"_build_tmp -dt "$image"
docker cp ./extras/ "$flavor"_"$mname"_build_tmp:"$wdir"
docker exec "$flavor"_"$mname"_build_tmp sh -c "cp ${wdir}/extras/* ${wdir}"
docker exec "$flavor"_"$mname"_build_tmp sh -c "cd ${wdir}/extras && pwd"
docker exec "$flavor"_"$mname"_build_tmp sh -c "cd ${wdir} && git clone ${repo} && cd idaes-ext && git checkout ${branch}"
docker exec "$flavor"_"$mname"_build_tmp sh -c "cd ${wdir}/idaes-ext && bash scripts/compile_solvers.sh ${flavor}"
docker exec "$flavor"_"$mname"_build_tmp sh -c "cd ${wdir}/idaes-ext && bash scripts/compile_libs.sh ${flavor}"
docker stop "$flavor"_"$mname"_build_tmp

docker cp "$flavor"_"$mname"_build_tmp:"$wdir"/idaes-ext/dist-lib/idaes-lib-"$flavor"-"$mname".tar.gz .
docker cp "$flavor"_"$mname"_build_tmp:"$wdir"/idaes-ext/dist-solvers/idaes-solvers-"$flavor"-"$mname".tar.gz .
docker cp "$flavor"_"$mname"_build_tmp:"$wdir"/idaes-ext/dist-petsc/idaes-petsc-"$flavor"-"$mname".tar.gz .

docker rm /"$flavor"_"$mname"_build_tmp

mv *.tar.gz ./tarballs/
