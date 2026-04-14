flavor=$1

if [ "$flavor" = "windows" ]; then
  image="idaes-ext-windows-build:latest"
elif [ "$flavor" = "el8" ]; then
  image="idaes-ext-el8-build:latest"
elif [ "$flavor" = "el9" ]; then
  image="idaes-ext-el9-build:latest"
elif [ "$flavor" = "ubuntu2404" ]; then
  image="idaes-ext-ubuntu2404-build:latest"
elif [ "$flavor" = "ubuntu2004" ]; then
  image="idaes-ext-ubuntu2004-build:latest"
elif [ "$flavor" = "ubuntu2204" ]; then
  image="idaes-ext-ubuntu2204-build:latest"
else
  echo "Specify flavor in {el9, el8, ubuntu2404, ubuntu2004, ubuntu2204, windows}."
  exit 1
fi

cd "$flavor"
docker build -t "$image" .
cd ..
