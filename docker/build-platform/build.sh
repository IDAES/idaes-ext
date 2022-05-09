flavor=$1

if [ "$flavor" = "windows" ]; then
  image="idaes-ext-windows-build:latest"
elif [ "$flavor" = "el7" ]; then
  image="idaes-ext-el7-build:latest"
elif [ "$flavor" = "el8" ]; then
  image="idaes-ext-el8-build:latest"
elif [ "$flavor" = "ubuntu1804" ]; then
  image="idaes-ext-ubuntu1804-build:latest"
elif [ "$flavor" = "ubuntu2004" ]; then
  image="idaes-ext-ubuntu2004-build:latest"
elif [ "$flavor" = "ubuntu2204" ]; then
  image="idaes-ext-ubuntu2204-build:latest"
else
  echo "Specify flavor in {el7, el8, ubuntu1804, ubuntu2004, ubuntu2204, windows}."
  exit 1
fi

cd "$flavor"
docker build -t "$image" .
cd ..
