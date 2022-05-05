flavor=$1
image="idaes-ext-${flavor}-test:latest"
cd "$flavor"
docker build --no-cache -t "$image" .
cd ..
