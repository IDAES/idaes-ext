flavor=$1
image="idaes-ext-${flavor}-test:latest"
cd "$flavor"
docker build -t "$image" .
cd ..
