#!/usr/bin/env bash

flavor=$1
arch=$2
tag=$3

if [ -z "$tag" ]; then
  tag="latest"
fi

if [ -z "$flavor" ] || [ -z "$arch" ]; then
  echo "Usage: $0 <flavor> <arch> [tag]"
  echo "  flavor: {el9, el8, ubuntu2404, ubuntu2004, ubuntu2204, windows}"
  echo "  arch:   {x86_64, aarch64}"
  echo "  tag:    docker tag value (optional, defaults to latest)"
  exit 1
fi

if [ "$flavor" = "windows" ]; then
  image="idaes-ext-${flavor}:${tag}"
elif [ "$flavor" = "el8" ]; then
  image="idaes-ext-${flavor}:${tag}"
elif [ "$flavor" = "el9" ]; then
  image="idaes-ext-${flavor}:${tag}"
elif [ "$flavor" = "ubuntu2404" ]; then
  image="idaes-ext-${flavor}:${tag}"
elif [ "$flavor" = "ubuntu2004" ]; then
  image="idaes-ext-${flavor}:${tag}"
elif [ "$flavor" = "ubuntu2204" ]; then
  image="idaes-ext-${flavor}:${tag}"
else
  echo "Specify flavor in {el9, el8, ubuntu2404, ubuntu2004, ubuntu2204, windows}."
  exit 1
fi

if [ "$arch" = "x86_64" ]; then
  platform="linux/amd64"
elif [ "$arch" = "aarch64" ]; then
  platform="linux/arm64"
else
  echo "Specify arch as one of {x86_64, aarch64}."
  exit 1
fi

cd "$flavor" || exit 1

if [ "$flavor" = "windows" ]; then
  docker build -t "$image" .
else
  docker build --platform="$platform" -t "$image" .
fi

cd ..
