#!/usr/bin/env bash
# Build image idaes-ext-<flavor>:<arch>. Extra args are passed to `build`.
# Uses podman if available, else docker (override: CONTAINER_ENGINE).
# windows requires docker.

set -euo pipefail

FLAVORS="el8 el9 ubuntu2004 ubuntu2204 ubuntu2404 ubuntu2604 windows"

usage() {
  echo "Usage: $(basename "$0") <flavor> <arch> [extra build args...]"
  echo "  flavor: one of {${FLAVORS// /, }}"
  echo "  arch:   x86_64 or aarch64 (windows: x86_64 only)"
}

if [ $# -lt 2 ]; then
  usage
  exit 1
fi

flavor=$1
arch=$2
shift 2

case " $FLAVORS " in
  *" $flavor "*) ;;
  *) echo "Unknown flavor: $flavor" >&2; usage; exit 1 ;;
esac

case "$arch" in
  x86_64)  platform="linux/amd64" ;;
  aarch64) platform="linux/arm64" ;;
  *) echo "Unknown arch: $arch" >&2; usage; exit 1 ;;
esac

if [ "$flavor" = "windows" ]; then
  if [ "$arch" != "x86_64" ]; then
    echo "The windows flavor only supports x86_64." >&2
    exit 1
  fi
  if [ "${CONTAINER_ENGINE:-docker}" != "docker" ]; then
    echo "The windows flavor requires docker (CONTAINER_ENGINE=$CONTAINER_ENGINE)." >&2
    exit 1
  fi
  if ! command -v docker >/dev/null 2>&1; then
    echo "The windows flavor requires docker, which was not found." >&2
    exit 1
  fi
  engine=docker
  platform_args=""
else
  if [ -n "${CONTAINER_ENGINE:-}" ]; then
    engine=$CONTAINER_ENGINE
  elif command -v podman >/dev/null 2>&1; then
    engine=podman
  elif command -v docker >/dev/null 2>&1; then
    engine=docker
  else
    echo "Neither podman nor docker was found." >&2
    exit 1
  fi
  platform_args="--platform $platform"
fi

image="idaes-ext-${flavor}:${arch}"

echo "Building $image with $engine ${platform_args}"

cd "$(dirname "$0")/$flavor"
"$engine" build $platform_args -t "$image" ${1+"$@"} .
