# Building IDAES Binaries

The build scripts are not intended for general use. IDAES users do not
generally need to build IDAES binaries.

Where each platform is built:

| Platform | Host | How |
|---|---|---|
| Linux (x86_64, aarch64) | Apple Silicon Mac | podman or docker (x86_64 under emulation) |
| macOS (arm64) | Apple Silicon Mac | natively |
| Windows (x86_64) | Windows | docker, Windows containers |

Linux flavors: el8, el9, ubuntu2004, ubuntu2204, ubuntu2404, ubuntu2604.
Only Apple Silicon (arm64) is supported on macOS; Intel (x86_64) macOS builds
were dropped in 4.0.0.

## CoinHSL

If you have access to the official HSL source, package it as
`coinhsl.tar.gz`. Without it, pass `--without-hsl` (see
[Building without HSL](#building-without-hsl)).

## Linux

The scripts use podman if it is installed, otherwise docker. Set
`CONTAINER_ENGINE=docker` (or `podman`) to override. Make sure your podman
machine (or Docker Desktop) has enough memory and CPUs for the builds.

1. Build the image (once per flavor and arch):
   ```bash
   cd docker/main-images
   bash build-image.sh <flavor> <arch>   # e.g. bash build-image.sh ubuntu2604 aarch64
   ```
   This creates the image `idaes-ext-<flavor>:<arch>`. Extra arguments, such
   as `--no-cache`, are passed through to `build`.
2. Put `coinhsl.tar.gz` in `docker/build-scripts/extras/`.
3. Build the tarballs:
   ```bash
   cd docker/build-scripts
   bash build-solvers.sh <flavor> <arch> [repo] [branch] [image] [compile_solvers args...]
   ```
   `repo` defaults to `https://github.com/IDAES/idaes-ext.git`, `branch` to
   `main`, and `image` to `idaes-ext-<flavor>:<arch>`. Pass `""` to keep a
   default. The tarballs are written to `docker/build-scripts/tarballs/`.

**NOTE**: If you are behind certain corporate firewalls/proxies, you may need
to add your SSL certificates to the beginning of the Dockerfiles.

## Windows

Docker must be in Windows-container mode; see
[the official documentation](https://learn.microsoft.com/en-us/virtualization/windowscontainers/manage-docker/configure-docker-daemon).

1. Build the image from Git Bash:
   ```bash
   cd docker/main-images
   bash build-image.sh windows x86_64
   ```
2. Put `coinhsl.tar.gz` in `docker/build-scripts/extras/`.
3. Build the tarballs from PowerShell:
   ```powershell
   cd docker\build-scripts
   .\build.ps1 [-Repo <url>] [-Branch <branch>] [-Image <image>] [-SolverArgs "--without-hsl"] [-Memory 8g]
   ```
   The tarballs are written to `docker\build-scripts\tarballs\`.

## macOS

### Environment

1. Install the Xcode command line tools.
2. Install Homebrew and the build dependencies. Homebrew's bash is required
   because the macOS default is too old for coinbrew:
   ```bash
   brew install cmake gcc pkgconfig boost bash wget
   ```
   `compile_solvers.sh` uses `HOMEBREW_PREFIX` (set by `brew shellenv`) to
   find gcc, g++, and gfortran.
3. Build PETSc. The newest version that works with the code in this
   repository is 3.20.6.
   ```bash
   mkdir -p ~/repo && cd ~/repo
   wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.20.6.tar.gz
   tar -zxf petsc-3.20.6.tar.gz && mv petsc-3.20.6 petsc && cd petsc
   ./configure --with-debug=0 --with-shared=0 --with-mpi=0 --with-fortran-bindings=0 \
       --download-metis --download-mumps --with-mumps-serial=1 --prefix=$HOME/repo/petsc-dist \
       --download-metis-cmake-arguments='-DCMAKE_POLICY_VERSION_MINIMUM=3.5'
   make && make install
   ```
   **NOTE**: If configure has trouble downloading MUMPS, download it yourself
   and use `--download-mumps=../MUMPS_5.6.1.tar.gz`.
4. `export PETSC_DIR=$HOME/repo/petsc-dist` (the script's default is
   `$HOME/src/petsc-dist`).

### Building

1. From the idaes-ext checkout, copy the sources to a build directory:
   ```bash
   bash ./scripts/build_directory.sh ~/repo/ext-build
   ```
2. Put `coinhsl.tar.gz` in the *parent* of the build directory
   (`~/repo/coinhsl.tar.gz` here).
3. Build:
   ```bash
   cd ~/repo/ext-build
   bash ./scripts/compile_solvers.sh darwin
   bash ./scripts/compile_libs.sh darwin
   ```
   The tarballs are `dist/idaes-solvers-darwin-arm64.tar.gz` and
   `dist-functions/idaes-functions-darwin-arm64.tar.gz`.

## Building without HSL

For testing, pass `--without-hsl` to `compile_solvers.sh`. HSL-dependent
solvers (`k_aug`, `dot_sens`) are skipped. For example:

```bash
bash build-solvers.sh ubuntu2204 x86_64 "" "" "" --without-hsl   # Linux
.\build.ps1 -SolverArgs "--without-hsl"                       # Windows
bash ./scripts/compile_solvers.sh darwin --without-hsl        # macOS
```

## Releasing

1. Set the version in `VERSION.md`.
2. Build all platforms and collect the tarballs in one directory.
3. Generate the hash file and copy it into `releases/`:
   ```bash
   bash scripts/hash_extensions.sh <release> <directory>
   ```
   idaes-pse downloads this file to verify the binaries.
4. Publish the tarballs as a GitHub release, and set `RELEASE_REPO` and
   `RELEASE_TAG` in `.github/workflows/test_binaries.yml` to test them.
