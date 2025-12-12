# README

This repository contains scripts and source code for compiling IDAES extensions.
These consist of IDAES property libraries and third party solvers. IDAES users
generally do not need to access this repo. The build scripts are specific
to the IDAES build environment and not intended for general use.

## Documentation

- More detailed build docs: [doc/build.md](doc/build.md)
- For more information on building binaries: [doc/build.md](doc/build.md)
- For more information about the contents of binary releases:
  https://idaes-pse.readthedocs.io/en/stable/tutorials/getting_started/binaries.html#binary-packages
- For installation instructions:
  https://idaes-pse.readthedocs.io/en/stable/tutorials/getting_started/index.html

## Platforms

### Build targets

We currently build binaries for:

- Linux
  - Distributions: el8, el9, ubuntu2004, ubuntu2204, ubuntu2404
  - Architectures: x86_64 and aarch64 (for all listed Linux distributions)
- Windows
  - Architecture: x86_64
- macOS
  - Architectures: Intel (x86_64) and Apple Silicon (arm64)
  - Note: Metis is **not** built/used on macOS arm64 due to compatibility issues.

### Test targets

We currently run tests against:

- Ubuntu: 20.04, 22.04, 24.04
- macOS: Intel (x86_64), Apple Silicon (arm64)
- Windows: 2022, 2025
- Rocky Linux: 8, 9
- Debian: 11, 12, 13
