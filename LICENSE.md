# IDAES Binary Distributions — Licenses and Third‑Party Notices

This document provides license and notice information for software
distributed in the IDAES binary tarballs:

- `idaes-functions-<platform>.tar.gz`
- `idaes-solvers-<platform>.tar.gz`

These tarballs are intended to be unpacked into the same installation
prefix so their `bin/`, `lib/`, `include/`, and `share/` directories can
be used together.

## Tarball contents (high level)

### `idaes-functions-<platform>.tar.gz`
Typically contains:
- `lib/` — IDAES "functions" shared libraries (platform extensions vary: `.so`, `.dylib`, `.dll`)
- `lib/helm_data/` — Helmholtz/EOS data and helpers (e.g., `*.json`, `*.py`, `*.nl`)
- `license_functions.txt` or `LICENSE.md` (this file)
- `version_functions.txt`

### `idaes-solvers-<platform>.tar.gz`
Typically contains:
- `bin/` — solver executables (e.g., `ipopt`, `bonmin`, `couenne`, `cbc`, `clp`, `petsc`, and platform variants like `.exe`)
- `lib/` — shared libraries and runtime dependencies (e.g., `libipopt*`, `libsipopt*`, `libpynumero_ASL*`, `petscpy/`)
- `include/` — headers (notably COIN-OR / Ipopt headers under `include/coin-or/`)
- `lib/pkgconfig/` — pkg-config metadata (e.g., `ipopt.pc`)
- `share/` — documentation installed by upstream projects (e.g., `share/doc/ipopt/`)
- `license.txt` or `LICENSE.md` (this file)
- `version_solvers.txt`

> **Platform variability**: Not all components are present on all platforms; some
> optional features (notably HSL and METIS) depend on availability, compatibility,
> and licensing constraints.

---

## IDAES components

### IDAES functions libraries and data
Examples (names vary by platform):
- `functions.*`
- `cubic_roots.*`
- `general_helmholtz_external.*`
- `lib/helm_data/*`

**License**: BSD 3-Clause (IDAES)
**IDAES license reference**: https://github.com/IDAES/idaes-pse/blob/main/LICENSE.md
**Source / build & packaging**: https://github.com/IDAES/idaes-ext

---

## COIN-OR solvers and related components

### Ipopt
- **License**: Eclipse Public License 2.0 (EPL-2.0)
- **Source**: https://github.com/coin-or/Ipopt

### Bonmin
- **License**: Eclipse Public License 1.0 (EPL-1.0)
- **Source**: https://github.com/coin-or/Bonmin

### Couenne
- **License**: Eclipse Public License 1.0 (EPL-1.0)
- **Source**: https://github.com/coin-or/Couenne

### Cbc
- **License**: Eclipse Public License 2.0 (EPL-2.0)
- **Source**: https://github.com/coin-or/Cbc

### Clp
- **License**: Eclipse Public License 2.0 (EPL-2.0)
- **Source**: https://github.com/coin-or/Clp

---

## AMPL Solver Library (ASL)

- **Source**: https://ampl.com/netlib/ampl/solvers/index.html
- **License**:
  Copyright (C) 2016 AMPL Optimization, Inc.; written by David M. Gay.
  
  Permission to use, copy, modify, and distribute this software and its
  documentation for any purpose and without fee is hereby granted, provided
  that the above copyright notice appear in all copies and that both that
  the copyright notice and this permission notice and warranty disclaimer
  appear in supporting documentation.
  
  The author and AMPL Optimization, Inc. disclaim all warranties with regard
  to this software, including all implied warranties of merchantability and
  fitness. In no event shall the author be liable for any special, indirect
  or consequential damages or any damages whatsoever resulting from loss of
  use, data or profits, whether in an action of contract, negligence or
  other tortious action, arising out of or in connection with the use or
  performance of this software.

---

## Pyomo / PyNumero

- **License**: BSD 3-Clause
- **Source**: https://github.com/Pyomo/pyomo

---

## PETSc

- **License**: BSD 2-Clause
https://petsc.org/main/install/license/#doc-license
- **PETSc Source**: https://petsc.org/

---

## MUMPS

- **License**: CeCILL-C
https://cecill.info/licences/Licence_CeCILL-C_V1-en.html
- **Source**: http://mumps.enseeiht.fr/

---

## METIS (platform-dependent)

- **METIS 5.x License**: Apache 2.0
https://github.com/KarypisLab/METIS/blob/master/LICENSE
**Source**: https://github.com/KarypisLab/METIS

> **Note**: METIS support is not included on macOS arm64 (Apple Silicon) in
> this distribution due to compatibility issues.

---

## HSL

- **Source**: https://www.hsl.rl.ac.uk/
- **License**: Proprietary (HSL terms apply)

**Acknowledgement required by HSL terms**:
> All technical papers, sales and publicity material resulting from use of the HSL codes within IPOPT must contain the following acknowledgement:
> "HSL, a collection of Fortran codes for large-scale scientific computation. See http://www.hsl.rl.ac.uk/."

---

## GCC Libraries

GCC runtime libraries are covered by the GCC runtime exception, so they can
be linked and distributed. The libraries were used unmodified.

- **Source**: https://gcc.gnu.org/
- **License**:
  - https://www.gnu.org/licenses/gpl-3.0.en.html
  - https://www.gnu.org/licenses/gcc-exception-3.1-faq.en.html

---

## MinGW

For Windows, MinGW runtime libraries are included.

**MinGW-w64**: https://www.mingw-w64.org/
**MSYS2 (common build environment)**: https://www.msys2.org/

---

## Additional notes

- Upstream projects may provide additional license texts and notices within the
  distribution (for example under `share/doc/`) and/or in their source repositories.
- This file is informational and does not replace or modify any upstream license terms.