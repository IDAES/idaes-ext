# IDAES Extensions

This repository hosts IDAES binary extensions.  These extensions include solvers and function libraries.  These binary extensions are intended for use with the IDAES Prosess Modeling Framework https://github.com/IDAES/idaes-pse. 

## Getting Extensions

If you have the IDAES framework installed, you can get the extensions by running the following command:

```sh
idaes get-extensions
```

You can also download a particular release of the extensions by providing a URL:

```sh
idaes get-extensions --url <url>
```

If you wish to get a particular version for the idaes extensions

## Contents

The extension contain a version of the IPOPT solver compiled the HSL linear solver library, for which IDAES has obtained 
a distribution license. All technical papers, sales and publicity material resulting from use of the HSL codes within IPOPT 
must contain the following acknowledgement: HSL, a collection of Fortran codes for large-scale scientific computation. See http://www.hsl.rl.ac.uk.

See https://github.com/IDAES/idaes-ext/blob/master/license.txt for additional information on thrird-party code contained 
in the binaries.

The binries also include property librairy functions compiled from source at https://github.com/IDAES/idaes-pse/tree/master/src.
