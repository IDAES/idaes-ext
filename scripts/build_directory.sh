builddir=$1; shift

mkdir ${builddir}
cp -r ./scripts ${builddir}/
cp -r ./petsc ${builddir}/
cp -r ./src ${builddir}/
cp LICENSE.md ${builddir}/
cp VERSION.md ${builddir}/
