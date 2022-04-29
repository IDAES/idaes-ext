builddir=$1; shift

mkdir ${builddir}
cp -r ./scripts ${builddir}/
cp -r ./petsc ${builddir}/
cp -r ./src ${builddir}/
cp license.txt ${builddir}/
cp version.txt ${builddir}/
