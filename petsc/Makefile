ALL: petsc

CFLAGS = -I${ASL_INC}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

bindir = $(prefix)/bin

# Default do NOT set the static runtime; create the scaffolding
# for Windows which needs some more special handling
STATIC_RUNTIME =

petsc: petsc.o fg_nl.o fg_dae.o printing.o
	${CLINKER} ${STATIC_RUNTIME} -o petsc petsc.o fg_nl.o fg_dae.o printing.o $(ASL_LIB) ${PETSC_LIB}

py: petscpy
	${PYTHON} get_petsc_py.py

install: petsc
	${CP} petsc ${bindir}
