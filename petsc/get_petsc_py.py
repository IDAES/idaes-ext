#!/usr/bin/env python3

import os
import sys
import json
import shutil

petsc_dir = os.environ["PETSC_DIR"]
petsc_py = os.path.join(petsc_dir, "lib/petsc/bin/")

sys.path.append(petsc_py)

from petsc_conf import get_conf

if __name__ == "__main__":
    c = get_conf()
    with open(os.path.join("petscpy","petsc_conf.json"), "w") as f:
        json.dump(list(c), f)
    for f in ["PetscBinaryIO.py", "PetscBinaryIOTrajectory.py"]:
        shutil.copyfile(os.path.join(petsc_py, f), os.path.join(os.getcwd(), "petscpy", f))
