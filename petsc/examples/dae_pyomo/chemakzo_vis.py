"""
Visualize the trajectory from example01.py. Assuming it was run with these
options:

    --ts_save_trajectory=1
    --ts_trajectory_type=visualization

That should put results from each time step in the directory Visualization-data

To run this you need either:

Option 1) $PETSC_DIR/lib/petsc/bin/ in PYTHONPATH

Option 2) Get these files and put them in your PYTHONPATH
  * https://github.com/petsc/petsc/blob/main/lib/petsc/bin/PetscBinaryIO.py
  * https://github.com/petsc/petsc/blob/main/lib/petsc/bin/PetscBinaryIOTrajectory.py
"""

import PetscBinaryIOTrajectory as pbt

if __name__ == '__main__':
    (t,v,names) = pbt.ReadTrajectory("Visualization-data")
    with open('vars.col') as f:
        names = list(map(str.strip, f.readlines()))
    with open('vars.typ') as f:
        typ = list(map(int,f.readlines()))
    names = [name for i, name in enumerate(names) if typ[i] in [0,1]]
    pbt.PlotTrajectories(t,v,names,["y[1]", "y[3]", "y[6]"])
    pbt.PlotTrajectories(t,v,names,["y[2]"])
    pbt.PlotTrajectories(t,v,names,["y[4]"])
    pbt.PlotTrajectories(t,v,names,["y[5]"])
