"""
Visualize the trajectory from example01.py. Assuming it was run with these
options:

    --ts_save_trajectory=1
    --ts_trajectory_type=visualization

That should put results from each time step in the directory Visualization-data

To run this you need either:

Option 1) $PETSC_DIR/lib/petsc/bin/ in PYTHONPATH

Option 2) Get these files and put them in PYTHONPATH
  * https://github.com/petsc/petsc/blob/main/lib/petsc/bin/PetscBinaryIO.py
  * https://github.com/petsc/petsc/blob/main/lib/petsc/bin/PetscBinaryIOTrajectory.py
"""
import sys
import PetscBinaryIOTrajectory as pbt

def main():
    try:
        stub = sys.argv[1]
    except:
        print("Specify nl file stub")
        return
    (t,v,names) = pbt.ReadTrajectory("Visualization-data")
    with open(f'{stub}.col') as f:
        names = list(map(str.strip, f.readlines()))
    with open(f'{stub}.typ') as f:
        typ = list(map(int,f.readlines()))
    names = [name for i, name in enumerate(names) if typ[i] in [0,1]]
    skip = ["xdummy", "timedummy"]
    for name in names:
        if name in skip: continue
        pbt.PlotTrajectories(t,v,names,[name])

if __name__ == '__main__':
    main()
