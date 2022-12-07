#!/bin/bash

sbatch scripts/job.sh
rm -rf build
mkdir build
cd build
module load vtk/9.0.1
pwd
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
mpirun -n $1 src/numsim ../parameters/$2
