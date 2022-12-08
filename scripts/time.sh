#!/bin/bash

module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

rm -rf build
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
make -j

cd ..

rm time.txt

for processes in 1 2 4 6 8 10 12 16
do
printf "\ntime $processes process(es):\n" >> time.txt
{ time mpirun -n $processes ./build/src/numsim_parallel parameters/lid_driven_cavity.txt ; } 2>> time.txt
done