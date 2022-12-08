#!/bin/bash
#
#SBATCH --job-name=submission
#SBATCH --output=result.txt
#
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=4
#SBATCH --time=10:00

module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

cd build
cmake -DPROFILE=ON -DCMAKE_BUILD_TYPE=Release ..
make -j

cd ..

# enable custom build openmpi 3.1 that works with slurm
export CPATH=/scratch-nfs/maierbn/openmpi/install-3.1/include
export PATH=/scratch-nfs/maierbn/openmpi/install-3.1/bin:$PATH

srun -n $1 ./build/src/numsim_parallel parameters/$2

gprof ./build/src/numsim_parallel > ./build/out1.txt
gprof -l ./build/src/numsim_parallel > ./build/out2.txt

