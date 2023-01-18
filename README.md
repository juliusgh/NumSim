# NumSim

Build the code and run the simulation using MPI on `n=4` cores:

```shell
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
mpirun -n 4 ./src/numsim_parallel ../parameters/lid_driven_cavity.txt
```
