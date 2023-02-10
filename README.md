# NumSim: Towards a Real-World Fluid Solver

This is the code associated with the project of the Numerical Simulation lecture
at the University of Stuttgart developed during the winter semester 2022-2023.
The code was written by Julius Herb, Niklas Hornischer and Torben Schiz.

We decided to included the code parts we attempted to implement but could not finish,
due to time limitations or structural difficulties in combination with other simulation
components. These parts are namely `preconditioned CG`, `multigrid` and `free surface`.
We do this, to show that we spend a lot of time on them and did not make false statements
in our presentation of the our project.

## Build the code

To run the simulation, the project must be built first. We provided a debug
and a release version for your simulation.

The following targets are built:

- `numsim`: **contains the code of the project and runs in serial**
- `numsim_parallel`: contains only the basic functionality of the second submission and runs in parallel
- `tests`: contains test cases, which can be uncommented (only built if GTest installed)
- `compare_output`: compare vtk files against a reference solution

### Release

```shell
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

### Debug

```shell
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j
```

## Run Simulation Scenarios

We provided multiple example applications for our project. The following sections
provide instructions on how to run them. New domains can be created using the
`convert_image_to_domain.py` script. The script is located in `parameters/script`.
It takes images as inputs. Some lines have to be adapted accordingly to the situation.

All simulations must be run from within the `build` directory.

### Convection

```shell
./src/numsim ../parameters/convection.txt
```

### Lid driven Cavity

![ldc](/fig/lid_driven_cavity.png)

```shell
./src/numsim ../parameters/lid_driven_cavity.txt
```

### Karman Vortex Street

![kvs](/fig/kvs_long.png)

```shell
./src/numsim ../parameters/kvs.txt
```

### Karman Vortex Street with Heat

![kvsh](/fig/kvs_long.png)

```shell
./src/numsim ../parameters/kvs_heat.txt
```

### Catamaran

![catamaran](/fig/catamaran_color_narrow.png)

Image:
By Obsidian Soul - Own work, CC0, <https://commons.wikimedia.org/w/index.php?curid=75764398.>

```shell
./src/numsim ../parameters/catamaran.txt
```

### Catamaran with heat

![catamaran](/fig/catamaran_color_narrow.png)

Image:
By Obsidian Soul - Own work, CC0, <https://commons.wikimedia.org/w/index.php?curid=75764398>.

```shell
./src/numsim ../parameters/catamaran_heat.txt
```

### Channel with Step

![cws](/fig/channel_with_step.png)

```shell
./src/numsim ../parameters/channel_with_step.txt
```

### Nozzle

![nozzle](/fig/nozzle.png)

```shell
./src/numsim ../parameters/nozzle.txt
```

### Two Inflows

![twochannels](/fig/two_channels.png)

```shell
./src/numsim ../parameters/two_channels.txt
```

### Z Channel

![z](/fig/z.png)

```shell
./src/numsim ../parameters/z.txt
```

### NumSim Obstacle

![numsimobstacle](/fig/numsim_obstacle.png)

```shell
./src/numsim ../parameters/numsim_obstacle.txt
```

### SimTech Obstacle

![simtechobstacle](/fig/simtech_obstacle.png)

```shell
./src/numsim ../parameters/numsim_obstacle.txt
```

## Memes

![meme1](/fig/meme1.jpeg)

![meme2](/fig/meme2.jpg)

![meme3](/fig/meme3.jpeg)

![meme4](/fig/meme4.png)
