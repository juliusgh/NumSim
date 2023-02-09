# NumSim

## Build the code

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

all simulations must be run from within the `build` directory.

### Convection

```shell
./scr/numsim ../parameters/convection.txt
```

### Lid driven Cavity

![ldc](/fig/lid_driven_cavity.png)

```shell
./scr/numsim ../parameters/lid_driven_cavity.txt
```

### Karman Vortex Street

![kvs](/fig/kvs_long.png)

```shell
./scr/numsim ../parameters/kvs.txt
```

### Karman Vortex Street with Heat

![kvsh](/fig/kvs_long.png)

```shell
./scr/numsim ../parameters/kvs_heat.txt
```

### Catamaran

![catamaran](/fig/catamaran_color_narrow.png)

Image:
By Obsidian Soul - Own work, CC0, <https://commons.wikimedia.org/w/index.php?curid=75764398.>

```shell
./scr/numsim ../parameters/catamaran.txt
```

### Catamaran with heat

![catamaran](/fig/catamaran_color_narrow.png)

Image:
By Obsidian Soul - Own work, CC0, <https://commons.wikimedia.org/w/index.php?curid=75764398>.

```shell
./scr/numsim ../parameters/catamaran_heat.txt
```

### Channel with Step

![cws](/fig/channel_with_step.png)

```shell
./scr/numsim ../parameters/channel_with_step.txt
```

### Nozzle

![nozzle](/fig/nozzle.png)

```shell
./scr/numsim ../parameters/nozzle.txt
```

### NumSim Obstacle

![numsimobstacle](/fig/numsim_obstacle.png)

```shell
./scr/numsim ../parameters/numsim_obstacle.txt
```

### SimTech Obstacle

![simtechobstacle](/fig/simtech_obstacle.png)

```shell
./scr/numsim ../parameters/numsim_obstacle.txt
```

### Two Inflows

![twochannels](/fig/two_channels.png)

```shell
./scr/numsim ../parameters/two_channels.txt
```

### Z Channel

![z](/fig/z.png)

```shell
./scr/numsim ../parameters/z.txt
```

## Memes

![meme1](/fig/meme1.jpeg)

![meme2](/fig/meme2.jpg)

![meme3](/fig/meme3.jpeg)

![meme4](/fig/meme4.png)
