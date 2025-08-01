# Stove Heat Simulation

This project demonstrates a minimal C++ program that reads a `stove.msh` mesh with the
[Gmsh](https://gmsh.info) library and computes a simple transient heat diffusion.
Four circular burners are kept at 200°C while the rest of the plate starts at
20°C. Results are written as Gmsh `.pos` files which can be animated directly in
the Gmsh GUI.

## Build

```bash
sudo apt-get install libgmsh-dev gmsh libeigen3-dev   # once
mkdir build && cd build
cmake ..
make
```

## Run

Copy the mesh next to the executable and run the solver:

```bash
cp ../src/stove.msh .
./tdk
```

This produces files `temperature_*.pos`. Open any of them in Gmsh and use
`Tools -> Animate` to view the temperature evolution.
