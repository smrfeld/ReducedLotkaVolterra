# Generate stochastic simulations of Lotka-Volterra system

## Contents

* `calculate_moments.cpp` - calculate higher order nearest neighbor (NN) and next-nearest neighbor (NNN) moments.
* `create_data_dirs.py` - create data directories.
* `main.cpp` - main script for running stochastic simulations.
* `mathematica` - Mathematica scripts
	* `figures` - generate the figures used in the paper. These appear into to the `figures` directory.
	* `viz_lattice` - misc. code for visualizing the lattice.

## Running the simulations

Compile:
```
g++ -std=c++14 -llattgillespie main.cpp -o main.o
```
and run:
```
./main.o
```
The output is in the `data` directories.

## Calculating extra moments

To calculate extra moments such as NNs and NNNs, compile:
```
g++ -std=c++14 -lbmla calculate_moments.cpp -o calculate_moments.o
```
and run:
```
./calculate_moments.o
```