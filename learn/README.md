# Learn moment closure model for Lotka-Volterra system

## Contents

* `create_data_dirs.py` - create data directories.
* `mathematica` - contains the following Mathematica notebooks:
	* `diagnose_sampling` - for visualizing Gibbs sampling steps.
	* `figures` - for generating the figures used in the paper.
	* `ixn_params_lp` - for lowpass filtering the learned interactions.
	* `monitor` - for monitoring progress during learning.
* `src` - contains the source code:
	* `learn_centered.cpp` and `dbm_centered.hpp` - these are the main codes used for training the moment closure model.
	* `sample_traj.cpp` and `dbm_static.hpp` - these are the codes used for testing the learned model by sampling at each timepoint.
	* `diagnose_sampling.cpp` - this is used to visualize the Gibbs sampling steps, if needed.

More directions follow.

## Building

Use the convenient Makefile:
```
mkdir build
cd build
cmake ..
```
then make the appropriate code as described below.

## Learning

This requires the stochastic simulations to have been generated as training data in the `stoch_sims` directory.

To learn the moment closure model, run from the build directory:
```
make learn_centered
cd ../bin
./learn_centered
```

The output will be in the `data/learn_centered` directories, specifically `moments`, `ixn_params`, and at the last optimization step `diff_eq_rhs`.

To visualize progress during learning, use the `mathematica/monitor` notebook.

## Lowpass filter

The interaction parameters can be lowpass filtered before sampling. Use the `mathematica/ixn_params_lp` notebook to do this.

## Sampling

To sample the learned model at each timepoint, run from the build directory:
```
make sample_traj
cd ../bin
./sample_traj
```
The output will be in the `data/learn_centered/sample_traj_lp` directory (or possibly `sample_traj` - make sure you specify in the code whether to read the lowpass filtered interaction parameters or the raw ones).

## Diagnosing the sampling

** Caution: this can produce a lot of data (sampled lattices) **

To diagnose the sampling, run from the build directory:
```
make diagnose_sampling
cd ../bin
./diagnose_sampling
```
The output will be in the `data/diagnose_sampling` directory. Use the `mathematice/diagnose_sampling` notebook to visualize the lattices, and compare the moments in the awake and asleep phases.

## Generating figures

Figures can be generated using the `mathematice/figures` notebook.