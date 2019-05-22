# Learn moment closure model for Lotka-Volterra system

## Contents

* `create_data_dirs.py` - create data directories.
* `mathematica` contains Mathematica scripts.
* ``

## Running

This requires the stochastic simulations to have been generated as training data in the `stoch_sims` directory.

Use the convenient Makefile:
```
mkdir build
cd build
cmake ..
```
then make the appropriate code as described below.

## Learning

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