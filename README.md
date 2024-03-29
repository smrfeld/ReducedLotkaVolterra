# Code used in paper "Deep Learning Moment Closure Approximations using Dynamic Boltzmann Distributions"

This contains code used to generate figures in the paper "Deep Learning Moment Closure Approximations using Dynamic Boltzmann Distributions":

[arXiv 1905.12122](https://arxiv.org/abs/1905.12122)

## Requirements

The dependencies are automatically managed using the [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) package manager. **There is no need to manually install dependencies; just proceed to `Contents` below**.

For completeness, the dependencies downloaded automatically are:
* DynamicBoltzmann library v4.5 [here](https://github.com/smrfeld/DynamicBoltzmann/releases/tag/4.0).
* Q3 C1 Finite Elements library v3.0 [here](https://github.com/smrfeld/Q3-C1-Finite-Elements/releases/tag/3.0).
* LatticeGillespie C++ library v2.0 [here](https://github.com/smrfeld/LatticeGillespieCpp/releases/tag/2.0).
* Armadillo library v9.300.2 [here](http://arma.sourceforge.net/download.html).

## Contents

* [stoch_sims](stoch_sims) contains code to generate the stochastic simulations.
* [learn](learn) contains code to train the dynamic Boltzmann distribution.

Other (not used in paper) contents:
* `ode` is a Mathematica notebook to solve the ODE system for a well-mixed Lotka-Volterra system.
* `ssa` is a Mathematica notebook to generate stochastic simulations using the Gillespie algorithm for the well-mixed Lotka-Volterra system. This requires the Gillesipe module for Mathematica, available [here](https://github.com/smrfeld/GillespieMathematica).
