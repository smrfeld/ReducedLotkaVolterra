# Learn moment closure model for the Lotka-Volterra system using dynamic Boltzmann distributions (DBDs)

Code used in the paper XXX.

## Requirements

* DynamicBoltzmann library v4.0 [here](https://github.com/smrfeld/DynamicBoltzmann/releases/tag/4.0).
* Q3 C1 Finite Elements library v3.0 [here](https://github.com/smrfeld/Q3-C1-Finite-Elements/releases/tag/3.0).
* LatticeGillespie C++ library v2.0 [here](https://github.com/smrfeld/LatticeGillespieCpp/releases/tag/2.0).
* Armadillo library v9.300.2 [here](http://arma.sourceforge.net/download.html).

## Contents

* [stoch_sims](stoch_sims) contains code to generate the stochastic simulations.
* [learn](learn) contains code to train the dynamic Boltzmann distribution.

Other contents:
* `ode` is a Mathematica notebook to solve the ODE system for a well-mixed Lotka-Volterra system.
* `ssa` is a Mathematica notebook to generate stochastic simulations using the Gillespie algorithm for the well-mixed Lotka-Volterra system. This requires the Gillesipe module for Mathematica, available [here](https://github.com/smrfeld/GillespieMathematica).