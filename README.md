# gamma-dde

Stan scripts for fitting ODE models to epidemic and viral dynamics data.
The method is described in the preprint

> T. Cassidy, P. Gillich, A.R. Humphries, and C.H. van Dorp, *Numerical methods and hypoexponential approximations for gamma distributed delay differential equations.* **todo: add link**

## Requirements

For using the Stan models, Stan version 2.4 or higher is required. 
See the [Stan website](www.mc-stan.org) for installation details. To get the latest version of Stan,
install `cmdstan` instad of `rstan` or `pystan`. The `cmdstan` interface for python can be installed with pip using

```bash
$ pip install cmdstanpy
```

and from python, `cmdstan` can be installed automatically

```py
import cmdstanpy
cmdstanpy.install_cmdstan()
```

For using the notebooks, [Jupyter](https://jupyter.org/) is required.


## Stan models

The folder `stan-models` contains the example SIR models.
Both the *fixed* and *smoothed* hypoexponential approximations are implemented.


## Notebooks

The `notebooks` folder contains a Jupyter notebook to simulate epidemic data and fit the Stan model
to this data.



_______________________________________________________________

FCI Open Source Copyright Assertion: C20017

copyright 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.
