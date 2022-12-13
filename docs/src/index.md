# SciMLWorkshop: Workshop Materials for Training in Scientific Computing and Scientific Machine Learning (SciML)

These exercises teach common workflows which involve SciML's tools like
DifferentialEquations.jl, ModelingToolkit.jl DiffEqFlux.jl, and the connections to parts
like stochastic differential equations and Bayesian estimation.

## Difficulty Levels

The designation (B) is for "Beginner", meaning that a user new to the package
should feel comfortable trying this exercise. An exercise designated (I) is
for "Intermediate", meaning the user may want to have some previous background
in DifferentialEquations.jl or try some (B) exercises first. The additional
(E) designation is for "Experienced", which are portions of exercises which may
take some work.

## Overview

The exercises are described as follows:

- Exercise 1 takes the user through solving a stiff ordinary differential equation
  and using the ModelingToolkit.jl to automatically convert the function to a
  symbolic form to derive the analytical Jacobian to speed up the solver. The
  same biological system is then solved with stochasticity, utilizing
  EnsembleProblems to understand 95% bounds on the solution. Finally,
  probabilistic programming is employed to perform Bayesian parameter estimation
  of the parameters against data.
- Exercise 2 takes the user through defining hybrid delay differential equation,
  that is a differential equation with events, and using differentiable programming
  techniques (automatic differentiation) to to perform gradient-based parameter
  estimation.
- Exercise 3 takes the user through differential-algebraic equation (DAE)
  modeling, the concept of index, and using both mass-matrix and implicit
  ODE representations. This will require doing a bit of math, but the student
  will understand how to change their equations to make their DAE numerically
  easier for the integrators.
- Exercise 4 takes the user through optimizing a PDE solver, utilizing
  automatic sparsity pattern recognition, automatic conversion of numerical
  codes to symbolic codes for analytical construction of the Jacobian,
  preconditioned GMRES, and setting up a solver for IMEX and GPUs, and compute
  adjoints of PDEs.
- Exercise 5 focuses on a chaotic orbit, utilizing parallel ensembles across
  supercomputers and GPUs to quickly describe phase space.
- Exercise 6 takes the user through training a neural stochastic differential
  equation, using GPU-accleration and adjoints through Flux.jl's neural
  network framework to build efficient training codes.

This exercise worksheet is meant to be a living document leading new users through
a deep dive of the SciML feature set. If you further suggestions
or want to contribute new problems, please 
[open an issue or PR at the SciMLWorkshop.jl repository](https://github.com/SciML/SciMLWorkshop.jl).

## Reproducibility
```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```
```@example
using Pkg # hide
Pkg.status() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>and using this machine and Julia version.</summary>
```
```@example
using InteractiveUtils # hide
versioninfo() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```
```@example
using Pkg # hide
Pkg.status(;mode = PKGMODE_MANIFEST) # hide
```
```@raw html
</details>
```
```@raw html
You can also download the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Manifest.toml"
```
```@raw html
">manifest</a> file and the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Project.toml"
```
```@raw html
">project</a> file.
```
