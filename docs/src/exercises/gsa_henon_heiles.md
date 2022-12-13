# Global Parameter Sensitivity and Optimality with GPU and Distributed Ensembles (B)

In this example we will investigate how the parameters "generally" effect the
solution in the chaotic Henon-Heiles system. By "generally" we will use global
sensitivity analysis methods to get an average global characterization of the
parameters on the solution. In addition to a global sensitivity approach, we
will generate large ensembles of solutions with different parameters using
a GPU-based parallelism approach.

## Part 1: Implementing the Henon-Heiles System (B)

The Henon-Heiles Hamiltonian system is described by the ODEs:

```math
\begin{align*}
\frac{dp_1}{dt} &= -q_1 (1 + 2q_2)\\
\frac{dp_2}{dt} &= -q_2 - (q_1^2 - q_2^2)\\
\frac{dq_1}{dt} &= p_1\\
\frac{dq_2}{dt} &= p_2
\end{align*}
```

with initial conditions $u_0 = [0.1,0.0,0.0,0.5]$.
Solve this system over the timespan $t\in[0,1000]$

## (Optional) Part 2: Alternative Dynamical Implmentations of Henon-Heiles (B)

The Henon-Heiles defines a Hamiltonian system with certain structures which
can be utilized for a more efficient solution. Use [the Dynamical problems page](https://docs.sciml.ai/dev/types/dynamical_types)
to define a `SecondOrderODEProblem` corresponding to the acceleration terms:

```math
\begin{align*}
\frac{d^2q_1}{dt^2} &= -q_1 (1 + 2q_2)\\
\frac{d^2q_2}{dt^2} &= -q_2 - (q_1^2 - q_2^2)
\end{align*}
```

Solve this with a method that is specific to dynamical problems, like `DPRKN6`.

The Hamiltonian can also be directly described:

```math
H(p,q) = \frac{1}{2}(p_1^2 + p_2^2) + \frac{1}{2}(q_1^2+q_2^2+2q_1^2 q_2 - \frac{2}{3}q_2^3)
```

Solve this problem using the `HamiltonianProblem` constructor from DiffEqPhysics.jl.

## Part 3: Parallelized Ensemble Solving

To understand the orbits of the Henon-Heiles system, it can be useful to solve
the system with many different initial conditions. Use the
[ensemble interface](https://docs.sciml.ai/dev/features/ensemble)
to solve with randomized initial conditions in parallel using threads with
`EnsembleThreads()`. Then, use `addprocs()` to add more cores and solve using
`EnsembleDistributed()`. The former will solve using all of the cores on a
single computer, while the latter will use all of the cores on which there
are processors, which can include thousands across a supercomputer! See
[Julia's parallel computing setup page](https://docs.julialang.org/en/v1/manual/parallel-computing/index.html)
for more details on the setup.

## Part 4: Parallelized GPU Ensemble Solving

Setup the CUDAnative.jl library and use the `EnsembleGPUArray()` method to
parallelize the solution across the thousands of cores of a GPU. Note that
this will efficiency solve for hundreds of thousands of trajectores.
