# Performance Optimizing and Parallelizing Semilinear PDE Solvers (I)

This problem will focus on implementing and optimizing the solution of the
2-dimensional Brusselator equations. The BRUSS equations are a well-known
highly stiff oscillatory system of partial differential equations which are
used in stiff ODE solver benchmarks. In this tutorial we will walk first
through a simple implementation, then do allocation-free implementations and
looking deep into solver options and benchmarking.

## Part 1: Implementing the BRUSS PDE System as ODEs

The Brusselator PDE is defined as follows:

$$\begin{align}
\frac{\partial u}{\partial t} &= 1 + u^2v - 4.4u + \alpha(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}) + f(x, y, t)\\
\frac{\partial v}{\partial t} &= 3.4u - u^2v + \alpha(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2})\end{align}$$

where

$$f(x, y, t) = \begin{cases}
5 & \quad \text{if } (x-0.3)^2+(y-0.6)^2 ≤ 0.1^2 \text{ and } t ≥ 1.1 \\
0 & \quad \text{else}\end{cases}$$

and the initial conditions are

$$\begin{align}
u(x, y, 0) &= 22\cdot y(1-y)^{3/2} \\
v(x, y, 0) &= 27\cdot x(1-x)^{3/2}\end{align}$$

with the periodic boundary condition

$$\begin{align}
u(x+1,y,t) &= u(x,y,t) \\
u(x,y+1,t) &= u(x,y,t)\end{align}$$

on a timespan of $t \in [0,22]$.

To solve this PDE, we will discretize it into a system of ODEs with the finite
difference method. We discretize `u` and `v` into arrays of the values at each
time point: `u[i,j] = u(i*dx,j*dy)` for some choice of `dx`/`dy`, and same for
`v`. Then our ODE is defined with `U[i,j,k] = [u v]`. The second derivative
operator, the Laplacian, discretizes to become a tridiagonal matrix with
`[1 -2 1]` and a `1` in the top right and bottom left corners. The nonlinear functions
are then applied at each point in space (they are broadcast). Use `dx=dy=1/32`.

You will know when you have the correct solution when you plot the solution
at `x=y=0` and see a periodic orbit, e.g., `ts=0:0.05:22; plot(ts, sol1.(ts,
idxs=1))`.

If you are not familiar with this process, see
[the Gierer-Meinhardt example from the SciMLTutorials.](http://tutorials.sciml.ai/html/introduction/03-optimizing_diffeq_code.html)

Note: Start by doing the simplest implementation!

## Part 2: Optimizing the BRUSS Code

PDEs are expensive to solve, and so we will go nowhere without some code
optimizing! Follow the steps described in the
[the Gierer-Meinhardt example from the SciMLTutorials](http://tutorials.sciml.ai/html/introduction/03-optimizing_diffeq_code.html)
to optimize your Brusselator code. Try other formulations and see what ends
up the fastest! Find a trade-off between performance and simplicity that suits
your needs.

## Part 3: Exploiting Jacobian Sparsity with Color Differentiation

Use the `sparsity!` function from [SparseDiffTools](https://github.com/JuliaDiffEq/SparseDiffTools.jl)
to generate the sparsity pattern for the Jacobian of this problem. Follow
the documentations [on the DiffEqFunction page](https://docs.sciml.ai/dev/features/performance_overloads)
to specify the sparsity pattern of the Jacobian. Generate an add the color
vector to speed up the computation of the Jacobian.

## (Optional) Part 4: Structured Jacobians

Specify the sparsity pattern using a BlockBandedMatrix from
[BlockBandedMatrices.jl](https://github.com/JuliaMatrices/BlockBandedMatrices.jl)
to accelerate the previous sparsity handling tricks.

## (Optional) Part 5: Automatic Symbolicification and Analytical Jacobian

Use the `modelingtoolkitize` function from ModelingToolkit.jl to convert your
numerical ODE function into a symbolic ODE function and use that to compute and
solve with an analytical sparse Jacobian.

## Part 6: Utilizing Preconditioned-GMRES Linear Solvers

Use the [linear solver specification page](https://docs.sciml.ai/dev/features/linear_nonlinear)
to solve the equation with `TRBDF2` with GMRES. Use the Sundials documentation
to solve the equation with `CVODE_BDF` with Sundials' special internal GMRES.
To both of these, use the [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
to add a preconditioner to the GMRES solver.

## Part 7: Exploring IMEX and Exponential Integrator Techniques (E)

Instead of using the standard `ODEProblem`, define a [`SplitODEProblem`](https://docs.sciml.ai/dev/types/split_ode_types)
to move some of the equation to the "non-stiff part". Try different splits
and solve with `KenCarp4` to see if the solution can be accelerated.

Next, use `MatrixFreeOperator` and `DiffEqArrayOperator` to define part of the equation as linear, and
use the `ETDRK4` exponential integrator to solve the equation. Note that this
technique is not appropriate for this equation since it relies on the
nonlinear term being non-stiff for best results.

## Part 8: Work-Precision Diagrams for Benchmarking Solver Choices

Use the `WorkPrecisionSet` method from
[DiffEqDevTools.jl](https://github.com/JuliaDiffEq/DiffEqDevTools.jl) to
benchmark multiple different solver methods and find out what combination is
most efficient.
[Take a look at DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl)
for usage examples.

## Part 9: GPU-Parallelism for PDEs (E)

Fully vectorize your implementation of the ODE and use a `CuArray` from
[CuArrays.jl](https://github.com/JuliaGPU/CuArrays.jl) as the initial condition
to cause the whole solution to be GPU accelerated.

## Part 10: Adjoint Sensitivity Analysis for Gradients of PDEs

In order to optimize the parameters of a PDE, you need to be able to compute
the gradient of the solution with respect to the parameters. This is done
through sensitivity analysis. For PDEs, generally the system is at a scale
where forward sensitivity analysis (forward-mode automatic differentiation)
is no longer suitable, and for these cases one uses adjoint sensitivity analysis.

Rewrite the PDE so the constant terms are parameters, and use the
[adjoint sensitivity analysis](https://docs.sciml.ai/latest/analysis/sensitivity/#Adjoint-Sensitivity-Analysis-1)
documentation to solve for the solution gradient with a cost function being the
L2 distance of the solution from the value 1. Solve with interpolated and
checkpointed adjoints. Play with using reverse-mode automatic differentiation
vs direct computation of vector-Jacobian products using the `autojacvec` option
of the `SensitivityAlg`. Find the set of options most suitable for this PDE.

If you have compute time, use this adjoint to optimize the parameters of the
PDE with respect to this cost function.
