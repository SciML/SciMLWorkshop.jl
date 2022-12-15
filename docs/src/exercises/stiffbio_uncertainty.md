# Investigating Sources of Randomness and Uncertainty in a Stiff Biological System (B)

In this problem we will walk through the basics of simulating models with
DifferentialEquations.jl. Let's take the
[Oregonator model of the Belousov-Zhabotinskii chemical reaction system](https://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Demos_Pitagora/DemoOrego/demoorego.pdf).
This system describes a classical example in non-equilibrium thermodynamics
and is a well-known natural chemical oscillator.

## Part 1: Simulating the Oregonator ODE model

When modeling, usually one starts off by investigating the deterministic model.
The deterministic ODE formulation of the Oregonator is
given by the equations

```math
\begin{align*}
\frac{dx}{dt} &= s(y-xy + x - qx^2)\\
\frac{dy}{dt} &= (-y - xy + z)/s\\
\frac{dz}{dt} &= w(x - z)
\end{align*}
```

with parameter values $s=77.27$, $w=0.161$, and $q=8.375 \times 10^{-6}$, and
initial conditions $x(0)=1$, $y(0)=2$, and $z(0)=3$. Use
[the tutorial on solving ODEs](https://docs.sciml.ai/DiffEqDocs/stable/getting_started/)
to solve this differential equation on the
timespan of $t\in[0,360]$ with the default ODE solver. To investigate the result,
plot the solution of all components over time, and plot the phase space plot of
the solution (hint: use `vars=(1,2,3)`). What shape is being drawn in phase space?

## Part 2: Investigating Stiffness

Because the reaction rates of `q` vs `s` is very large, this model has a "fast"
system and a "slow" system. This is typical of ODEs which exhibit a property
known as stiffness. Stiffness changes the ODE solvers which can handle the
equation well.
[Take a look at the ODE solver page](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
and investigate solving the equation using methods for non-stiff equations
(ex: `Tsit5`) and stiff equations (ex: `Rodas5`).

Benchmark using $t\in[0,50]$ using `@btime` from BenchmarkTools.jl. What
happens when you increase the timespan?

## (Optional) Part 3: Specifying Analytical Jacobians (I)

Stiff ODE solvers internally utilize the Jacobian of the ODE system in order
to improve the stepsizes in the solution. However, computing and factorizing
the Jacobian is costly, and thus it can be beneficial to provide the analytical
solution.

Use the
[ODEFunction definition page](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction)
to define an `ODEFunction` which holds both the OREGO ODE and its Jacobian, and solve using `Rodas5`.

## (Optional) Part 4: Automatic Symbolicification and Analytical Jacobian Calculations

Deriving Jacobians by hand is tedious. Thankfully symbolic mathematical systems
can do the work for you. And thankfully, DifferentialEquations.jl has tools
to automatically convert numerical problems into symbolic problems to perform
the analysis on!

follow the [ModelingToolkit.jl `modelingtoolkitize` tutorial](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/modelingtoolkitize/)
to automatically convert your ODE definition
to its symbolic form using `modelingtoolkitize` and calculate the analytical
Jacobian. Use the compilation functions to build the `ODEFunction` with the
embedded analytical solution.

## Part 5: Adding stochasticity with stochastic differential equations

How does this system react in the presence of stochasticity? We can investigate
this question by using stochastic differential equations. A stochastic
differential equation formulation of this model is known as the multiplicative
noise model, is created with:

```math
\begin{align*}
dx &= s(y-xy + x - qx^2)dt + \sigma_1 x dW_1\\
dy &= \frac{-y - xy + z}{s}dt + \sigma_2 y dW_2\\
dz &= w(x - z)dt + \sigma_3 z dW_3
\end{align*}
```

with $\sigma_i = 0.1$ where the `dW` terms describe a Brownian motion, a
continuous random process with normally distributed increments. Use the
[tutorial on solving SDEs](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/)
to solve simulate this model. Then,
[use the `EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/)
to generate and plot 100 trajectories of the stochastic model, and use
`EnsembleSummary` to plot the mean and 5%-95% region over time.

Try solving with the `ImplicitRKMil` and `SOSRI` methods. Notice that it isn't
stiff every single time!

(For fun, see if you can make the Euler-Maruyama `EM()` method solve this equation.
This requires a choice of `dt` small enough to be stable. This is the "standard"
method!)

## Part 6: Gillespie jump models of discrete stochasticity

When biological models have very few particles, continuous models no longer
make sense, and instead using the full discrete formulation can be required
to accuracy describe the dynamics. A discrete differential equation, or
Gillespie model, is a continuous-time Markov chain with Poisson-distributed
jumps. A discrete description of the Oregonator model is given by a chemical
reaction systems:

```
A+Y -> X+P
X+Y -> 2P
A+X -> 2X + 2Z
2X  -> A + P (note: this has rate kX^2!)
B + Z -> Y
```

where reactions take place at a rate which is proportional to its components,
i.e. the first reaction has a rate `k*A*Y` for some `k`.
Use the [tutorial on Gillespie SSA models](https://docs.sciml.ai/JumpProcesses/stable/tutorials/discrete_stochastic_example/)
to implement the `JumpProblem` for this model, and use the `EnsembleProblem`
and `EnsembleSummary` to characterize the stochastic trajectories.

For what rate constants does the model give the oscillatory dynamics for the
ODE approximation? For information on the true reaction rates, consult
[the original paper](https://pubs.acs.org/doi/abs/10.1021/ja00780a001).

## Part 7: Probabilistic Programming / Bayesian Parameter Estimation with DiffEqBayes.jl + Turing.jl (I)

In many cases, one comes to understand the proper values for their model's
parameters by utilizing data fitting techniques. In this case, we will use
the DiffEqBayes.jl library to perform a Bayesian estimation of the parameters.
For our data we will the following potential output:

```julia
t = 0.0:1.0:30.0
data = [1.0 2.05224 2.11422 2.1857 2.26827 2.3641 2.47618 2.60869 2.7677 2.96232 3.20711 3.52709 3.97005 4.64319 5.86202 9.29322 536.068 82388.9 57868.4 1.00399 1.00169 1.00117 1.00094 1.00082 1.00075 1.0007 1.00068 1.00066 1.00065 1.00065 1.00065
        2.0 1.9494 1.89645 1.84227 1.78727 1.73178 1.67601 1.62008 1.56402 1.50772 1.45094 1.39322 1.33366 1.2705 1.19958 1.10651 0.57194 0.180316 0.431409 251.774 591.754 857.464 1062.78 1219.05 1335.56 1419.88 1478.22 1515.63 1536.25 1543.45 1539.98
        3.0 2.82065 2.68703 2.58974 2.52405 2.48644 2.47449 2.48686 2.52337 2.58526 2.67563 2.80053 2.9713 3.21051 3.5712 4.23706 12.0266 14868.8 24987.8 23453.4 19202.2 15721.6 12872.0 10538.8 8628.66 7064.73 5784.29 4735.96 3877.66 3174.94 2599.6]
```

[Follow the exmaples on the parameter estimation page](https://docs.sciml.ai/DiffEqBayes/dev/examples/)
to perform a Bayesian parameter estimation. What are the most likely parameters
for the model given the posterior parameter distributions?

Use the `ODEProblem` to perform the fit. If you have time, use the `EnsembleProblem`
of `SDEProblem`s to perform a fit over averages of the SDE solutions. Note that
the SDE fit will take significantly more computational resources! See the GPU
parallelism section for details on how to accelerate this.

## (Optional) Part 8: Using Catalyst's Reaction Network DSL

Catalyst.jl is a helper library for the DifferentialEquations.jl
ecosystem for defining chemical reaction systems at a high level for easy
simulation in these various forms. Use the descrption
[from the Chemical Reaction Networks documentation page](https://docs.sciml.ai/Catalyst/stable/)
to build a reaction network and generate the ODE/SDE/jump equations, and
compare the result to your handcoded versions.
