# Fitting Hybrid Delay Pharmacokinetic Models with Automated Responses (B)

Hybrid differential equations are differential equations with events, where
events are some interaction that occurs according to a prespecified condition.
For example, the bouncing ball is a classic hybrid differential equation given
by an ODE (Newton's Law of Gravity) mixed with the fact that, whenever the
ball hits the floor (`x=0`), then the velocity of the ball flips (`v=-v`).

In addition, many models incorporate delays, that is the driving force of the
equation is dependent not on the current values, but values from the past.
These delay differential equations model how individuals in the economy act
on old information, or that biological processes take time to adapt to a new
environment.

In this equation we will build a hybrid delayed pharmacokinetic model and
use the parameter estimation techniques to fit this it to a data.

## Part 1: Defining an ODE with Predetermined Doses

First, let's define the simplest hybrid ordinary differential equation: an ODE
where the events take place at fixed times. The ODE we will use is known as
the one-compartment model:

```math
\begin{align*}
\frac{d[Depot]}{dt} &= -K_a [Depot] + R\\
\frac{d[Central]}{dt} &= K_a [Depot] - K_e [Central]
\end{align*}
```

with $t \in [0,90]$, $u_0 = [100.0,0]$, and $p=[K_a,K_e]=[2.268,0.07398]$.

With this model, use [the event handling documentation page](https://docs.sciml.ai/dev/features/callback_functions)
to define a `DiscreteCallback` which fires at `t ∈ [24,48,72]` and adds a
dose of 100 into `[Depot]`. (Hint: you'll want to set `tstops=[24,48,72]` to
force the ODE solver to step at these times).

## Part 2: Adding Delays

Now let's assume that instead of there being one compartment, there are many
transit compartment that the drug must move through in order to reach the
central compartment. This effectively delays the effect of the transition from
`[Depot]` to `[Central]`. To model this effect, we will use the delay
differential equation which utilizes a fixed time delay $\tau$:

```math
\begin{align*}
\frac{d[Depot]}{dt} &= -K_a [Depot](t)\\
\frac{d[Central]}{dt} &= K_a [Depot](t-\tau) - K_e [Central]
\end{align*}
```

where the parameter $τ = 6.0$.
[Use the DDE tutorial](https://docs.sciml.ai/dev/tutorials/dde_example)
to define and solve this delayed version of the hybrid model.

## Part 3: Automatic Differentiation (AD) for Optimization (I)

In order to fit parameters $(K_a,K_e,\tau)$ we will want to be able to calculate
the gradient of the solution with respect to the initial conditions. One way to
do this is via Automatic Differentiation (AD). For small numbers of parameters
(<100), it is fastest to use Forward-Mode Automatic Differentiation
(even faster than using adjoint sensitivity analysis!). Thus for this problem
we will make use of ForwardDiff.jl to use Dual number arithmetic to retrieve
both the solution and its derivative w.r.t. parameters in a single solve.

[Use the information from the page on local sensitvity analysis](https://docs.sciml.ai/dev/analysis/sensitivity)
to define the input dual numbers, solve the equation, and plot both the solution
over time and the derivative of the solution w.r.t. the parameters.

## Part 4: Fitting Known Quantities with DiffEqParamEstim.jl + Optim.jl

Now let's fit the delayed model to a dataset. For the data, use the array

```julia
t = 0.0:12.0:90.0
data = [100.0 0.246196 0.000597933 0.24547 0.000596251 0.245275 0.000595453 0.245511
        0.0 53.7939 16.8784 58.7789 18.3777 59.1879 18.5003 59.2611]
```

Use [the parameter estimation page](https://docs.sciml.ai/dev/analysis/parameter_estimation)
to define a loss function with `build_loss_objective` and optimize the parameters
against the data. What parameters were used to generate the data?

## Part 5: Implementing Control-Based Logic with ContinuousCallbacks (I)

Now that we have fit our delay differential equation model to the dataset, we
want to start testing out automated treatment strategies. Let's assume that
instead of giving doses at fixed time points, we invent a wearable which
monitors the patient and administers a dose whenever the internal drug
concentration falls below 25. To model this effect, we will need to use
`ContinuousCallbacks` to define a callback that triggers when `[Central]` falls
below the threshold value.

[Use the documentation on the event handling page](https://docs.sciml.ai/dev/features/callback_functions) to define such a callback,
and plot the solution over time. How many times does the auto-doser administer
a dose? How much does this change as you change the delay time $\tau$?

## Part 6: Global Sensitivity Analysis with the Morris and Sobol Methods

To understand how the parameters effect the solution in a global sense, one
wants to use Global Sensitivity Analysis. Use the
[GSA documentation page](https://docs.sciml.ai/dev/analysis/global_sensitivity)
perform global sensitivity analysis and quantify the effect of the various
parameters on the solution.
