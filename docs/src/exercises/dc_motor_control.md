# Controlling a DC Motor (E)

## Part 1 (Optional): Build the DC Motor System using the ModelingToolkit Standard Library

![](https://user-images.githubusercontent.com/1814174/206985282-d752cca3-cf3b-4551-922e-e99914050b04.png)

where $R = 0.5\Omega, L = 4.5\times 10^{-3} H, k = 0.5 N\cdot m/A, J = 0.02 \text{kg} \cdot m^2, f = 0.01 N\cdot m\cdot s/\text{rad}$. $R$ is the armature resistance, $L$ is the armature inductance, $k$ is the motor constant, $J$ is inertia, and $f$ is the friction factor.

To test the DC motor, hook it up to a constant voltage signal $5 V$ and a constant load $1 Nm$.

### Example Solution Plot

![](https://user-images.githubusercontent.com/1814174/206985435-dd80270f-7642-42ec-a4f8-04f503c8df5a.png)

# Part 2: Add a controller

Make use of the [Blocks.LimPI and Blocks.Feedback](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/#ModelingToolkitStandardLibrary.Blocks.LimPI) from ModelingToolkitStandardLibrary to close the loop around the motor model using a PI controller. The following is a block diagram of the construction:

```
                                 +-----+     +-----+
                            r    |     |  u  |     | y
                            --+->|  C  +---->|  P  +-+->
                             -|  |     |     |     | |
                              |  +-----+     +-----+ |
                              |                      |
                              +----------------------+
```

Simulate the system and plot the response to a unit step in the reference and, if you feel ambitious, a step in load disturbance appearing on the plant input as well.

## Part 3: Tune/detune the controller

ModelingToolkit supports inserting *analysis points* into models to facilitate *linear analysis*. An analysis point can be viewed as a named connection of particular interest for causal linear analysis. An example use case for analysis points is to derive [sensitivity functions](https://en.wikipedia.org/wiki/Sensitivity_(control_systems)). An analysis point is created automatically if a connection is named using the syntax `connect(output, :name, input)`.

The documentation for analysis points is [available here](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/#Linear-Analysis)

!!! note

    the analysis point connections are inherently *causal*, the order of the connected signals in an analysis-point connection is thus important. The wrong causality may cause a sensitivity function to turn into a complementary sensitivity function.

Let the linearized system model be described by the transfer function $P(s)$ and the linearized controller by $C(s)$, derive the sensitivity function:

```math
S(s) = \dfrac{1}{1 + P(s)C(s)}
```

at the plant output and plot the [Bode diagram](https://juliacontrol.github.io/ControlSystems.jl/latest/lib/plotting/#ControlSystemsBase.bodeplot).

To accomplish this, see functions:

* `using ControlSystemsBase`
* matrices, simplified_system = get_sensitivity(sys, analysis_point_name)`. This function is nice enough to linearize the system automatically.
* `ss` to create a linear state-space system. The matrices to create the state-space model are given by the call above.
* `bodeplot`

Try to tune the controller by changing the controller gain $k$ and the integration time $T$ in the PI-controller (`LimPI`) such that the closed-loop system becomes unstable. Watch how the sensitivity function changes as you reduce the stability margin. The documentation for the Block library of ModelingToolkitStandardLibrary is available [here](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/)

## Part 4 (Optional): Plot the Nyquist curve for the loop-transfer function

The loop-transfer function (at the input or output of the plant), may be computed using the function `get_looptransfer(sys, analysis_point_name)`. Create a linear state-space model representing the loop-transfer function in a similar fashion to how the sensitivity function was created above, and plot the Nyquist curve using the command `nyquistplot`.

### Example Solution Plot

With the PI-controller tuned like $k = 1.1, T = 0.05, u_{max} = 10, T_a = 0.035$, the plots should look like:

![](https://user-images.githubusercontent.com/1814174/206988013-0d4e7454-4e64-4448-9ee5-037bd1028149.png)
