# RLC Circuit Acausal Model (I)

```@example rlc
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant
using DifferentialEquations, Plots

# Create components with @named
@named resistor = Resistor(R = 1.0)
@named inductor = Inductor(L = 1.0)
@named capacitor = Capacitor(C = 1.0)
@named source = Voltage()
@named ground = Ground()
@named constant = Constant(k = 1.0)

# Define connections
eqs = [
    connect(constant.output, source.V)
    connect(source.p, resistor.p)
    connect(resistor.n, inductor.p)
    connect(inductor.n, capacitor.p)
    connect(capacitor.n, source.n)
    connect(capacitor.n, ground.g)
]

# Create and compile the system
@named sys = System(eqs, t; systems=[resistor, inductor, capacitor, source, ground, constant])
sys = mtkcompile(sys)

# Create and solve the ODE problem
prob = ODEProblem(sys, [capacitor.v => 0.0, inductor.i => 0.0], (0.0, 10.0))
sol = solve(prob)
plot(sol)
```
