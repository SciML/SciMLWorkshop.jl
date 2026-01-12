# RLC Circuit Acausal Model (I)

```@example rlc
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, @mtkmodel, @mtkcompile
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant
using DifferentialEquations, Plots

# RLC Circuit using standard library components
@mtkmodel RLCCircuit begin
    @structural_parameters begin
        R_val = 1.0
        L_val = 1.0
        C_val = 1.0
        V_val = 1.0
    end
    @components begin
        resistor = Resistor(R = R_val)
        inductor = Inductor(L = L_val)
        capacitor = Capacitor(C = C_val)
        source = Voltage()
        ground = Ground()
        constant = Constant(k = V_val)
    end
    @equations begin
        connect(constant.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, inductor.p)
        connect(inductor.n, capacitor.p)
        connect(capacitor.n, source.n)
        connect(capacitor.n, ground.g)
    end
end

# Create and compile the model
@mtkcompile sys = RLCCircuit()

# Create and solve the ODE problem with initial conditions
# Capacitor voltage starts at 0, inductor current starts at 0
prob = ODEProblem(sys, [sys.capacitor.v => 0.0, sys.inductor.i => 0.0], (0.0, 10.0))
sol = solve(prob)
plot(sol)
```
