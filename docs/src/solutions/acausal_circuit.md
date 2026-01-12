# RLC Circuit Acausal Model (I)

```@example rlc
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations, Plots

# Define the Pin connector
@connector Pin begin
    v(t)
    i(t), [connect = Flow]
end

# Ground component
@mtkmodel Ground begin
    @components begin
        g = Pin()
    end
    @equations begin
        g.v ~ 0
    end
end

# OnePort base component (two-pin device)
@mtkmodel OnePort begin
    @components begin
        p = Pin()
        n = Pin()
    end
    @variables begin
        v(t) = 0.0
        i(t) = 0.0
    end
    @equations begin
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
    end
end

# Resistor component
@mtkmodel Resistor begin
    @extend OnePort()
    @parameters begin
        R = 1.0
    end
    @equations begin
        v ~ i * R
    end
end

# Capacitor component
@mtkmodel Capacitor begin
    @extend OnePort()
    @parameters begin
        C = 1.0
    end
    @equations begin
        D(v) ~ i / C
    end
end

# Inductor component
@mtkmodel Inductor begin
    @extend OnePort()
    @parameters begin
        L = 1.0
    end
    @equations begin
        D(i) ~ v / L
    end
end

# Constant voltage source
@mtkmodel ConstantVoltage begin
    @extend OnePort()
    @parameters begin
        V = 1.0
    end
    @equations begin
        V ~ v
    end
end

# RLC Circuit model
@mtkmodel RLCCircuit begin
    @structural_parameters begin
        R = 1.0
        L = 1.0
        C = 1.0
        V = 1.0
    end
    @components begin
        resistor = Resistor(R = R)
        inductor = Inductor(L = L)
        capacitor = Capacitor(C = C)
        source = ConstantVoltage(V = V)
        ground = Ground()
    end
    @equations begin
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
