# RLC Circuit Acausal Model (I)

```@example rlc
module RLC

using ModelingToolkit

@variables t
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name=name)
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], []; name=name), g)
end

function OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t)=1.0 i(t)=1.0
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
          ]
    compose(ODESystem(eqs, t, sts, []; name=name), p, n)
end

function Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function Capacitor(;name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
           D(v) ~ i / C
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function Inductor(; name, L = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters L = L
    D = Differential(t)
    eqs = [
        D(i) ~ v / L,
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

function ConstantVoltage(;name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V
    eqs = [
           V ~ v
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

end

using .RLC
using ModelingToolkit, Plots, DifferentialEquations
L = C = V = R = 1.0
@named inductor = RLC.Inductor(L=L)
@named resistor = RLC.Resistor(R=R)
@named capacitor = RLC.Capacitor(C=C)
@named source = RLC.ConstantVoltage(V=V)
@named ground = RLC.Ground()

rlc_eqs = [
          connect(source.p, resistor.p)
          connect(resistor.n, inductor.p)
          connect(inductor.n, capacitor.p)
          connect(capacitor.n, source.n)
          connect(capacitor.n, ground.g)
         ]

@named _rlc_model = ODESystem(rlc_eqs, RLC.t)
@named rlc_model = compose(_rlc_model,
                          [inductor, resistor, capacitor, source, ground])
sys = structural_simplify(rlc_model)
u0 = [
      capacitor.v => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob)
plot(sol)
```
