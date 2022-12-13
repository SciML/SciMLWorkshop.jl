# Controlling a DC Motor (A)

```@example dc_motor
using ModelingToolkit, DifferentialEquations, Plots, ControlSystemsBase
using ModelingToolkitStandardLibrary, ModelingToolkitStandardLibrary.Electrical, ModelingToolkitStandardLibrary.Mechanical.Rotational, ModelingToolkitStandardLibrary.Blocks
@variables t
# R = 0.5 # [Ohm] armature resistance
# L = 4.5e-3 # [H] armature inductance
# k = 0.5 # [N.m/A] motor constant
# J = 0.02 # [kg.m²] inertia
# f = 0.01 # [N.m.s/rad] friction factor
function Motor(; name, R = 0.5, L = 4.5e-3, k = 0.5, J = 0.02, f = 0.01)
    @named ground = Ground()
    @named source = Voltage()
    @named R1 = Resistor(R = R)
    @named L1 = Inductor(L = L)
    @named emf = EMF(k = k)
    @named fixed = Fixed()
    @named load = Torque(use_support = false)
    @named inertia = Inertia(J = J)
    @named friction = Damper(d = f)

    connections = [connect(fixed.flange, emf.support, friction.flange_b)
                   connect(emf.flange, friction.flange_a, inertia.flange_a)
                   connect(inertia.flange_b, load.flange)
                   connect(source.p, R1.p)
                   connect(R1.n, L1.p)
                   connect(L1.n, emf.p)
                   connect(emf.n, source.n, ground.g)]
    subcomps = [
        ground,
        source,
        R1,
        L1,
        emf,
        fixed,
        load,
        inertia,
        friction,
    ]
    @named model = ODESystem(connections, t)
    compose(model, subcomps)
end

pi_k = 1.1
pi_T = 0.05
@named motor = Motor();
tau_L_step = -0.3 # [N.m] amplitude of the load torque step
@named ref = Blocks.Step(height = 1, start_time = 0)
@named pi_controller = Blocks.LimPI(k = pi_k, T = pi_T, u_max = 10, Ta = 0.035)
@named feedback = Blocks.Feedback()
@named load_step = Blocks.Step(height = tau_L_step, start_time = 3)
@named speed_sensor = SpeedSensor()

connections = [
            connect(motor.load.flange, speed_sensor.flange)
            connect(ref.output, feedback.input1)
            connect(speed_sensor.w, :y, feedback.input2)
            connect(load_step.output, motor.load.tau)
            connect(feedback.output, pi_controller.err_input)
            connect(pi_controller.ctr_output, :u, motor.source.V)]

subcomps = [
    motor,
    ref,
    pi_controller,
    feedback,
    load_step,
    speed_sensor,
]
@named model = ODESystem(connections, t)
model = compose(model, subcomps)

sys = structural_simplify(model)
prob = ODEProblem(sys, [], (0, 6.0))
sol = solve(prob, Rodas4())

p1 = plot(sol.t, sol[motor.inertia.w], ylabel = "Angular Vel. in rad/s",
                label = "Measurement", title = "DC Motor with Speed Controller")
Plots.plot!(sol.t, sol[ref.output.u], label = "Reference")
p2 = plot(sol.t, sol[motor.load.tau.u], ylabel = "Disturbance in Nm", label = "")

mat, simplified_sys = get_sensitivity(model, :y);
S = ss(mat...);
bplot = bodeplot(S, plotphase=false)
nplot = nyquistplot(-ss(get_looptransfer(model, :u)[1]...))
plot(p1, p2, bplot, nplot, layout = (2, 2))
```
