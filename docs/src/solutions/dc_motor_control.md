# Controlling a DC Motor (E)

```@example dc_motor
using ModelingToolkit
using ModelingToolkit: @mtkmodel, @mtkcompile
const t = ModelingToolkit.t_nounits
const D = ModelingToolkit.D_nounits
using DifferentialEquations, Plots, ControlSystemsBase
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks

# DC Motor model using @mtkmodel
@mtkmodel Motor begin
    @structural_parameters begin
        R = 0.5      # [Ohm] armature resistance
        L = 4.5e-3   # [H] armature inductance
        k = 0.5      # [N.m/A] motor constant
        J = 0.02     # [kg.m²] inertia
        f = 0.01     # [N.m.s/rad] friction factor
    end
    @components begin
        ground = Ground()
        source = Voltage()
        R1 = Resistor(R = R)
        L1 = Inductor(L = L)
        emf = EMF(k = k)
        fixed = Fixed()
        load = Torque(use_support = false)
        inertia = Inertia(J = J)
        friction = Damper(d = f)
    end
    @equations begin
        connect(fixed.flange, emf.support, friction.flange_b)
        connect(emf.flange, friction.flange_a, inertia.flange_a)
        connect(inertia.flange_b, load.flange)
        connect(source.p, R1.p)
        connect(R1.n, L1.p)
        connect(L1.n, emf.p)
        connect(emf.n, source.n, ground.g)
    end
end

# PI Controller parameters
pi_k = 1.1
pi_T = 0.05
tau_L_step = -0.3  # [N.m] amplitude of the load torque step

# Full system model with motor and controller
@mtkmodel DCMotorControlSystem begin
    @components begin
        motor = Motor()
        ref = Blocks.Step(height = 1, start_time = 0)
        pi_controller = Blocks.LimPI(k = pi_k, T = pi_T, u_max = 10, Ta = 0.035)
        feedback = Blocks.Feedback()
        load_step = Blocks.Step(height = tau_L_step, start_time = 3)
        speed_sensor = SpeedSensor()
    end
    @equations begin
        connect(motor.load.flange, speed_sensor.flange)
        connect(ref.output, feedback.input1)
        connect(speed_sensor.w, :y, feedback.input2)
        connect(load_step.output, motor.load.tau)
        connect(feedback.output, pi_controller.err_input)
        connect(pi_controller.ctr_output, :u, motor.source.V)
    end
end

@mtkcompile sys = DCMotorControlSystem()

prob = ODEProblem(sys, [], (0, 6.0))
sol = solve(prob, Rodas4())

p1 = plot(sol.t, sol[sys.motor.inertia.w], ylabel = "Angular Vel. in rad/s",
                label = "Measurement", title = "DC Motor with Speed Controller")
Plots.plot!(sol.t, sol[sys.ref.output.u], label = "Reference")
p2 = plot(sol.t, sol[sys.motor.load.tau.u], ylabel = "Disturbance in Nm", label = "")

plot(p1, p2, layout = (2, 1))
```
