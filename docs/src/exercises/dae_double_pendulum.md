# Differential-Algebraic Equation Modeling of a Double Pendulum (B)

Differential-Algebraic Equation (DAE) systems are like ODEs but allow for adding
constraints into the models. This problem will look at solving the double
pendulum problem with enforcement of the rigid body constraints, requiring that
the total distance `L` is constant throughout the simulation. While these
equations can be rewritten in an ODE form, in many cases it can be simpler
to solve the equation directly with the constraints. This tutorial will
cover both the idea of index, how to manually perform index reduction,
and how to make use of mass matrix and implicit ODE solvers to handle these
problems.

## Part 1: Simple Introduction to DAEs: Mass-Matrix Robertson Equations

A mass-matrix ordinary differential equation (ODE) is an ODE where the
left-hand side, the derivative side, is multiplied by a matrix known as the
mass matrix. This is described as:

```math
Mu' = f(u,p,t)
```

where $M$ is the mass matrix. When $M$ is invertible, there is an ODE which is
equivalent to this formulation. When $M$ is not invertible, this can have a
distinctly different behavior and is as Differential-Algebraic Equation (DAE).

Solve the Robertson DAE:

```math
\begin{align*}
\frac{dy_1}{dt} &= -0.04y_1 + 10^4 y_2y_3\\
\frac{dy_2}{dt} &=  0.04y_1 - 10^4 y_2y_3 - 3\times 10^7 y_2^2\\
1 &= y_1 + y_2 + y_3
\end{align*}
```

with $y(0) = [1,0,0]$ and $dy(0) = [-0.04,0.04,0.0]$ using the mass-matrix
formulation and `Rodas5()`. Use the
[ODEProblem page](https://docs.sciml.ai/dev/types/ode_types)
to find out how to declare a mass matrix.

(Hint: what if the last row has all zeros?)

## Part 2: Solving the Implicit Robertson Equations with IDA

Use the [DAE Tutorial](https://docs.sciml.ai/dev/tutorials/dae_example)
to define a DAE in its implicit form and solve the Robertson equation with IDA.
Why is `differential_vars = [true,true,false]`?

## Part 3: Manual Index Reduction of the Single Pendulum

The index of a DAE is a notion used to measure distance from
its related ODE. There are many different definitions of index,
but we're going to stick to the idea of differential index:
the number of differentiations required to convert a system
of DAEs into explicit ODE form. DAEs of high index are
usually transformed via a procedure called index reduction.
The following example will demonstrate this.

Consider the index 3 DAE system of the cartesian pendulum.
After writing down the force equations in both directions,
we arrive at the following DAE:

```math
\begin{align*}
m\ddot{x} &= \frac{x}{L}T \\
m\ddot{y} &= \frac{y}{L}T - mg \\
x^2 + y^2 &= L
\end{align*}
```

Notice that we don't have an equation describing the
behaviour of `T`. Let us now perform index reduction to
extract an equation for `T`

Differentiate this third equation twice with respect to time
to reduce it from index 3 to index 1.

## Part 4: Single Pendulum Solution with IDA

Write these equations in implicit form and solve the system using
IDA.

## Part 5: Solving the Double Pendulum DAE System

The following equations describe a double
pendulum system:

```math
\begin{align*}
m_2\ddot{x_2} &= \frac{x_2}{L_2}T_2 \\
m_2\ddot{y_2} &= \frac{y_2}{L_2}T_2 - m_2g \\
{x_2}^2 + {y_2}^2 &= L_2 \\
m_1\ddot{x_1} &= \frac{x_1}{L_1}T_1 - \frac{x_2}{L_2}T_2 \\
m_2\ddot{y_1} &= \frac{y_1}{L_1}T_2 - m_1g - \frac{y_2}{L_2}T_2 \\
{x_1}^2 + {y_1}^2 &= L_1 \\
\end{align*}
```

Perform index reduction and solve it like in the previous example.
