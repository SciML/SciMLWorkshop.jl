# Training Neural Stochastic Differential Equations with GPU acceleration (I)

In the previous models we had to define a model. Now let's shift the burden of
model-proofing onto data by utilizing neural differential equations. A neural
differential equation is a differential equation where the model equations
are replaced, either in full or in part, by a neural network. For example, a
neural ordinary differential equation is an equation $u^\prime = f(u,p,t)$
where $f$ is a neural network. We can learn this neural network from data using
various methods, the easiest of which is known as the single shooting method,
where one chooses neural network parameters, solves the equation, and checks
the ODE's solution against data as a loss.

In this example we will define and train various forms of neural differential
equations. Note that all of the differential equation types are compatible with
neural differential equations, so this is only going to scratch the surface of
the possibilites!

## Part 1: Constructing and Training a Basic Neural ODE

Use the [DiffEqFlux.jl README](https://github.com/JuliaDiffEq/DiffEqFlux.jl) to
construct a neural ODE to train against the training data:

```{julia;eval=false}
u0 = Float32[2.; 0.]
datasize = 30
tspan = (0.0f0,1.5f0)

function trueODEfunc(du,u,p,t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u.^3)'true_A)'
end
t = range(tspan[1],tspan[2],length=datasize)
prob = ODEProblem(trueODEfunc,u0,tspan)
ode_data = Array(solve(prob,Tsit5(),saveat=t))
```

## Part 2: GPU-accelerating the Neural ODE Process

Use the `gpu` function from Flux.jl to transform all of the calculations onto
the GPU and train the neural ODE using GPU-accelerated `Tsit5` with adjoints.

## Part 3: Defining and Training a Mixed Neural ODE

Gather data from the Lotka-Volterra equation:

```{julia;eval=false}
function lotka_volterra(du,u,p,t)
  x, y = u
  α, β, δ, γ = p
  du[1] = dx = α*x - β*x*y
  du[2] = dy = -δ*y + γ*x*y
end
u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0]
prob = ODEProblem(lotka_volterra,u0,tspan,p)
sol = Array(solve(prob,Tsit5())(0.0:1.0:10.0))
```

Now use the
[mixed neural section of the documentation](https://github.com/JuliaDiffEq/DiffEqFlux.jl#mixed-neural-des)
to define the mixed neural ODE where the functional form of $\frac{dx}{dt}$ is
known, and try to derive a neural formulation for $\frac{dy}{dt}$ directly from
the data.

## Part 4: Constructing a Basic Neural SDE

Generate data from the Lotka-Volterra equation with multiplicative noise

```{julia;eval=false}
function lotka_volterra(du,u,p,t)
  x, y = u
  α, β, δ, γ = p
  du[1] = dx = α*x - β*x*y
  du[2] = dy = -δ*y + γ*x*y
end
function lv_noise(du,u,p,t)
  du[1] = p[5]*u[1]
  du[2] = p[6]*u[2]
end
u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0,0.1,0.1]
prob = SDEProblem(lotka_volterra,lv_noise,u0,tspan,p)
sol = [Array(solve(prob,SOSRI())(0.0:1.0:10.0)) for i in 1:20] # 20 solution samples
```

Train a neural stochastic differential equation $dX = f(X)dt + g(X)dW_t$ where
both the drift ($f$) and the diffusion ($g$) functions are neural networks.
See if constraining $g$ can make the problem easier to fit.

## Part 5: Optimizing the training behavior with minibatching (E)

Use minibatching on the data to improve the training procedure. An example
[can be found at this PR](https://github.com/FluxML/model-zoo/pull/88).
