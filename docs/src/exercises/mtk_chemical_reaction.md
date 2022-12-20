# ModelingToolkit Simple Chemical Reaction ODE (B)

Implement the following system in ModelingToolkit and solve it with the initial condition:

```math
\begin{align*}
y1 &= 1\\
y2,\;...,\; y7 &= 0 \\
y8 &= 0.0057,
\end{align*}
```

And parameter values:

```math
\begin{align*}
k1 &= 1.71\\
k2 &= 280\\
k3 &= 8.32\\
k4 &= 0.69\\
k5 &= 0.43\\
k6 &= 1.81,
\end{align*}
```

solve on the interval $t\in [0, 321.8122]$, and plot the solution.

```math
\begin{align*}
\frac{\mathrm{d} \mathrm{y1}\left( t \right)}{\mathrm{d}t} =& 0.0007 + k5 \mathrm{y2}\left( t \right) + k3 \mathrm{y3}\left( t \right) - k1 \mathrm{y1}\left( t \right) \\
\frac{\mathrm{d} \mathrm{y2}\left( t \right)}{\mathrm{d}t} =&  - 8.75 \mathrm{y2}\left( t \right) + k1 \mathrm{y1}\left( t \right) \\
\frac{\mathrm{d} \mathrm{y3}\left( t \right)}{\mathrm{d}t} =&  - 10.03 \mathrm{y3}\left( t \right) + k5 \mathrm{y4}\left( t \right) + 0.035 \mathrm{y5}\left( t \right) \\
\frac{\mathrm{d} \mathrm{y4}\left( t \right)}{\mathrm{d}t} =& k1 \mathrm{y3}\left( t \right) + k3 \mathrm{y2}\left( t \right) - 1.12 \mathrm{y4}\left( t \right) \\
\frac{\mathrm{d} \mathrm{y5}\left( t \right)}{\mathrm{d}t} =&  - 1.745 \mathrm{y5}\left( t \right) + k5 \mathrm{y6}\left( t \right) + k5 \mathrm{y7}\left( t \right) \\
\frac{\mathrm{d} \mathrm{y6}\left( t \right)}{\mathrm{d}t} =& k1 \mathrm{y5}\left( t \right) + k4 \mathrm{y4}\left( t \right) + k4 \mathrm{y7}\left( t \right) - k5 \mathrm{y6}\left( t \right) - k2 \mathrm{y6}\left( t \right) \mathrm{y8}\left( t \right) \\
\frac{\mathrm{d} \mathrm{y7}\left( t \right)}{\mathrm{d}t} =&  - k6 \mathrm{y7}\left( t \right) + k2 \mathrm{y6}\left( t \right) \mathrm{y8}\left( t \right) \\
\frac{\mathrm{d} \mathrm{y8}\left( t \right)}{\mathrm{d}t} =& k6 \mathrm{y7}\left( t \right) - k2 \mathrm{y6}\left( t \right) \mathrm{y8}\left( t \right)
\end{align*}
```

## Example Solution Plot

An example of the what the solution should look like is shown below:

![](https://user-images.githubusercontent.com/1814174/206983500-fef38c1e-b35d-4457-a08b-b2d64d3a59eb.png)
