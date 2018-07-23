
# SDE density tracking

The implementation consists of one single class that tracks the density distribution of the solution to a stochastic differential equation (SDE). 
Starting from a generic SDE
$$
dX_t = \mu(X_t)~ dt + \sigma(X_t) ~dW_t
$$
with some time-independent drift and volatility functions $\mu$ and $\sigma$, and initial position $X_0$, 
the class calculates the probability density $p(x,t)$, to be at position $x$ at time $t$. 
The class can furthermore deal with absorbing or reflective boundary conditions. 

Technically, we are not solving the Fokker-Planck equation, but the Chapman-Kolmogorov forward equation. 
This implementation has the additional advantage that the range of $x$-values is dynamic (no previous grid must be specified), and $\mu$ and $\sigma$ must not necessarily be differentiable. 

See description.pdf and documentation inside SDE_DensityTracking.py for details. 