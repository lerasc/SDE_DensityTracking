
# SDE density tracking

The implementation consists of one single class that tracks the density distribution of the solution to a stochastic differential equation (SDE). 
Starting from a generic SDE, with some time-independent drift and volatility functions, the class calculates the probability density to be at 
a certain position at a given time. The class can furthermore deal with absorbing or reflective boundary conditions. 

Technically, we are not solving the Fokker-Planck equation, but the Chapman-Kolmogorov forward equation. 
This implementation has the additional advantage that no discretization grid must be specified, and drift and volatility functions must not necessarily be differentiable. 

See description.pdf and documentation inside SDE_DensityTracking.py for details. 

![Exaple Output](https://github.com/slera90/SDE_DensityTracking/blob/master/description/PDE_solution.png)