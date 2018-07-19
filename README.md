
This class provides a numerical solution of the Chapman-Kolmogorov equation / Fokker-Planck equation 
Consider a stochastic differential equation of the form
        dX_t 	= mu(X_t) dt + sigma(X_t) dW_t,    X_0 = x0 												   (1)
where mu and sigma are two functions of X, and x0 is some initial value. In this class, we calculate numerically
the probability p(x,t) of being at position x at time t. Instead of solving the Fokker-Planck equation associated
with equation (1), we use a method called density tracking by quadrature (DTQ) [1], which first disretizes (1) in
time as
    X_{i+1} = X_i + mu(X_t) h + sigma(X_i) sqrt(h) Z_{i+1},													(2)
with h the discretization time-step, and Z_{i+1} a Gaussian random variable with zero mean and unit variance.
Then,  we apply the discretized Chapman-Kolmogorov equation to get
    p(x, t_{i+1}) = k sum_{j=-M}^M G(x, y_j) p(y_j, t_i) 														(3)
where k is the spatial discretization width, G(x,y) is the density of being at position x at time i+1 when at
position y at time t, and M is some finite size cut-off. Details are found in [1] and a similar implementation
exists for R [2].

references:
----------
[1]	Bhat, H.S., and Madushani, R.W.M.A, "Density tracking by quadrature for stochastic differential equations."
    arXiv preprint arXiv:1610.09572 (2016).
[2]	Rdtq: Density Tracking by Quadrature (https://cran.r-project.org/web/packages/Rdtq/index.html)
[3] Veestraeten, D. "The conditional probability density function for a reflected Brownian motion."
    Computational Economics 24.2 (2004): 185-207.
[4] Molini, A., et al. "First passage time statistics of Brownian motion with purely time dependent drift and
    diffusion." Physica A: Statistical Mechanics and its Applications 390.11 (2011): 1841-1852.      
