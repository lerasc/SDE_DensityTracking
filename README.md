
This class provides a numerical solution of the Chapman-Kolmogorov equation / Fokker-Planck equation 
Consider a stochastic differential equation of the form
$$
\mathrm{d}X_t   = \mu(X_t) \mathrm{d}t + \sigma(X_t) \mathrm{d}W_t, ~~   X_0 = x0   
$$       
where $\mu$ and $\sigma$ are two functions of $X$, and $x0$ is some initial value. 
In this class, we calculate numerically the probability $p(x,t)$ of being at position $x$ at time $t$. 
Instead of solving the Fokker-Planck equation associated
with the above equation we use a method called density tracking by quadrature (DTQ) [1], which first disretizes 
the equation in time as
$$
X_{i+1} = X_i + mu(X_t) h + sigma(X_i) sqrt(h) Z_{i+1},
$$
with $h$ the discretization time-step, and $Z_{i+1}$ a Gaussian random variable with zero mean and unit variance.
We can see immediately that the position of $X_{i+1}$ is distributed according to a Gaussian 
with mean $X_i +  \mu(X_i) h$ and standard deviation $\sigma(X_i) \sqrt{ h}$. 
Let us denote by $G(x,y)$ a Gaussian with mean $y+\mu(y) h$ and standard deviation $\sigma(y) \sqrt{h}$. 
Then, the density at position $x$ at time $t+1$ is just given by the Chapman-Kolmogorov equation, which, 
in discretized form, reads
$$
p(x, t_{i+1}) = k sum_{j=-M}^M G(x, y_j) p(y_j, t_i). 
$$
Here, $k$ is the spatial discretization width and $M$ is some finite size cut-off. 
Details are found in the documentation of SDE_DensityTracking.py, and in [1]. 
A similar implementation exists for R [2].

references:
----------
[1] Bhat, H.S., and Madushani, R.W.M.A, "Density tracking by quadrature for stochastic differential equations."
    arXiv preprint arXiv:1610.09572 (2016). <br>
[2] Rdtq: Density Tracking by Quadrature (https://cran.r-project.org/web/packages/Rdtq/index.html) <br>
[3] Veestraeten, D. "The conditional probability density function for a reflected Brownian motion."
    Computational Economics 24.2 (2004): 185-207. <br>
[4] Molini, A., et al. "First passage time statistics of Brownian motion with purely time dependent drift and
    diffusion." Physica A: Statistical Mechanics and its Applications 390.11 (2011): 1841-1852.      
