
# SDE\_DensityTracking:  A Python Package

The implementation consists of one single class that tracks the density distribution of the solution to a stochastic differential equation (SDE). Starting from a generic SDE
$$
d X_t = \mu(X_t)~  d t + \sigma(X_t)  ~d W_t
$$
with some time-independent drift and volatility functions $\mu$ and $\sigma$, and initial position $X_0$, the class calculates the probability density $p(x,t)$, to be at position $x$ at time $t$. The class can furthermore deal with absorbing or reflective boundary conditions.

#### Input and Output

Upon initialization of the class with appropriate arguments, one starts tracking the evolution of the density by calling the *track\_density* method. Once finished, the density $p = p(x,t)$ can be retrieved by calling the *get\_density* method. The return value can be either a functional or discretized object, and it can be either a function of $x$, of $t$, or both. See documentation inside *SDE\_DensityTracking.py* for a detailed description of all input and output values.

#### Theoretical Background

In our implementation, we have followed the approach of [1], which is straight forward. We first discretize the above SDE in time, giving

$$
X_{t+1} = X_t +  \mu(X_t) \Delta t + \sigma(X_t) \sqrt{ \Delta t} Z_t
$$
with $\Delta t$ the discretization step width and $Z_t$ a normally distributed random variable. We can see immediately that the position of $X_{t+1}$ is just a Gaussian with mean $X_t +  \mu(X_t) \Delta t$ and standard deviation $\sigma(X_t) \sqrt{ \Delta t}$.

Let us denote by $G(x,y)$ a Gaussian with mean $y+\mu(y) \Delta t$ and standard deviation $\sigma(y) \sqrt{ \Delta t}$. Then, the density at position $x$ at time $t+1$ is just given by the Chapman-Kolmogorov equation
$$
p(x, t+1) = \int_{\mathbb{R}} d y ~ G(x, y) ~ p(y, t).
$$
Discretizing the integral yields
$$
p(x, t+1) \approx k  \sum_{i=-\infty}^\infty  \dd y ~ G(x, y) ~ p(y_i, t).
$$
where $k$ is the spatial discretization step, ie. $k = y_{i+1}-y_i ~ \forall i$.

For the numerical implementation, the infinite sum in is truncated naturally by just summing over all values $y_i$ calculated at time $t$. The considered value range $x_i$ at time $t+1$ is determined dynamically at runtime as the minimum and maximum value below and above which the density is less than some pre-specified precision value \texttt{tol}, that we can safely approximate as zero. (In principle, this implementation can cause some problem in case of a a bimodal initial distribution with two unconnected peaks, i.e. with regions of zero probability mass between the two peaks. Because this is a very specific case, we ignore it for here for simplicity.)

This new range of values at time $t+1$ naturally truncates the summation in \eqref{eq:CK_discrete} in the next iteration, and so on. If there is an absorbing boundary, we must just replace the Gaussian $G(x,y)$ by 
$$
p(s, t ~ ; ~s_0, \underline{s}=0) = \frac{1}{\sqrt{2 \pi \sigma^2 t}} = \left[  \exp \left( - \frac{ ( s_0 + \mu t - s )^2 }{2 \sigma^2 t} \right)  - \exp \left( - \frac{2 s_0 \mu}{\sigma^2} \right) \exp \left(- \frac{  ( s_0 - \mu t + s )^2 }{2 \sigma^2 t} \right) \right],
$$
which is the density of a random walk with drift $\mu$ and volatility $\sigma$, in presence of an absorbing boundary at $s = 0$. A formula similar to equation \eqref{eq:p_absorbing_explicit} exists for the reflective barrier. See for  instance [2] for a derivation of these results. 
Note that technically, we are not solving the Fokker-Planck equation, but the Chapman-Kolmogorov forward equation. This implementation has the additional advantage that the range of $x$-values is dynamic (no previous grid must be specified), and $\mu$ and $\sigma$ must not necessarily be differentiable.

#### Consistency Check

We test the implementation of our PDE solver by considering the case of constant drift and volatility, $\mu(x) =  -0.5, ~\sigma = 1$, initial condition $p(x, t=0) = \delta(x-10)$ and an absorbing barrier at $x=5$.  As can be seen in the figure plot, the black dashed lines on top of the colored lines indicates the analytical solution above, in excellent agreement with the numerical one.


#### References
[1] Bhat, H. S., & Madushani, R. W. M. A. (2016). Density tracking by quadrature for stochastic differential equations. arXiv preprint arXiv:1610.09572.

[2] Molini, A., et al. "First passage time statistics of Brownian motion with purely time dependent drift and
        diffusion." Physica A: Statistical Mechanics and its Applications 390.11 (2011): 1841-1852.


![Exaple Output](https://github.com/lerasc/SDE_DensityTracking/blob/master/PDE_solution.png)