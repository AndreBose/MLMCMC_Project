# Multi-Level Monte Carlo Markov Chain

Many natural phenomena and engineering applications are described by *Differential Models* expressed in form of *Partial Differential Equations* (PDE) or *Ordinary Differential Equations* (ODE). In several scenarios, models depend on unknown parameters. A problem of applicative interest is to estimate them from a noisy and partial observation of the solution. The Parameter Estimation problem can be tackled within a *Bayesian framework*, setting an initial prior probability distribution for the parameters and using observations to update the knowledge about them computing the posterior. However, the resulting posterior may be very complicated and not easy to sample from, since the likelihood requires the exact solution of the differential problem, that is usually impossible to have in closed form. Numerical solvers provide
accurate approximations but, when the complexity of the model and the required accuracy increase, the computation can become extremely time consuming, making standard
*Monte Carlo Markov Chain* (MCMC) procedure unfeasible.

In this project we explore the *Multi Level Monte Carlo Markov Chain* (MLMCMC), that applies the standard MCMC procedure solving the differential problem at different levels of accuracy, so that the *coarse levels* (faster but less accurate) propose samples to the *fine levels* (slower but more accurate), hoping to achieve more efficient and less correlated samplings and a reduction of the variance when estimating expected values.

## Code

1. ```Comet_equation```: we apply MLMCMC to an Advection-Diffusion PDE adopting nested levels of accuracy of the numerical solver and tuning the likelihood variance (see ```Comet_equation_M_MLDA_QoIs.ipynb```). We also introduce a Neural Networks-based Surrogate Model at coarse level and studythe impact of prior settings, data collection scheme and hyper-parameters tuning (see ```Comet_equation_M_MLDA_QoIs_red_mod.ipynb```).
2. ```Dynamical_systems```: we apply MLMCMC to a Compartmental ODE model for Epidemiology (SIR and SEIR) adopting nested levels of refinements; we tune the multilevel parameters and study the impact of prior settings.

## Report

For further details see ```/Deliverables/4 (15 02) Final Report```.

