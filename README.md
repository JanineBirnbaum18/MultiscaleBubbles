# Multiscale Vesiculation, Fluid flow, Failure, and Interaction Nonlinear model (MVFFIN)
Simulate bubble growth and multi-phase flow

## Technologies
Matlab software using the bubble-scale vesiculation model of Coumans et al. (2020) [https://doi.org/10.1016/j.jvolgeores.2020.107002]

## Contains
Coumans\_coupled.m: main function definition for suspension-scale dynamics and passing arguments to bubble-scale model
Julia\_clasts.m: sample input script for vesiculation of small grains
Numerical\_Model\_v2.m: bubble-scale model solver from Coumans et al. (2020)
Sample\_inputs.m: additional sample input script
create\_axes.m: support plotting function
getFunctions\_dynamic\_dimensionless.m: support functions with numerical solver for linearized Navier-Stokes flow at the suspension scale
getFunctions\_outgas.m: support functions with constitutive relationships and solvers for diffusive outgassing
getFunctions\_permeable.m: support functions with constitutive relationships and solvers for permeable gas flow
getFunctions\_thermal.m: support fuctions with constitutive relationships and solvers for thermal diffusion
get_Functions\_v2.m: support functions and constitutive relationships for bubble-scale model from Coumans et al. (2020)
