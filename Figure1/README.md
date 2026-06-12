# Figure 1

This folder contains the MATLAB code used to generate Figure 1c–e of the manuscript.

## Description

Figure 1 compares numerical simulations of the Gray–Scott reaction–diffusion model with experimentally observed metal-thiolate coordination polymer networks (MTCPs).

Panels c–e correspond to simulated spatial patterns obtained from the Gray–Scott model under different inhibitor diffusion coefficients.

## Model

Gray–Scott reaction–diffusion system on a two-dimensional periodic domain.

## Parameters

Fixed parameters:

* F = 0.047
* k = 0.060
* Du = 1.20

Variable parameter:

* Dv = 0.60 (Figure 1c)
* Dv = 0.36 (Figure 1d)
* Dv = 0.22 (Figure 1e)

## Numerical Method

* Spatial discretization: finite-difference Laplacian
* Boundary condition: periodic
* Time integration: fourth-order Runge–Kutta (RK4)

## Output

The generated patterns correspond to Figure 1c–e in the manuscript.

## Related Files

* Figure1_simulation.m
* Figure1_plot.m
