# PyCav

## Objective

* Learn integration technique for computing required statistical moments of the evolving bubble dynamics

## Todo

- [x] Verify Monte Carlo gets moments right
- [x] Implement Simpson's rule
- [x] Turn pressure from constant to time-dependent function
* [x] SSP-RK2
* [x] SSP-RK3
* [x] Adaptive Euler/RK2
- [x] Adaptive RK23
- [x] Linear bubble dynamics 
* [x] Keller Miksis
* [x] Stage 1 complete
* [x] Add $dpdt$ to Keller Miksis

## Development plan

Stage 1:
* Uncoupled polydisperse bubble cavitation 
  * Modular bubble model (e.g. RPE, KM, etc.)
  * Modular bubble distribution in $R_o$ coordinate (default to log-normal)
* Time dependence for bubble forcing (pressure $p(t)$
* Approximate exact moments via Monte Carlo

Stage 2:
* One-way coupled dynamics
  * Pressure affects bubble state (radii, void fraction, number density function, etc.)
  * Changing bubble state does not effect time-dependent pressure forcing
  * Use high-class-count to approximate exact solution

Stage 3:
* Fully-coupled dynamics
* Spatial dependencies
