# PyCav

## Todo

- [x] Verify Monte Carlo gets moments right
- [x] Implement Simpson's rule
- [ ] Turn pressure from constant to time-dependent function
- [ ] Adaptive RK23 
- [ ] Linear bubble dynamics 
* [ ] Keller Miksis

## Goals

Stage 1:
* Uncoupled polydisperse bubble cavitation 
  * Modular bubble model (e.g. RPE, KM, etc.)
  * Modular bubble distribution in $R_o$ coordinate (default to log-normal)
* Python code suitable for ML
* Learn integration technique for computing required statistical moments of the evolving bubble dynamics
* Time dependence for bubble forcing (pressure $p(t)$
* Approximate exact moments via Monte Carlo

Stage 2:
* One-way coupled dynamics
  * Has actual pressure waves, spatial dependence, etc.
  * Pressure affects bubble state (radii, void fraction, number density function, etc.)
  * Changing bubble state does not effect time-dependent pressure forcing
  * Use high-class-count to approximate exact solution

Stage 3:
* Fully-coupled dynamics
