# pubick
Numerical Schemes for Improved Accuracy in Convection-Dominated Transport Problems

This repository provides MATLAB implementations for the numerical solution of convection-dominated fluid-flow and transport problems. It includes the families of high-resolution schemes based on Bézier curves (PUBICK and CUBICK) and on Hermite interpolation (HPUS, FDHPUS, EDHPUS), together with several classical schemes from the literature (FOU, SMART, QUICK, ADBQUICKEST, Superbee, Minmod). A set of benchmark problems is also provided to facilitate the evaluation, comparison, and validation of these schemes.

**Schemes implemented:**
  - PUBICK
  - CUBICK
  - TOPUS
  - FOU
  - SMART
  - QUICK
  - ADBQUICKEST
  - HPUS
  - EDHPUS
  - FDHPUS
  - Superbee
  - Minmod


**Problems:**
  - Smooth linear advection (`01_linear_advection_smooth`)  
  - Linear advection with a boundary layer (`02_linear_advection_boundary_layer`)  
  - Total variation (TV) analysis for linear advection (`03_linear_advection_TV_analysis`)  
  - Viscous Burgers equation (`04_burgers_viscous`)  
  - Inviscid Burgers equation (`05_burgers_inviscid`)  
  - Buckley‑Leverett equation for two-phase flow (`06_buckley_leverett`)  
  - Euler equations for compressible flow (`07_euler_equations`)  


**References:**  
 - Oviedo, P. R.; Ribeiro, A. M.; dos Santos, R. W.; de Queiroz, R. A. B.; Rocha, B. M. (2025) *Piecewise Bézier Curve-Based Numerical Schemes for Improved Accuracy in Convection-Dominated Transport Problems.* Submitted manuscript.
 - Rojas Oviedo, P. C.; Magalhães, T. M.; Rocha, B. M.; Queiroz, R. A. B. (2024). *A bounded scheme based on Bézier curves for convection-dominated transport problems.* Computers & Mathematics with Applications. https://doi.org/10.1016/j.cam.2023.115502
