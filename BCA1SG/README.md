# BCA1SG: Block Coordinate Ascent with One-Step GR
## Goal
`BCA1SG` is an R package implementing the block coordinate ascent with one-step GR (BCA1SG) algorithm on the semiparametric models for panel count data, interval-censored survival data, and degradation data.

## Details
* A comprehensive description of the BCA1SG algorithm can be found in Wang et al. (2020).
* For panel count data, we focus on the semiparametric nonhomogeneous Poisson process model in Wellner and Zhang (2007).
* For interval-censored survival data, we focus on the semiparametric proportional hazard model in Section 4 of Huang and Wellner (1997).
* For degradation data, we focus on the semiparametric random-effects inverse Gaussian process model in Section 3 of Wang and Xu (2010).

## References
Wang Y., Ye, Z.-S., and Cao, H.(2020). On Computation of Semi-Parametric Maximum Likelihood Estimators with Shape Constraints. Submitted.

Wellner J.A. and Zhang Y.(2007). Two Likelihood-Based Semiparametric Estimation Methods for Panel Count Data with Covariates. The Annals of Statistics, 35(5), 2106-2142.

Huang J. and Wellner, J.A.(1997). Interval-Censored Survival Data: A Review of Recent Progress. Proceedings of the Fifth Seattle Symposium in Biostatistics, 123-169.

Wang X. and Xu, D.(2010). An Inverse Gaussian Process Model for Degradation Data. Technometrics, 52(2), 188-197.
