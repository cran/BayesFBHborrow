# Package: BayesFBHborrow
- Version: 2.0.9
- Title: Bayesian Borroing for TTE data from a Flexible Baseline Hazard Function
- Authors: Darren Scott, Sophia Axillus and Grant Izmirlian
- Imports: rlang, dplyr, invgamma, mvtnorm, checkmate, magrittr, ggplot2, patchwork,
  	   kableExtra, stats, survival, survminer, extraDistr, bayestestR

- Description: 

  + Allows Bayesian borrowing from a historical dataset for time-to-
    event data. A flexible baseline hazard function is achieved via a piecewise
    exponential likelihood with time varying split points and smoothing prior on the
    historic baseline hazards. The method is described in Scott and Lewin (2026)
    <doi:10.1093/biostatistics/kxag006>, and a paper focused on the software is
    in Scott, Axillus, Lewin and Izmirlian (2026) <doi:10.48550/arXiv.2408.04327>.

- License: Apache License (>=2)
