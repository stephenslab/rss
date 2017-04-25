---
layout: default
---

### General information

> Q: I can't afford `MATLAB`. It would be nice if RSS methods could be made to work with some open source language so everyone can run it.

- A: At the beginning of this project, I merely used the set of `MATLAB` codes here to test some research ideas about statistical modelling and Bayesian computation, and I did not expect anyone else in the world would use these codes. I am extremely happy to see that I was wrong! Hence, to further improve the accessibility of RSS methods, I am currently leading an `R` package development for RSS methods: https://github.com/stephenslab/rssr.

### RSS methods based on Markov chain Monte Carlo (MCMC)

### RSS methods based on Variational Bayes (VB)

### Estimate linkage disequilibrium (LD)

> Q: The shrinkage estimator of LD in [Wen and Stephens (2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081) requires the "scaled population recombination rate" (`rho_ij`). How do I compute it?

- A: One can easily figure out how to calculate this quantity based on [Li and Stephens (2003)](https://www.ncbi.nlm.nih.gov/pubmed/14704198) and [Wen and Stephens (2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081). In case you do not have time to read these two papers, I create a short [tutorial](Recombination) to illustrate the calculation of scaled population recombination rate using HapMap genetic map. 

### Applications 
