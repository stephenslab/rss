---
layout: default
---

[Wen and Stephens (*Ann. Appl. Stat.*, 2010)]: https://www.ncbi.nlm.nih.gov/pubmed/21479081
[Li and Stephens (*Genetics*, 2003)]: https://www.ncbi.nlm.nih.gov/pubmed/14704198

## General information

> Q: I can't afford `MATLAB`. It would be nice if RSS methods could be made to work with some open source language so everyone can run it.

- A: To further improve the accessibility of RSS methods,
I am working on an `R` package `rssr` for RSS methods:
[https://github.com/stephenslab/rssr](https://github.com/stephenslab/rssr).

## RSS methods based on Markov chain Monte Carlo (MCMC)

## RSS methods based on Variational Bayes (VB)

## Estimation of linkage disequilibrium (LD)

> Q: The shrinkage estimator of LD in [Wen and Stephens (*Ann. Appl. Stat.*, 2010)][]
requires the "scaled population recombination rate" (`rho_ij`). How do I compute it?

- A: One can easily figure out how to calculate this quantity based on
[Li and Stephens (*Genetics*, 2003)][] and [Wen and Stephens (*Ann. Appl. Stat.*, 2010)][].
In case you do not have time to read these two great papers,
I create a short [tutorial](Recombination) to illustrate the calculation of
scaled population recombination rate using HapMap genetic map.

> Q: Can I use unphased genotype data to compute the shrinkage estimator?

- A: Yes. You need to modify [Line 46](https://github.com/stephenslab/rss/blob/master/misc/get_corr.m#L46)
of [`get_corr.m`](https://github.com/stephenslab/rss/blob/master/misc/get_corr.m):
`S = 0.5*cov(Gpanel)`, where `Gpanel` is the genotype matrix from an external reference panel.
Please see Section 2.4 of [Wen and Stephens (*Ann. Appl. Stat.*, 2010)][] for more details.

## Applications

> Q: If I want to use RSS, how should I prepare for the input data?

- A: Thanks for your interest in trying RSS methods!
To help you prepare your own data for RSS-based analyses,
I create a short [tutorial](Input-Data-Formats) about input data formats. 
