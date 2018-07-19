---
layout: default
---

[Zhu and Stephens (*Ann. Appl. Stat.*, 2017)]: https://projecteuclid.org/euclid.aoas/1507168840
[Zhu and Stephens (*bioRxiv*, 2017)]: https://doi.org/10.1101/160770 

We illustrate how to use RSS methods and software through the following examples.

## [Example 1](Example-1): Fit RSS models via MCMC, and estimate SNP heritability (PVE).

This example illustrates how to fit RSS models using MCMC
algorithms. Three types of prior distributions are considered:

- BVSR in [Guan and Stephens (*Ann. Appl. Stat.*, 2011)](https://projecteuclid.org/euclid.aoas/1318514285),
- BSLMM in [Zhou, Carbonetto and Stephens (*PLoS Genet.*, 2013)](https://doi.org/10.1371/journal.pgen.1003264),
- ASH in [Stephens (*Biostatistics*, 2017)](https://doi.org/10.1093/biostatistics/kxw041).

The MCMC output is further used to estimate the SNP heritability.
This example is closely related to Section 4.2 of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][]. 

## [Example 2](Example-2): Fit RSS-BVSR via MCMC with three types of LD matrices.

This example illustrates the impact of different LD estimates on the
RSS results. Three types of estimated LD matrices are considered:
cohort sample LD, panel sample LD and shrinkage panel sample LD in
[Wen and Stephens, (*Ann. Appl. Stat.*, 2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081)
This example is closely related to Section 4.1 of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][]. 

## [Example 3](Example-3): Fit RSS-BVSR via MCMC with two types of standard error (SE) vectors.

This example illustrates the impact of two definitions of
SE vector (`se`) on the RSS results.
This example is closely related to Section 2.1 of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][].  

## [Example 4](Example-4): Fit RSS-BVSR via variational Bayes (VB) methods.

This example illustrates how to fit an RSS-BVSR model using
variational Bayes (VB) approximation, and compares the results with
previous work based on individual-level data
[(Carbonetto and Stephens, *Bayesian Anal.*, 2012)](https://projecteuclid.org/euclid.ba/1339616726).
This example is closely related to Section "Connection with enrichment analysis of individual-level data" 
of [Zhu and Stephens (*bioRxiv*, 2017)][].

## [Example 5](Example-5): Enrichment analysis of GWAS summary statistics using RSS.

This example illustrates how to perform enrichment analysis of GWAS summary statistics
based on variational Bayes (VB) inference of RSS-BVSR model.
This example is closely related to the actual data analyses in [Zhu and Stephens (*bioRxiv*, 2017)][].
