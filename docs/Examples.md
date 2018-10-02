---
layout: default
---

[Zhu and Stephens (2017)]: https://projecteuclid.org/euclid.aoas/1507168840
[Zhu and Stephens (2018)]: https://doi.org/10.1101/160770 

We illustrate how to use RSS methods and software through the following examples.

## [Example 1](Example-1): Fit RSS models via MCMC, and estimate SNP heritability (PVE).

This example illustrates how to fit RSS models using MCMC
algorithms. Three types of prior distributions are considered:
BVSR [(Guan and Stephens, 2011)](https://projecteuclid.org/euclid.aoas/1318514285),
BSLMM [(Zhou et al, 2013)](https://doi.org/10.1371/journal.pgen.1003264),
and ASH [(Stephens, 2017)](https://doi.org/10.1093/biostatistics/kxw041).
The MCMC output is further used to estimate the SNP heritability.
This example is closely related to Section 4.2 of [Zhu and Stephens (2017)][]. 

## [Example 2](Example-2): Fit RSS-BVSR via MCMC with three types of LD matrix.

This example illustrates the impact of different LD estimates on the
RSS results. Three types of estimated LD matrices are considered:
sample LD based on cohort individuals, sample LD based on panel individuals,
and shrinkage LD estimate based on panel individuals
[(Wen and Stephens, 2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081).
This example is closely related to Section 4.1 of [Zhu and Stephens (2017)][]. 

## [Example 3](Example-3): Fit RSS-BVSR via MCMC with two types of standard error (SE).

This example illustrates the impact of two definitions of SE on RSS results.
This example is closely related to Section 2.1 of [Zhu and Stephens (2017)][].  

## [Example 4](Example-4): Fit RSS-BVSR via variational Bayes (VB) methods.

This example illustrates how to fit an RSS-BVSR model using
variational Bayes (VB) approximation, and compares the results with
previous work based on individual-level data
[(Carbonetto and Stephens, 2012)](https://projecteuclid.org/euclid.ba/1339616726).
This example is closely related to the section entitled
"Connection with enrichment analysis of individual-level data" of [Zhu and Stephens (2018)][].

## [Example 5](Example-5): Enrichment and prioritization analysis based on RSS.

This example illustrates how to perform enrichment and prioritization analysis
of GWAS summary statistics based on VB inference of RSS-BVSR model.
This example consists of two parts:

### [Example 5, Part A](Example-5A)

This part shows enrichment and prioritization analysis of a synthetic dataset
used in simulation studies of [Zhu and Stephens (2018)][].
Part A gives users a quick view of how RSS works in enrichment and prioritization analysis.

### [Example 5, Part B](Example-5B)

This part shows an end-to-end enrichment and prioritization analysis of
inflammatory bowel disease GWAS summary statistics
[(Liu et al, 2015)](https://www.ncbi.nlm.nih.gov/pubmed/26192919)
and a gene set named *IL23-mediated signaling events* (Pathway Commons 2, PID, 37 genes).
Part B illustrates the actual data analyses performed in [Zhu and Stephens (2018)][].
