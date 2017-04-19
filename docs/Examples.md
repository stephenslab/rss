We illustrate how to use RSS software through the following examples.

- [**Example 1**](Example-1): Fit the RSS models via MCMC, and
  estimate the SNP heritability (PVE).

> This example illustrates how to fit RSS models using MCMC
> algorithms. Three types of prior distributions are considered: BVSR,
> BSLMM and ASH. The output of MCMC is further used to estimate the
> SNP heritability.

- [**Example 2**](Example-2): Fit RSS-BVSR via MCMC with three types
  of LD matrices.

> This example illustrates the impact of different LD estimates on the
> RSS results. Three types of estimated LD matrices are considered:
> cohort sample LD, shrinkage panel sample LD [(Wen and Stephens, 2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081)
> and panel sample LD.

- [**Example 3**](Example-3): Fit RSS-BVSR via MCMC with two types of
  SE vectors.

> This example illustrates the impact of two definitions of `se` on
> the RSS results.

- [**Example 4**](Example-4): Fit RSS-BVSR via variational Bayes (VB)
  methods.

> This example illustrates how to fit an RSS-BVSR model using
> variational Bayes (VB) approximation, compare the results with
> previous work based on individual-level data [(Carbonetto and Stephens, 2012)](https://projecteuclid.org/euclid.ba/1339616726).
