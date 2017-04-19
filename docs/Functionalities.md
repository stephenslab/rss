---
layout: default
---

### Model Fitting based on Markov chain Monte Carlo (MCMC)

- [**`rss_bvsr.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_bvsr.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Bayesian variable selection regression" (BVSR; [Guan and Stephens, 2011](https://projecteuclid.org/euclid.aoas/1318514285)) prior.
- [**`rss_bslmm.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_bslmm.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Bayesian sparse linear mixed model" (BSLMM; [Zhou, Carbonetto and Stephens, 2013](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003264)) prior.
- [**`rss_ash.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_ash.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Adaptive shrinkage" (ASH; [Stephens, 2017](https://doi.org/10.1093/biostatistics/kxw041)) prior.

Details of MCMC algorithms are available in this [document](http://www.stat.uchicago.edu/~xiangzhu/rss_mcmc.pdf). 

### Model Fitting based on Variational Bayes (VB)

- [**`rss_varbvsr.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the BVSR prior. This can be viewed as an extension of [Carbonetto and Stephens (2012)](https://projecteuclid.org/euclid.ba/1339616726) for the analysis of summary-level data.

### Miscellaneous

- [**`import_1000g_vcf.sh`**](https://github.com/stephenslab/rss/blob/master/misc/import_1000g_vcf.sh) <br> Output [1000 Genomes](http://www.internationalgenome.org/data) phased haplotypes of a given list of SNPs in [IMPUTE reference-panel format](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#input_options). 
- [**`compute_pve.m`**](https://github.com/stephenslab/rss/blob/master/src/compute_pve.m) <br> Compute the estimated PVE (or SNP heritability, defined in [Guan and Stephens (2011)](https://projecteuclid.org/euclid.aoas/1318514285)) based on GWAS summary data.
- [**`band_storage.m`**](https://github.com/stephenslab/rss/blob/master/misc/band_storage.m) <br> Convert a symmetric, banded matrix to a compact matrix in such a way that only the main diagonal and the nonzero superdiagonals are stored.
- [**`find_bandwidth.m`**](https://github.com/stephenslab/rss/blob/master/misc/find_bandwidth.m) <br> Find the bandwidth of a symmetric, banded matrix.
- [**`get_corr.m`**](https://github.com/stephenslab/rss/blob/master/misc/get_corr.m) <br> Compute (LD) matrix using the shrinkage estimator proposed in [Wen and Stephens (2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081).
- [**`enrich_datamaker.m`**](https://github.com/stephenslab/rss/blob/master/misc/enrich_datamaker.m) <br> Simulate phenotype data from the genetic association enrichment model described in [Carbonetto and Stephens (2013)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003770), and then compute the single-SNP summary statistics for each SNP. 
