### Model Fitting based on Markov chain Monte Carlo

- [**`rss_bvsr.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_bvsr.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Bayesian variable selection regression" (BVSR; [Guan and Stephens, 2011](https://projecteuclid.org/euclid.aoas/1318514285)) prior.
- [**`rss_bslmm.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_bslmm.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Bayesian sparse linear mixed model" (BSLMM; [Zhou, Carbonetto and Stephens, 2013](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003264)) prior.
- [**`rss_ash.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_ash.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Adaptive shrinkage" (ASH; [Stephens, 2016](http://biorxiv.org/content/early/2016/01/29/038216)) prior.

The details of the model fitting algorithms are available [here](http://www.stat.uchicago.edu/~xiangzhu/rss_mcmc.pdf). 

### Model Fitting based on Variational Bayes

- [**`rss_varbvsr.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the BVSR prior. This can be viewed as an extension of [Carbonetto and Stephens (2012)](https://projecteuclid.org/euclid.ba/1339616726) for the analysis of summary-level data.

### Miscellaneous 

- [**`compute_pve.m`**](https://github.com/stephenslab/rss/blob/master/src/compute_pve.m) <br> compute the estimated PVE (or SNP heritability, defined in [Guan and Stephens (2011)](https://projecteuclid.org/euclid.aoas/1318514285)) based on GWAS summary data.

- [**`band_storage.m`**](https://github.com/stephenslab/rss/blob/master/misc/band_storage.m) <br> convert a symmetric, banded matrix to a compact matrix in such a way that only the main diagonal and the nonzero superdiagonals are stored.

