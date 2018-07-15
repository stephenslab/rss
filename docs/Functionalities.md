---
layout: default
---

[Zhu and Stephens (*Ann. Appl. Stat.*, 2017)]: https://projecteuclid.org/euclid.aoas/1507168840
[Zhu and Stephens (*bioRxiv*, 2017)]: https://doi.org/10.1101/160770 

## Model Fitting based on Markov chain Monte Carlo (MCMC)

- [**`rss_bvsr.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_bvsr.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Bayesian variable selection regression" (BVSR; [Guan and Stephens, 2011](https://projecteuclid.org/euclid.aoas/1318514285)) prior using a Metropolis-Hastings algorithm.
- [**`rss_bslmm.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_bslmm.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Bayesian sparse linear mixed model" (BSLMM; [Zhou, Carbonetto and Stephens, 2013](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003264)) prior using a component-wise MCMC algorithm.
- [**`rss_ash.m`**](https://github.com/stephenslab/rss/blob/master/src/rss_ash.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the "Adaptive shrinkage" (ASH; [Stephens, 2017](https://doi.org/10.1093/biostatistics/kxw041)) prior using a component-wise MCMC algorithm.

Details of MCMC algorithms are available in
[Supplementary Appendix B](http://stephenslab.uchicago.edu/assets/papers/Zhu2017-supplement.pdf)
of [][]
and [here](http://www.stat.uchicago.edu/~xiangzhu/rss_mcmc.pdf). 

## Model Fitting based on Variational Bayes (VB)

- [**`rss_varbvsr.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the BVSR prior using a mean-field VB algorithm. This can be viewed as an extension of [Carbonetto and Stephens (2012)](https://projecteuclid.org/euclid.ba/1339616726) for the analysis of summary-level data.
- [**`rss_varbvsr_squarem.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_squarem.m) <br> Fit the Bayesian model that consists of the RSS likelihood and the BVSR prior using a mean-field VB algorithm and a [SQUAREM (Varadhan and Roland, 2008)](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9469.2007.00585.x/abstract) accelerator.
- [**`rss_varbvsr_parallel.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_parallel.m) <br> Parallel implementation of [`rss_varbvsr.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr.m).
- [**`rss_varbvsr_pasquarem.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_pasquarem.m) <br> Parallel implementation of [`rss_varbvsr_squarem.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_squarem.m)
- [**`rss_varbvsr_bigmem.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem.m) <br> "Big-data" implementation of [`rss_varbvsr_parallel.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_parallel.m).
- [**`rss_varbvsr_bigmem_squarem.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem_squarem.m) <br> "Big-data" implementation of [`rss_varbvsr_pasquarem.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_pasquarem.m).

## Miscellaneous

- [**`import_1000g_vcf.sh`**](https://github.com/stephenslab/rss/blob/master/misc/import_1000g_vcf.sh) <br> Output [1000 Genomes](http://www.internationalgenome.org/data) phased haplotypes of a given list of SNPs in [IMPUTE reference-panel format](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#input_options). 
- [**`compute_pve.m`**](https://github.com/stephenslab/rss/blob/master/src/compute_pve.m) <br> Compute the estimated PVE (or SNP heritability, defined in [Guan and Stephens (2011)](https://projecteuclid.org/euclid.aoas/1318514285)) based on GWAS summary data.
- [**`band_storage.m`**](https://github.com/stephenslab/rss/blob/master/misc/band_storage.m) <br> Convert a symmetric, banded matrix to a compact matrix in such a way that only the main diagonal and the nonzero superdiagonals are stored.
- [**`find_bandwidth.m`**](https://github.com/stephenslab/rss/blob/master/misc/find_bandwidth.m) <br> Find the bandwidth of a symmetric, banded matrix.
- [**`get_corr.m`**](https://github.com/stephenslab/rss/blob/master/misc/get_corr.m) <br> Compute (LD) matrix using the shrinkage estimator proposed in [Wen and Stephens (2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081).
- [**`data_maker.m`**](https://github.com/stephenslab/rss/blob/master/misc/data_maker.m) <br> Simulate phenotype data from the genome-wide multiple-SNP model described in [Zhou, Carbonetto and Stephens (2013)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003264), and then compute the single-SNP summary statistics for each SNP.
- [**`enrich_datamaker.m`**](https://github.com/stephenslab/rss/blob/master/misc/enrich_datamaker.m) <br> Simulate phenotype data from the genetic association enrichment model described in [Carbonetto and Stephens (2013)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003770), and then compute the single-SNP summary statistics for each SNP.
- [**`null_single.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/null_single.m) & [**`null_template.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/null_template.m) <br> Run genome-wide multiple-SNP analysis of single-SNP summary data, using `rss_varbvsr*` functions.
- [**`gsea_wrapper.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/gsea_wrapper.m) & [**`gsea_template.m`**](https://github.com/stephenslab/rss/blob/master/src_vb/gsea_template.m) <br> Run genome-wide enrichment analysis of single-SNP summary data, using `rss_varbvsr*` functions.
