---
layout: default
---

[Zhu and Stephens (*Ann. Appl. Stat.*, 2017)]: https://projecteuclid.org/euclid.aoas/1507168840
[Zhu and Stephens (*bioRxiv*, 2017)]: https://doi.org/10.1101/160770 
[`rss_bvsr.m`]: https://github.com/stephenslab/rss/blob/master/src/rss_bvsr.m
[Guan and Stephens (*Ann. Appl. Stat.*, 2011)]: https://projecteuclid.org/euclid.aoas/1318514285
[`rss_bslmm.m`]: https://github.com/stephenslab/rss/blob/master/src/rss_bslmm.m
[Zhou, Carbonetto and Stephens (*PLoS Genet.*, 2013)]: https://doi.org/10.1371/journal.pgen.1003264
[Stephens (*Biostatistics*, 2017)]: https://doi.org/10.1093/biostatistics/kxw041
[`rss_varbvsr.m`]: https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr.m
[`rss_varbvsr_squarem.m`]: https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_squarem.m
[`rss_varbvsr_bigmem_squarem.m`]: https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem_squarem.m

## Markov chain Monte Carlo (MCMC)

#### [`rss_bvsr.m`][]

Fit a Bayesian model that consists of RSS likelihood
and prior introduced in "Bayesian variable selection regression"
(BVSR; [Guan and Stephens (*Ann. Appl. Stat.*, 2011)][])
using a Metropolis-Hastings algorithm.

#### [`rss_bslmm.m`][]

Fit a Bayesian model that consists of RSS likelihood
and prior introduced in "Bayesian sparse linear mixed model"
(BSLMM; [Zhou, Carbonetto and Stephens (*PLoS Genet.*, 2013)][])
using a component-wise MCMC algorithm.

#### [`rss_ash.m`](https://github.com/stephenslab/rss/blob/master/src/rss_ash.m)

Fit a Bayesian model that consists of RSS likelihood
and prior introduced in "Adaptive shrinkage" (ASH; [Stephens (*Biostatistics*, 2017)][])
using a component-wise MCMC algorithm.

Details of MCMC algorithms are available in
[Supplementary Appendix B](http://stephenslab.uchicago.edu/assets/papers/Zhu2017-supplement.pdf)
of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][]
(or [here](http://www.stat.uchicago.edu/~xiangzhu/rss_mcmc.pdf)).

Note that [`rss_bvsr.m`][] and [`rss_bslmm.m`][] were used to
generate results in [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][].   

## Variational Bayes (VB)

#### [`rss_varbvsr.m`][]

Fit a Bayesian model that consists of RSS likelihood and BVSR prior
using a mean-field VB algorithm. The VB algorithm largely follows
[Carbonetto and Stephens (*Bayesian Anal.*, 2012)](https://projecteuclid.org/euclid.ba/1339616726).

#### [`rss_varbvsr_squarem.m`][]

Fit a Bayesian model that consists of RSS likelihood and BVSR prior
using a mean-field VB algorithm and a SQUAREM accelerator introduced in
[Varadhan and Roland (*Scand. Stat. Theory. Appl.*, 2008)](https://doi.org/10.1111/j.1467-9469.2007.00585.x).

#### [`rss_varbvsr_parallel.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_parallel.m)

This is a parallel implementation of [`rss_varbvsr.m`][].

#### [`rss_varbvsr_pasquarem.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_pasquarem.m)

This is a parallel implementation of [`rss_varbvsr_squarem.m`][].

#### [`rss_varbvsr_bigmem.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem.m)

This is a memory-efficient implementation of
[`rss_varbvsr_parallel.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_parallel.m).

#### [`rss_varbvsr_bigmem_squarem.m`][]

This is a memory-efficient implementation of
[`rss_varbvsr_pasquarem.m`](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_pasquarem.m).

Details of VB algorithms are available in
[Supplementary Notes](https://www.biorxiv.org/content/biorxiv/suppl/2018/07/16/160770.DC2/160770-1.pdf)
of [Zhu and Stephens (*bioRxiv*, 2017)][].

Note that [`rss_varbvsr_squarem.m`][] and [`rss_varbvsr_bigmem_squarem.m`][]
were used to generate results in [Zhu and Stephens (*bioRxiv*, 2017)][]. 

## Miscellaneous

#### [`import_1000g_vcf.sh`](https://github.com/stephenslab/rss/blob/master/misc/import_1000g_vcf.sh)

Output [1000 Genomes](http://www.internationalgenome.org/data)
phased haplotypes of a given list of SNPs in
[IMPUTE reference-panel format](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#input_options). 

#### [`compute_pve.m`](https://github.com/stephenslab/rss/blob/master/src/compute_pve.m)

Use GWAS summary data to estimate PVE (or SNP heritability),
a quantity defined in [Guan and Stephens (*Ann. Appl. Stat.*, 2011)][].
This function corresponds to Equation 3.7 in [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][].

#### [`band_storage.m`](https://github.com/stephenslab/rss/blob/master/misc/band_storage.m)

Convert a symmetric, banded matrix to a compact matrix in such a way
that only the main diagonal and the nonzero super-diagonals are stored.

#### [`find_bandwidth.m`](https://github.com/stephenslab/rss/blob/master/misc/find_bandwidth.m)

Find the bandwidth of a symmetric, banded matrix.

#### [`get_corr.m`](https://github.com/stephenslab/rss/blob/master/misc/get_corr.m)

Compute linkage disequilibrium (LD) matrix using the shrinkage estimator proposed in
[Wen and Stephens (*Ann. Appl. Stat.*, 2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081).
This function is also implemented in an `R` package
[LDshrink](https://github.com/stephenslab/LDshrink).

#### [`data_maker.m`](https://github.com/stephenslab/rss/blob/master/misc/data_maker.m)

Simulate phenotype data from the genome-wide multiple-SNP model described in
[Zhou, Carbonetto and Stephens (*PLoS Genet.*, 2013)][],
and then compute the single-SNP summary statistics for each SNP.
This function was used in simulation studies of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][].  

#### [`enrich_datamaker.m`](https://github.com/stephenslab/rss/blob/master/misc/enrich_datamaker.m)

Simulate phenotype data from the genetic association enrichment model described in
[Carbonetto and Stephens (*PLoS Genet.*, 2013)](https://doi.org/10.1371/journal.pgen.1003770),
and then compute the single-SNP summary statistics for each SNP.
This function was used in simulation studies of [Zhu and Stephens (*bioRxiv*, 2017)][].

#### [`null_single.m`](https://github.com/stephenslab/rss/blob/master/src_vb/null_single.m) & [`null_template.m`](https://github.com/stephenslab/rss/blob/master/src_vb/null_template.m)

Fit genome-wide multiple-SNP "baseline model" to single-SNP summary data, using
[`rss_varbvsr*` functions](https://github.com/stephenslab/rss/tree/master/src_vb).
These scripts were used in data analyses of [Zhu and Stephens (*bioRxiv*, 2017)][].

#### [`gsea_wrapper.m`](https://github.com/stephenslab/rss/blob/master/src_vb/gsea_wrapper.m) & [`gsea_template.m`](https://github.com/stephenslab/rss/blob/master/src_vb/gsea_template.m)

Fit genome-wide multiple-SNP "enrichment model" to single-SNP summary data, using
[`rss_varbvsr*` functions](https://github.com/stephenslab/rss/tree/master/src_vb).
These scripts were used in data analyses in [Zhu and Stephens (*bioRxiv*, 2017)][].

#### [`ash_lrt_31traits.R`](https://github.com/stephenslab/rss/blob/master/misc/ash_lrt_31traits.R)

Compute an enrichment likelihood ratio as a sanity check for the Bayesian model-based
enrichment analysis method developed in [Zhu and Stephens (*bioRxiv*, 2017)][].
