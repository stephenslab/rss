---
layout: default
---

[Zhu and Stephens (2017)]: https://projecteuclid.org/euclid.aoas/1507168840
[Zhu and Stephens (2018)]: https://www.nature.com/articles/s41467-018-06805-x 
[rss_bvsr.m]: https://github.com/stephenslab/rss/blob/master/src/rss_bvsr.m
[(Guan and Stephens, 2011)]: https://projecteuclid.org/euclid.aoas/1318514285
[Guan and Stephens (2011)]: https://projecteuclid.org/euclid.aoas/1318514285
[rss_bslmm.m]: https://github.com/stephenslab/rss/blob/master/src/rss_bslmm.m
[(Zhou et al, 2013)]: https://doi.org/10.1371/journal.pgen.1003264
[Zhou et al (2013)]: https://doi.org/10.1371/journal.pgen.1003264
[rss_ash.m]: https://github.com/stephenslab/rss/blob/master/src/rss_ash.m
[(Stephens, 2017)]: https://doi.org/10.1093/biostatistics/kxw041
[Carbonetto and Stephens (2012)]: https://projecteuclid.org/euclid.ba/1339616726
[rss_varbvsr.m]: https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr.m
[rss_varbvsr_squarem.m]: https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_squarem.m
[(Varadhan and Roland, 2008)]: https://doi.org/10.1111/j.1467-9469.2007.00585.x
[rss_varbvsr_bigmem_squarem.m]: https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem_squarem.m

# Markov chain Monte Carlo (MCMC)

## [rss_bvsr.m][]

Fit the RSS-BVSR model that consists of RSS likelihood
and BVSR prior [(Guan and Stephens, 2011)][]:

$$
\begin{aligned}
\widehat{\boldsymbol\beta} &\sim {\cal N}({\bf SRS}^{-1}{\boldsymbol\beta},{\bf SRS}),\\
\beta_j &\sim \pi\cdot{\cal N}(0,\sigma_B^2) + (1-\pi)\cdot\delta_0,
\end{aligned}
$$

using a Metropolis-Hastings algorithm.

## [rss_bslmm.m][]

Fit the RSS-BSLMM model that consists of RSS likelihood
and BSLMM prior [(Zhou et al, 2013)][]:

$$
\begin{aligned}
\widehat{\boldsymbol\beta} &\sim {\cal N}({\bf SRS}^{-1}{\boldsymbol\beta},{\bf SRS}),\\
\beta_j &\sim \pi\cdot{\cal N}(0,\sigma_B^2+\sigma_P^2) + (1-\pi)\cdot{\cal N}(0,\sigma_P^2),
\end{aligned}
$$

using a component-wise MCMC algorithm.

## [rss_ash.m][]

Fit the RSS-ASH model that consists of RSS likelihood
and ASH prior [(Stephens, 2017)][]:

$$
\begin{aligned}
\widehat{\boldsymbol\beta} &\sim {\cal N}({\bf SRS}^{-1}{\boldsymbol\beta},{\bf SRS}),\\
\beta_j &\sim \pi_0 \cdot \delta_0 + {\textstyle\sum}_{k=1}^K \pi_k \cdot {\cal N}(0,\sigma_k^2),
\end{aligned}
$$

using a component-wise MCMC algorithm.

Details of MCMC algorithms for [rss_bvsr.m][] and [rss_bslmm.m][] are
available in Supplementary Appendix B of [Zhu and Stephens (2017)][]
Details of [rss_ash.m][] are available
[here](http://www.stat.uchicago.edu/~xiangzhu/rss_mcmc.pdf).
Note that only [rss_bvsr.m][] and [rss_bslmm.m][] were used to
generate results in [Zhu and Stephens (2017)][].   

# Variational Bayes (VB)

## [rss_varbvsr.m][]

Fit the following extended RSS-BVSR model

$$
\begin{aligned}
\widehat{\boldsymbol\beta} &\sim {\cal N}({\bf SRS}^{-1}{\boldsymbol\beta},{\bf SRS}),\\
\beta_j &\sim \pi_j\cdot{\cal N}(0,\sigma_j^2) + (1-\pi_j)\cdot\delta_0,
\end{aligned}
$$

using a mean-field VB algorithm.
The VB algorithm largely follows [Carbonetto and Stephens (2012)][].
This is an extended RSS-BVSR model because each SNP $j$ can have
its own hyper-parameters $$\{\pi_j,\sigma_j^2\}$$,
whereas the standard RSS-BVSR model assumes that all SNPs share
the same hyper-parameters $$\{\pi,\sigma_B^2\}$$.

## [rss_varbvsr_squarem.m][]

This is a variant of [rss_varbvsr.m][] with the SQUAREM
accelerator [(Varadhan and Roland, 2008)][] added.

## [rss_varbvsr_parallel.m](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_parallel.m)

This is a parallel implementation of [rss_varbvsr.m][].

## [rss_varbvsr_pasquarem.m](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_pasquarem.m)

This is a parallel implementation of [rss_varbvsr_squarem.m][].

## [rss_varbvsr_bigmem.m](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem.m)

This is a memory-efficient implementation of
[rss_varbvsr_parallel.m](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_parallel.m).

## [rss_varbvsr_bigmem_squarem.m][]

This is a memory-efficient implementation of
[rss_varbvsr_pasquarem.m](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_pasquarem.m).

Details of VB algorithms, SQUAREM accelerator and parallel implementation
are available in Supplementary Notes of [Zhu and Stephens (2018)][].
Note that only [rss_varbvsr_squarem.m][] and [rss_varbvsr_bigmem_squarem.m][]
were used to generate results in [Zhu and Stephens (2018)][].
The other functions were developed merely for testing and benchmarking. 

# Miscellaneous

## [import_1000g_vcf.sh](https://github.com/stephenslab/rss/blob/master/misc/import_1000g_vcf.sh)

Output [1000 Genomes](http://www.internationalgenome.org/data)
phased haplotypes of a given list of SNPs in
[IMPUTE reference-panel format](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#input_options). 

## [compute_pve.m](https://github.com/stephenslab/rss/blob/master/src/compute_pve.m)

Use GWAS summary data to estimate PVE (or SNP heritability),
a quantity defined by Equation 2.10 in [Guan and Stephens (2011)][].
This function corresponds to Equation 3.7 in [Zhu and Stephens (2017)][].

## [band_storage.m](https://github.com/stephenslab/rss/blob/master/misc/band_storage.m)

Convert a symmetric, banded matrix to a compact matrix in such a way
that only the main diagonal and the nonzero super-diagonals are stored.
This function is used to reduce the file size of a large LD matrix.

## [find_bandwidth.m](https://github.com/stephenslab/rss/blob/master/misc/find_bandwidth.m)

Find the bandwidth of a symmetric, banded matrix.

## [get_corr.m](https://github.com/stephenslab/rss/blob/master/misc/get_corr.m)

Compute linkage disequilibrium (LD) matrix using the shrinkage estimator proposed in
[Wen and Stephens (2010)](https://www.ncbi.nlm.nih.gov/pubmed/21479081).
This function is also implemented in an R package
[ldshrink](https://github.com/stephenslab/ldshrink).

## [data_maker.m](https://github.com/stephenslab/rss/blob/master/misc/data_maker.m)

Simulate phenotype data from the genome-wide multiple-SNP model
described in [Zhou et al (2013)][],
and then compute the single-SNP summary statistics for each SNP.
This function was used in some simulation studies of [Zhu and Stephens (2017)][].  

## [enrich_datamaker.m](https://github.com/stephenslab/rss/blob/master/misc/enrich_datamaker.m)

Simulate phenotype data from the genetic association enrichment model described in
[Carbonetto and Stephens (2013)](https://doi.org/10.1371/journal.pgen.1003770),
and then compute the single-SNP summary statistics for each SNP.
This function was used in some simulation studies of [Zhu and Stephens (2018)][].

## [null_single.m](https://github.com/stephenslab/rss/blob/master/src_vb/null_single.m) & [null_template.m](https://github.com/stephenslab/rss/blob/master/src_vb/null_template.m)

Fit genome-wide multiple-SNP "baseline model" to single-SNP summary data, using
[rss_varbvsr functions](https://github.com/stephenslab/rss/tree/master/src_vb).
These scripts were used in data analyses of [Zhu and Stephens (2018)][].

## [gsea_wrapper.m](https://github.com/stephenslab/rss/blob/master/src_vb/gsea_wrapper.m) & [gsea_template.m](https://github.com/stephenslab/rss/blob/master/src_vb/gsea_template.m)

Fit genome-wide multiple-SNP "enrichment model" to single-SNP summary data, using
[rss_varbvsr functions](https://github.com/stephenslab/rss/tree/master/src_vb).
These scripts were used in data analyses of [Zhu and Stephens (2018)][].

## [null_wrapper_fixsb.m](https://github.com/stephenslab/rss/blob/master/src_vb/null_wrapper_fixsb.m) & [gsea_wrapper_fixsb.m](https://github.com/stephenslab/rss/blob/master/src_vb/gsea_wrapper_fixsb.m)

Fit genome-wide multiple-SNP "baseline model" and "enrichment model" to single-SNP summary data,
using a fixed prior variance of causal genetic effects (`sigb` square) in
[rss_varbvsr functions](https://github.com/stephenslab/rss/tree/master/src_vb).
These scripts were used in simulation studies of [Zhu and Stephens (2018)][].

## [ash_lrt_31traits.R](https://github.com/stephenslab/rss/blob/master/misc/ash_lrt_31traits.R)

Compute a simple likelihood ratio as a sanity check for the more
complicated enrichment analysis method developed in [Zhu and Stephens (2018)][].
This likelihood ratio calculation is based on an R package
[ashr](https://cran.r-project.org/web/packages/ashr/index.html).
