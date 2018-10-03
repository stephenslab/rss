---
layout: default
---

[Wen and Stephens (2010)]: https://www.ncbi.nlm.nih.gov/pubmed/21479081
[Li and Stephens (2003)]: https://www.ncbi.nlm.nih.gov/pubmed/14704198
[get_corr.m]: https://github.com/stephenslab/rss/blob/master/misc/get_corr.m
[Nick]: https://github.com/CreRecombinase
[rssr]: https://github.com/stephenslab/rssr
[LDshrink]: https://github.com/stephenslab/LDshrink
[Zhu and Stephens (2018)]: https://doi.org/10.1101/160770
[ash_lrt_31traits.R]: https://github.com/stephenslab/rss/blob/master/misc/ash_lrt_31traits.R
[Supplementary Figure 17]: https://www.biorxiv.org/content/biorxiv/suppl/2018/07/16/160770.DC2/160770-3.pdf
[compute_pip.m]: https://github.com/stephenslab/rss/blob/master/src_vb/compute_pip.m  

## General information

#### Q: RSS is cool, but it is implemented in MATLAB, which is a proprietary programming language. It would be nice if RSS methods could be made to work with some open source language so that everyone can run it freely.

A: To further improve the accessibility of RSS methods,
I am currently working on an R package [rssr][]
with [Nick][] (Nicholas Knoblauch) in my spare time.
Comments and suggestions are welcome!

## RSS methods & applications

#### Q: If I want to use RSS methods, how should I prepare for the input data?

A: Thanks for your interest in trying RSS methods!
To help you prepare your own data for RSS-based analyses,
I write a short [tutorial](Input-Data-Formats) about input data formats.
Please note that you need to store your input data in MAT-files if
you plan to use the MATLAB implementation of RSS methods.

#### Q: [Zhu and Stephens (2018)][] described a simple likelihood ratio calculation as a sanity check for the more sophisticated enrichment analysis based on RSS. Do you have software for this sanity check?

A: Yes. I write a stand-alone script [ash_lrt_31traits.R][] for this sanity check.
Please carefully read the instruction in this script.
For more details of this sanity check, please see the caption of
[Supplementary Figure 17][] in [Zhu and Stephens (2018)][].

#### Q: [Zhu and Stephens (2018)][] described an integrated framework based on RSS with both enrichment and prioritization components, and we are only interested in the gene prioritization component. How do we perform the gene-level association test using RSS (something like Figure 3 in the paper)? 

A: [Example 5 Part A](Example-5A) gives a step-by-step illustration of
gene-level association methods built on RSS, with input and output files included.
Since this example is based on a relatively small simulation,
it is highly recommended to run the codes and see how the methods work.

## Linkage disequilibrium (LD)

#### Q: RSS uses the shrinkage estimator of LD in [Wen and Stephens (2010)][], the computation of which requires the "scaled population recombination rate" (`rho_ij`). How do I compute `rho_ij`?

A: Details of calculating `rho_ij` are provided in these two great papers:
[Li and Stephens (2003)][], [Wen and Stephens (2010)][].
I also write a short [tutorial](Recombination) to illustrate the calculation of
scaled population recombination rate (`rho_ij`) using HapMap genetic map.

#### Q: Most published RSS analyses used phased haplotypes to estimate LD. Can I use unphased genotype data to compute this shrinkage estimator of LD?

A: Yes. You can feed [get_corr.m][] an unphased genotype matrix and specify that `isgeno=true`.
There is a "0.5" correction term for genotype data, which comes from a random mating assumption.
Please see Section 2.4 of [Wen and Stephens (2010)][] for more details.

#### Q: I know of [Wen and Stephens (2010)][], but I'm not interested in imputation. Do you have software that just computes the banded covariance or LD matrix?

A: Yes. If you have MATLAB available, you can use [get_corr.m][].
If you prefer open source language such as R, [Nick][] and I are working
on an R package [LDshrink][] that can produce the banded LD matrix.

## Preprocessed data
