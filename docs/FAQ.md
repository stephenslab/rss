---
layout: default
---

[Wen and Stephens (*Ann. Appl. Stat.*, 2010)]: https://www.ncbi.nlm.nih.gov/pubmed/21479081
[Li and Stephens (*Genetics*, 2003)]: https://www.ncbi.nlm.nih.gov/pubmed/14704198
[`get_corr.m`]: https://github.com/stephenslab/rss/blob/master/misc/get_corr.m
[Nick]: https://github.com/CreRecombinase
[`rssr`]: https://github.com/stephenslab/rssr
[`LDshrink`]: https://github.com/stephenslab/LDshrink
[Zhu and Stephens (*bioRxiv*, 2017)]: https://doi.org/10.1101/160770
[`ash_lrt_31traits.R`]: https://github.com/stephenslab/rss/blob/master/misc/ash_lrt_31traits.R

## General information

#### Q: I can't afford [`MATLAB`](https://www.mathworks.com/). It would be nice if RSS methods could be made to work with some open source language so everyone can run it.

A: To further improve the accessibility of RSS methods,
I am working on an `R` package [`rssr`][]
with [Nick][] (Nicholas Knoblauch) in my spare time.
Comments and suggestions are welcome!

## Linkage disequilibrium (LD)

#### Q: The shrinkage estimator of LD in [Wen and Stephens (*Ann. Appl. Stat.*, 2010)][] requires the "scaled population recombination rate" (`rho_ij`). How do I compute it?

A: Details of calculating `rho_ij` are provided in these two great papers:
[Li and Stephens (*Genetics*, 2003)][] and [Wen and Stephens (*Ann. Appl. Stat.*, 2010)][].
I also write a short [tutorial](Recombination) to illustrate the calculation of
scaled population recombination rate (`rho_ij`) using HapMap genetic map.

#### Q: Can I use unphased genotype data to compute the shrinkage estimator?

A: Yes. You can input an unphased genotype matrix as `Hpanel` in [`get_corr.m`][]
and specify that `isgeno=true`.
With `isgeno=true`, [`get_corr.m`][] will execute the following code chunk:

```matlab
if isgeno
  disp('Hpanel is an unphased genotype matrix.');
  S = 0.5*S;
end
```

The "0.5" correction term comes from a random mating assumption. 
Please see Section 2.4 of [Wen and Stephens (*Ann. Appl. Stat.*, 2010)][] for more details.

#### Q: I know of [Wen and Stephens (*Ann. Appl. Stat.*, 2010)][], but I'm not interested in imputation. Do you have software that just computes the banded covariance matrix?

A: Yes. If you have `MATLAB` available, you can use [`get_corr.m`][].
If you prefer open source language like `R`, [Nick][] and I are working
on an `R` package [`LDshrink`][] to compute the banded LD matrix. 

## Applications

#### Q: If I want to use RSS methods, how should I prepare for the input data?

A: Thanks for your interest in trying RSS methods!
To help you prepare your own data for RSS-based analyses,
I write a short [tutorial](Input-Data-Formats) about input data formats.

#### Q: [Zhu and Stephens (*bioRxiv*, 2017)][] described a simple likelihood ratio calculation as a sanity check. Do you have software for this sanity check?

A: Yes. I write a stand-alone script [`ash_lrt_31traits.R`][] for this sanity check.
Please carefully read the instruction in this script.
For more details of this sanity check, please see the caption of
[Supplementary Figure 17](https://www.biorxiv.org/content/biorxiv/suppl/2018/07/16/160770.DC2/160770-3.pdf)
in [Zhu and Stephens (*bioRxiv*, 2017)][].
