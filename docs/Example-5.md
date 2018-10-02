[Zhu and Stephens (2018)]: https://doi.org/10.1101/160770
[Carbonetto and Stephens (2013)]: https://doi.org/10.1371/journal.pgen.1003770
[example5]: https://github.com/stephenslab/rss/tree/master/examples/example5

# Example 5: Enrichment and prioritization analysis of GWAS summary statistics using RSS.

## Overview

This example illustrates how to perform enrichment and prioritization analysis of
GWAS summary statistics based on variational Bayes (VB) inference of RSS-BVSR model.
This example consists of:

- [Part A](Example-5A): enrichment and prioritization analysis of
a synthetic dataset used in simulation studies of [Zhu and Stephens (2018)][];
- [Part B](Example-5B): enrichment and prioritization analysis of
inflammatory bowel disease GWAS summary statistics
[(Liu et al, 2015)](https://www.ncbi.nlm.nih.gov/pubmed/26192919)
and a gene set named *IL23-mediated signaling events* (Pathway Commons 2, PID, 37 genes).

[Part A](Example-5A) provides a quick view of how RSS works in enrichment and prioritization analysis.
[Part B](Example-5B) illustrates the actual data analyses performed in [Zhu and Stephens (2018)][].

## Details

The following figure provides a schematic overview of the method.

<img src="images/rss_gsea.png" width="200" />

As shown above, RSS fits two models for the enrichment and prioritization analysis.

- Baseline model (M0): SNPs across the whole genome are equally
likely to be associated with the phenotype of interest. 
- Enrichment model (M1): SNPs "inside" a gene set are more likely
(i.e. "enriched") to be associated with a target phenotype than remaining SNPs.

If the gene set is truly enriched, then the observed GWAS data
should show more support for the enrichment over baseline model,
that is, yielding a larger Bayes factor (BF).

In addition to identifying enrichments, RSS also automatically prioritizes
loci within an enriched set by comparing the posterior distributions of
genetic effects (beta) under M0 and M1.
Here we summarize the posterior of beta as P1,
the posterior probability that at least one SNP in a locus is trait-associated.
Differences between P1 estimated under M0 and M1 reflect
the influence of enrichment on genetic associations,
which can help identify new trait-associated genes. 

The key difference between RSS and previous work,
notably, [Carbonetto and Stephens (2013)][], is that
RSS uses publicly available GWAS summary data,
rather than individual-level genotypes and phenotypes.
To perform similar analysis of GWAS individual-level data,
please see <https://github.com/pcarbo/bmapathway>.

To reproduce results of Example 5,
please use scripts in the directory [example5][],
and follow the step-by-step guide below.

- Part A: <http://stephenslab.github.io/rss/Example-5A>
- Part B: <http://stephenslab.github.io/rss/Example-5B>

Before running either part, please make sure the
[VB subroutines](https://github.com/stephenslab/rss/tree/master/src_vb)
of RSS are installed. Please find installation instructions [here](RSS-via-VB).

It is advisable to go through the simulated example ([Part A](Example-5A))
before diving into the real data example ([Part B](Example-5B)).
