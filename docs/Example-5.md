[Zhu and Stephens (*bioRxiv*, 2017)]: https://doi.org/10.1101/160770
[Carbonetto and Stephens (*PLoS Genet.*, 2013)]: http://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1003770
[example5]: https://github.com/stephenslab/rss/tree/master/examples/example5

# Example 5: Enrichment analysis of GWAS summary statistics using RSS.

## Overview

This example illustrates how to perform enrichment analysis of
GWAS summary statistics based on variational Bayes (VB) inference of RSS-BVSR model.
This example consists of two parts:

- [Part A](Example-5A). An enrichment analysis of a synthetic dataset used in
simulation studies of [Zhu and Stephens (*bioRxiv*, 2017)][].
This part gives you a quick view of how RSS works in an enrichment analysis.
- [Part B](Example-5B). An end-to-end enrichment analysis of inflammatory bowel disease GWAS summary statistics
([Liu et al, *Nat Genet.*, 2015](https://www.ncbi.nlm.nih.gov/pubmed/26192919)) and
a gene set named *IL23-mediated signaling events* (Pathway Commons 2, PID, 37 genes) using RSS.
This part illustrates the actual data analyses performed in [Zhu and Stephens (*bioRxiv*, 2017)][].

## Details

The following figure provides a schematic overview of the analysis method.
For more details, please see [Zhu and Stephens (*bioRxiv*, 2017)][].

![](images/rss_gsea.png)

Here the enrichment analysis consists of fitting the following two models:

- **baseline model**: SNPs across the genome are equally likely to be associated with a target phenotype. 
- **enrichment model**: SNPs "inside" a gene in a gene set are more likely (i.e. "enriched") to be
associated with a target phenotype than remaining SNPs. 

The key difference between RSS and previous work
notably, [Carbonetto and Stephens (*PLoS Genet.*, 2013)][], is that
RSS uses **publicly available** GWAS summary data, rather than individual-level genetic data.

To reproduce results of Example 5,
please use scripts in the directory [example5][],
and follow the step-by-step guide in each part.

- Part A: [http://stephenslab.github.io/rss/Example-5A](Example-5A)
- Part B: [http://stephenslab.github.io/rss/Example-5B](Example-5B)

Before running either part, please make sure the
[VB subroutines](https://github.com/stephenslab/rss/tree/master/src_vb)
of RSS are installed. Please find installation instructions [here](RSS-via-VB).

It is advisable to go through the simulated example in [Example 5 Part A](Example-5A)
before diving into the real data example in [Example 5 Part B](Example-5B).
