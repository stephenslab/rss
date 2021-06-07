# Regression with Summary Statistics (RSS)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1473797.svg)](https://doi.org/10.5281/zenodo.1473797)

### Overview

Multiple regression analyses often assume that the response and
covariates of each individual are observed, and use them to infer the
regression coefficients. Here, motivated by the applications in
genetics, we assume that these individual-level data are not
available, but instead the summary statistics of univariate regression
(essentially, the effect size estimates and their standard errors) are
provided. We also assume that information on the correlation structure
among covariates is available. The aim is to infer the multiple
regression coefficients using the marginal regression summary
statistics.

This work is motivated by applications in genome-wide association studies
([GWAS](https://en.wikipedia.org/wiki/Genome-wide_association_study)).
When fitting the multiple regression model to individual-level data of GWAS,
the covariates are the genotypes typed at different genetic variants
(typically [SNPs](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism)),
the response is the quantitative phenotype (e.g. height or blood lipid level),
and the regression coefficients are the effects of each SNP on phenotype.
Due to privacy and logistical issues, the individual-level data are often not easily available.
In contrast, the GWAS summary statistics (from standard single-SNP analysis)
are widely available in the public domain
(e.g. [GIANT](https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files)
and [PGC](https://www.med.unc.edu/pgc/results-and-downloads/downloads)).
Moreover, the correlation among covariates (genotypes of SNPs),
known as [linkage disequilibrium](https://en.wikipedia.org/wiki/Linkage_disequilibrium),
also can be obtained from public databases
(e.g. the [1000 Genomes Project](https://www.1000genomes.org/home)).
When the protected individual-level data are not available,
can we perform "multiple-SNP" analysis using these public assets?

Here we provide a generally-applicable framework for the
multiple-SNP analyses using GWAS single-SNP summary data.
Specifically, we introduce a “Regression with Summary Statistics” (RSS) likelihood,
which relates the multiple regression coefficients to univariate regression results.
We then combine the RSS likelihood with suitable priors to
perform Bayesian inference for the regression coefficients.

### License

The repository is licensed under the [MIT License](LICENSE).

### Support

1. Get started from some short [tutorials](http://stephenslab.github.io/rss).
2. Refer to [FAQ](https://stephenslab.github.io/rss/faq.html) for answers to some common questions.
3. Create a new [issue](https://github.com/stephenslab/rss/issues) to report bugs and/or request features.
4. Send an email to `xiangzhu[at]psu.edu`.

### Citation

- **The Regression with Summary Statistics (RSS) likelihood** <br>
Xiang Zhu and Matthew Stephens (2017).
Bayesian large-scale multiple regression with
summary statistics from genome-wide association studies.
*Annals of Applied Statistics* 11(3): 1561-1592.
[[Article PDF](https://stephenslab.uchicago.edu/assets/papers/Zhu2017.pdf)]
[[Journal Page](https://dx.doi.org/10.1214/17-AOAS1046)]
[[bioRxiv Page](https://doi.org/10.1101/042457)]
[[Supplementary Information](https://stephenslab.uchicago.edu/assets/papers/Zhu2017-supplement.pdf)]
[[Software](https://github.com/stephenslab/rss/tree/master/src)]

- **RSS-E: Enrichment and prioritization analysis based on RSS likelihood** <br>
Xiang Zhu and Matthew Stephens (2018).
Large-scale genome-wide enrichment analyses identify new
trait-associated genes and pathways across 31 human phenotypes.
*Nature Communications* 9, 4361.
[[Article PDF](https://www.nature.com/articles/s41467-018-06805-x.pdf)]
[[Journal Page](https://doi.org/10.1038/s41467-018-06805-x)]
[[bioRxiv Page](https://doi.org/10.1101/160770)]
[[Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-06805-x/MediaObjects/41467_2018_6805_MOESM1_ESM.pdf)]
[[Online Results](https://xiangzhu.github.io/rss-gsea/)]
[[Software](https://github.com/stephenslab/rss/tree/master/src_vb)]

- **RSS-NET: Integrated analysis of regulatory networks based on RSS likelihood** <br>
Xiang Zhu, Zhana Duren and Wing Hung Wong (2021).
Modeling regulatory network topology improves
genome-wide analyses of complex human traits.
*Nature Communications* 12, 2851.
[[Article PDF](https://www.nature.com/articles/s41467-021-22588-0.pdf)]
[[Journal Page](https://doi.org/10.1038/s41467-021-22588-0)]
[[bioRxiv Page](https://doi.org/10.1101/2020.03.13.990010)]
[[Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-22588-0/MediaObjects/41467_2021_22588_MOESM1_ESM.pdf)]
[[Online Results](https://xiangzhu.github.io/rss-net-results/)]
[[Software](https://github.com/SUwonglab/rss-net)]

- Genetic architecture inference of complex traits based on RSS likelihood <br> TBA
- Simple and robust heritability estimation based on RSS likelihood <br> TBA
- Cross-population genetic analysis of complex traits based on RSS likelihood <br> TBA 

### Collaboration

Here we have developed a [likelihood function](http://dx.doi.org/10.1214/17-AOAS1046)
of multiple regression coefficients based on univariate regression summary data,
which opens the door to a wide range of statistical machinery for inference.
Using this likelihood, we have implemented Bayesian methods to estimate SNP heritability,
detect genetic association, assess gene set or network enrichment,
prioritize trait-associated genes and infer genetic architecture.
Please check our [progress updates](https://stephenslab.github.io/rss/news.html) regularly. 

If you have specific applications that use GWAS summary data as input,
and want to build new statistical methods based on the RSS likelihood,
please feel free to contact us. We are glad to help!  

### Contact

[Xiang Zhu](https://github.com/xiangzhu) <br>
[Matthew Stephens Lab](http://stephenslab.uchicago.edu) <br>
[Department of Statistics](https://stat.uchicago.edu) <br>
[University of Chicago](https://www.uchicago.edu) <br>

