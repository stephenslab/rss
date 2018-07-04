# Regression with Summary Statistics (RSS)

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
In contrast, the GWAS summary statistics (from standard single-SNP analysis) are widely available in the public domain
(e.g. [GIANT](https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files)
and [PGC](https://www.med.unc.edu/pgc/results-and-downloads/downloads)).
Moreover, the correlation among covariates (genotypes of SNPs),
known as [linkage disequilibrium](https://en.wikipedia.org/wiki/Linkage_disequilibrium), also can be obtained from public databases (e.g. the [1000 Genomes Project](http://www.1000genomes.org/home)). When the individual-level data are not available, can we perform "multiple-SNP" analysis using these public assets?

Here we provide a generally-applicable framework for the multiple-SNP analyses using GWAS single-SNP summary data.
Specifically, we introduce a “Regression with Summary Statistics” (RSS) likelihood,
which relates the multiple regression coefficients to univariate regression results.
We then combine the RSS likelihood with suitable priors to perform Bayesian inference for the regression coefficients.

### License 
Distributed under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

The repository is distributed in the hope that it will be useful,
but **without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**.
Please see [LICENSE](LICENSE) for more details.

### Support
1. Get started from some short [tutorials](http://stephenslab.github.io/rss).
2. Refer to the [FAQ](http://stephenslab.github.io/rss/FAQ) page for answers to some common questions.
3. Create a new [issue](https://github.com/stephenslab/rss/issues) to report bugs and/or request features.
4. Send an email to `xiangzhu[at]uchicago.edu`.

### Citation
- The Regression with Summary Statistics (RSS) likelihood <br> Xiang Zhu and Matthew Stephens (2017).
[Bayesian large-scale multiple regression with summary statistics from genome-wide association studies](http://stephenslab.uchicago.edu/assets/papers/Zhu2017.pdf).
[*Annals of Applied Statistics* 11(3): 1561-1592](http://dx.doi.org/10.1214/17-AOAS1046). 

- Enrichment and prioritization analysis of GWAS summary data based on RSS <br> Xiang Zhu and Matthew Stephens (2017). [A large-scale genome-wide enrichment analysis identifies new trait-associated genes, pathways and tissues across 31 human phenotypes](https://doi.org/10.1101/160770). bioRxiv. [[Supplementary Materials](http://www.biorxiv.org/content/early/2017/07/08/160770.figures-only)] [[Full Results](http://xiangzhu.github.io/rss-gsea/results)]

- Inferring genetic architecture of complex human traits based on RSS <br> TBA

### Collaboration
In this project, we have derived a [likelihood function](https://en.wikipedia.org/wiki/Likelihood_function)
of multiple regression coefficients based on univariate regression summary data,
which opens the door to a wide range of statistical machinery for inference.
Using this likelihood, we have implemented Bayesian methods to estimate SNP heritability,
detect genetic association, perform gene set enrichment analysis, infer genetic architecture, etc.
Please check our [progress updates](http://stephenslab.github.io/rss/News) regularly. 

If you have specific applications that use GWAS summary data as input,
and want to build new statistical methods based on the RSS likelihood,
please feel free to contact us. We are glad to help!  

### Contact
[Xiang Zhu](https://github.com/xiangzhu) <br>
[Matthew Stephens Lab](http://stephenslab.uchicago.edu) <br>
[Department of Statistics](https://galton.uchicago.edu) <br>
[University of Chicago](https://www.uchicago.edu) <br>

