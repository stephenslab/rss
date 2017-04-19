# Regression with Summary Statistics (RSS)

### Overview
Multiple regression analyses often assume that the response and covariates of each individual are observed, and use them to infer the regression coefficients. Here, motivated by the applications in genetics, we assume that these individual-level data are not available, but instead the summary statistics of univariate regression (essentially, the effect size estimates and their standard errors) are provided. We also assume that information on the correlation structure among covariates is available. The aim is to infer the multiple regression coefficients using the marginal regression summary statistics.

This work is motivated by applications in genome-wide association studies ([GWAS](https://en.wikipedia.org/wiki/Genome-wide_association_study)). When fitting the multiple regression model to GWAS individual-level data, the covariates are the genotypes typed at different genetic variants ([SNPs](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism)), the response is the quantitative phenotype (e.g. height), and the regression coefficients are the effects of each SNP on phenotype. Due to privacy and logistical issues, the individual-level data are often not easily available. In contrast, the GWAS summary statistics (from single-SNP analysis) are widely available in the public domain (e.g. [GIANT](https://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files) and [PGC](https://www.med.unc.edu/pgc/downloads)). Moreover, the correlation among covariantes (SNPs), known as [linkage disequilibrium](https://en.wikipedia.org/wiki/Linkage_disequilibrium), also can be obtained from public databases (e.g. the [1000 Genomes Project](http://www.1000genomes.org/home)). When the individual-level data are not available, can we perform "multiple-SNP" analysis using these public assets?

Here we provide a generally-applicable framework for the multiple-SNP analyses using GWAS summary data. Specifically, we introduce a “Regression with Summary Statistics” (RSS) likelihood, which relates the multiple regression coefficients to univariate regression results. We then combine the RSS likelihood with suitable priors to perform Bayesian inference for the regression coefficients.

### License 
Distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The repository is distributed in the hope that it will be useful, but **without any warranty**; without even the implied warranty of **merchantability** or **fitness for a particular purpose**. See [LICENSE](LICENSE) for more details.

### Support
1. Get started from some short [tutorials](http://stephenslab.github.io/rss).
2. Refer to the [FAQ](https://github.com/stephenslab/rss/wiki/FAQ) page for answers to some common questions.
3. Create a new [issue](https://github.com/stephenslab/rss/issues) to report bugs and/or request features.
4. Send an email to `xiangzhu[at]uchicago.edu`.

### Citation
- The Regression with Summary Statistics (RSS) likelihood <br> Xiang Zhu and Matthew Stephens (2016). [Bayesian large-scale multiple regression with summary statistics from genome-wide association studies](https://doi.org/10.1101/042457). bioRxiv. (*Annals of Applied Statistics* To appear.)

- Gene set enrichment analysis based on RSS <br> Manuscript in preparation. [Online Notebook.](http://xiangzhu.github.io/rss-gsea/_book/)

### Collaboration
This work derives a likelihood of multiple regression coefficients based on univariate regression summary data, which opens the door to a wide range of statistical machinery for inference. Using this likelihood, we implement Bayesian methods to estimate SNP heritability, detect genetic association, perform gene set enrichment analysis, etc. The update on our progress can be found in [NEWS](NEWS.md). 

If you have specific applications that use GWAS summary data as input, and want to build new methods based on the RSS likelihood, please feel free to contact us. We are glad to help!  

### Contact
[Xiang Zhu](https://github.com/xiangzhu) <br>
[Matthew Stephens Lab](http://stephenslab.uchicago.edu) <br>
[Department of Statistics](https://galton.uchicago.edu) <br>
[University of Chicago](https://www.uchicago.edu) <br>

