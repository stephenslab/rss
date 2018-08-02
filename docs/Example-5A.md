[Zhu and Stephens (*bioRxiv*, 2017)]: https://doi.org/10.1101/160770
[`example5_simulated.m`]: https://github.com/stephenslab/rss/blob/master/examples/example5/example5_simulated.m
[`example5_simulated_data.mat`]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/example5_simulated_data.mat
[`example5_simulated_results.mat`]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/example5_simulated_results.mat
[Wellcome Trust Case Control Consortium, *Nature*, 2007]: https://www.ncbi.nlm.nih.gov/pubmed/17554300
[`null_wrapper_fixsb.m`]: https://github.com/stephenslab/rss/blob/master/src_vb/null_wrapper_fixsb.m
[`gsea_wrapper_fixsb.m`]: https://github.com/stephenslab/rss/blob/master/src_vb/gsea_wrapper_fixsb.m

# Example 5: Enrichment analysis of GWAS summary statistics using RSS (Part A).

## Overview

This is Part A of [Example 5](Example-5),
which illustrates how to perform enrichment analysis of
GWAS summary statistics based on variational Bayes (VB) inference of RSS-BVSR model.
This part describes an enrichment analysis of a synthetic dataset used in
simulation studies of [Zhu and Stephens (*bioRxiv*, 2017)][].
This part gives you a quick view of how RSS works in an enrichment analysis.

The input dataset here is simulated under the enrichment model of RSS,
using real genotypes of 12,758 SNPs on chromosome 16 from 1458 individuals
in the UK Blood Service Control Group ([Wellcome Trust Case Control Consortium, *Nature*, 2007][]).
Please see the caption of
[Supplementary Figure 1](https://www.biorxiv.org/content/biorxiv/suppl/2018/07/16/160770.DC2/160770-3.pdf)
in [Zhu and Stephens (*bioRxiv*, 2017)][] for the simulation details. 

To reproduce results of Example 5 Part B,
please use the script [`example5_simulated.m`][],
and follow the step-by-step guide below.
Before running [`example5_simulated.m`][], please make sure the
[VB subroutines](https://github.com/stephenslab/rss/tree/master/src_vb) of RSS are installed.
Please find installation instructions [here](RSS-via-VB).

Once the software is installed and the input data is downloaded,
one should be able to run this part by simply
typing the following line in a `Matlab` console:

```matlab
>> run example5_simulated.m
```

## Step-by-step illustration

**Step 1**. Download the input data file [`example5_simulated_data.mat`][],
which contains the GWAS summary statistics and LD matrix estimates.
Please contact me if you have trouble accessing this file.

The data file [`example5_simulated_data.mat`][] contains the following elements.

```matlab
>> example_data = matfile('example5_simulated_data.mat');
>> example_data

example_data =
  matlab.io.MatFile
  Properties:
      Properties.Source: 'example5_simulated_data.mat'
    Properties.Writable: false
                      R: [12758x12758 double]
                betahat: [12758x1     double]
                     se: [12758x1     double]
                   snps: [676x1       double]
```

- `R`: 12758 by 12758 matrix, LD matrix estimated from a reference panel
- `betahat`: 12758 by 1 vector, single-SNP effect size estimate for each SNP
- `se`: 12758 by 1 vector, standard errors of the single-SNP effect size estimates
- `snps`: 676 by 1 vector, indices of SNPs that are "inside" the target pathway

Typically RSS only requires these four input variables for an enrichment analysis.
To further reduce computation, RSS uses the matrix `SiRiS` instead of `R`:

```matlab
p     = length(betahat);
Si    = 1 ./ se(:);
SiRiS = sparse(repmat(Si,1,p) .* R .* repmat(Si',p,1));
clear Si R;
```

**Step 2**. Specify the grid for hyper-parameters.
The total computational cost of RSS is proportional to the grid size,
and thus please consider reducing the size of `theta0` and/or `theta`
if you want to finish running [`example5_simulated.m`][] faster.
(It takes 4 hours to complete the analysis on a single CPU using the grid below.)

```matlab
% specify hyper-parameters
theta0 = (-4.5:0.05:-3.5)';      % grid for the genome-wide log-odds (base 10)
theta  = (1.5:0.05:2.5)';        % grid for the log-fold enrichment (base 10)
sigb   = 1;                      % prior SD of genetic effects
```

**Step 3**. Initialize the variational parameters.
As in [Zhu and Stephens (*bioRxiv*, 2017)][], here we use a simple random start
to set initial values of variational parameters `{alpha,mu}` for baseline models.

```matlab
% initialize the variational parameters
myseed = 200;

rng(myseed, 'twister');
alpha0 = rand(p,1);
alpha0 = alpha0 ./ sum(alpha0);

rng(myseed+1, 'twister');
mu0 = randn(p,1);

n0         = length(theta0);
alpha0_rss = repmat(alpha0, [1 n0]);
mu0_rss    = repmat(mu0, [1 n0]);
```

**Step 4**. Fit the baseline model.
Since we set `sigb=1` in Step 2, we use a wrapper [`null_wrapper_fixsb.m`][]
that fixes the value of `sigb` for all elements in `theta0`.

```matlab
[b_logw,b_alpha,b_mu,b_s] = null_wrapper_fixsb('squarem',betahat,se,SiRiS,sigb,theta0,alpha0_rss,mu0_rss);
```

**Step 5**. Fit the enrichment model.
Since we set `sigb=1` in Step 2, we use a wrapper [`gsea_wrapper_fixsb.m`][]
that fixes the value of `sigb` for all elements in `theta0` and `theta`.

```matlab
[log10_bf,e_logw,e_alpha,e_mu,e_s] = gsea_wrapper_fixsb('squarem',betahat,se,SiRiS,snps,sigb,theta0,theta,b_logw,b_alpha,b_mu);
```

Note that here we use the variational parameter estimates
from the baseline model fitting `{b_alpha,b_mu}` to set
initial values of variational parameters for enrichment models.

**Step 6**. Save the analysis results.
The analysis results are saved as [`example5_simulated_results.mat`][].

```matlab
>> load example5_simulated_results.mat
>> whos
  Name              Size                  Bytes  Class     Attributes

  b_alpha       12758x21                2143344  double
  b_logw           21x1                     168  double
  b_mu          12758x21                2143344  double
  b_s           12758x21                2143344  double
  e_alpha       12758x21x21            45010224  double
  e_logw           21x21                   3528  double
  e_mu          12758x21x21            45010224  double
  e_s           12758x21x21            45010224  double
  log10_bf          1x1                       8  double
  rss_time          1x1                       8  double
```

There are three groups of output variables in the result file.

The first group consists of estimated variational parameters
`b_alpha` & `b_mu` & `b_s` under the baseline model,
where `b_alpha(:,i)` & `b_mu(:,i)` & `b_s(:,i)` correspond to
estimation under the hyper-parameter setting `[sigb, theta0(i)]`.

The second group consists of estimated variational parameters
`e_alpha` & `e_mu` & `e_s` under the enrichment model,
where `e_alpha(:,i,j)` & `e_mu(:,i,j)` & `e_s(:,i,j)` correspond to
estimation under the hyper-parameter setting `[sigb, theta0(i), theta(j)]`.

Note that contrasting `{b_alpha,b_mu}` with `{e_alpha,e_mu}` can further
help prioritize genetic associations in light of inferred enrichments.
Please see [Zhu and Stephens (*bioRxiv*, 2017)][] for more details.

The third group consists of estimated variational lower bounds and Bayes factor.

- `b_logw`: variational lower bounds under the baseline model.
- `e_logw`: variational lower bounds under the enrichment model.
- `log10_bf`: log 10 enrichment Bayes factor.

Since this dataset is simulated from the enrichment model,
RSS yields a large log 10 enrichment Bayes factor as expected:

```matlab
>> log10_bf

log10_bf =

   19.3335
```
