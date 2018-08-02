[Zhu and Stephens (*bioRxiv*, 2017)]: https://doi.org/10.1101/160770
[`rss_varbvsr_bigmem_squarem.m`]: https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem_squarem.m
[example5]: https://github.com/stephenslab/rss/tree/master/examples/example5
[`example5_null.m`]: https://github.com/stephenslab/rss/blob/master/examples/example5/example5_null.m
[`example5_null.sbatch`]: https://github.com/stephenslab/rss/blob/master/examples/example5/example5_null.sbatch
[`example5_gsea.m`]: https://github.com/stephenslab/rss/blob/master/examples/example5/example5_gsea.m
[`ibd2015_sumstat.mat`]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/ibd2015_sumstat.mat
[`ibd2015_sumstat_path2641.mat`]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/ibd2015_sumstat_path2641.mat
[`ibd2015_gsea_seed_459_path2641_squarem.mat`]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/ibd2015_gsea_seed_459_path2641_squarem.mat

# Example 5: Enrichment analysis of GWAS summary statistics using RSS (Part B).

## Overview

This is Part B of [Example 5](Example-5),
which illustrates how to perform enrichment analysis of
GWAS summary statistics based on variational Bayes (VB) inference of RSS-BVSR model.
This part describes an end-to-end enrichment analysis of inflammatory bowel disease GWAS summary statistics
([Liu et al, *Nat Genet.*, 2015](https://www.ncbi.nlm.nih.gov/pubmed/26192919)) and
a gene set named *IL23-mediated signaling events* (Pathway Commons 2, PID, 37 genes) using RSS.
This part illustrates the actual data analyses performed in [Zhu and Stephens (*bioRxiv*, 2017)][].

To reproduce results of Example 5 Part B,
please use scripts in the directory [example5][],
and follow the step-by-step guide below.
Before running any script in [example5][], please make sure the
[VB subroutines](https://github.com/stephenslab/rss/tree/master/src_vb) of RSS are installed.
Please find installation instructions [here](RSS-via-VB).

Since a genome-wide enrichment analysis is conducted here,
this part is more complicated than [Example 5 Part A](Example-5A).
It is advisable to go through [Example 5 Part A](Example-5A) before
diving into this real data example.

## Step-by-step illustration

### Fitting the baseline model with [`example5_null.m`][]

Since we have to perform multiple regression analyses of 1.1 million
common SNPs multiple times when fitting the baseline model,
we need to use the memory-efficient and parallel implementation [`rss_varbvsr_bigmem_squarem.m`][],
which requires [Parallel Computing Toolbox](https://www.mathworks.com/help/distcomp/index.html).
If you do not have this toolbox available, please skip this section,
and use our result files for the following enrichment analysis instead (see Step 4). 

**Step 1**. Download the input data file [`ibd2015_sumstat.mat`][],
which contains the GWAS summary statistics and LD matrix estimates.
This file is large (17G) because it has a LD matrix of 1.1 million common SNPs.
Please contact me if you have trouble accessing this file.

Before proceeding to the next step, let's look at the contents of `ibd2015_sumstat.mat`.

```matlab
>> sumstat = matfile('ibd2015_sumstat.mat');
>> sumstat

sumstat = 
  matlab.io.MatFile
  Properties:
      Properties.Source: 'ibd2015_sumstat.mat'
    Properties.Writable: false                                               
                  SiRiS: [22x1 cell]                                         
                betahat: [22x1 cell]                                         
                    chr: [22x1 cell]                                         
                    pos: [22x1 cell]                                         
                     se: [22x1 cell]
```

Here GWAS summary statistics and LD estimates are saved as
[cell arrays](https://www.mathworks.com/help/matlab/cell-arrays.html).
For each Chromosome `j`,

- `betahat{j,1}` stores single-SNP effect size estimates of all SNPs on Chromosome `j`;
- `se{j,1}` stores standard errors of `betahat{j, 1}`;
- `chr{j,1}` and `pos{j, 1}` store physical positions of these SNPs (GRCh37);
- `SiRiS{j,1}` stores a [sparse](https://www.mathworks.com/help/matlab/ref/sparse.html) matrix,
defined as `repmat((1./se),1,p) .* R .* repmat((1./se)',p,1)`,
where `R` is the estimated LD matrix of these `p` SNPs.

**Step 2**. Specify several dataset-specific variables that are required by
[`null_template.m`](https://github.com/stephenslab/rss/blob/master/src_vb/null_template.m),
a template script that fits a baseline model.

```matlab
% specify trait-specific information
trait_name  = 'ibd2015';
sample_size = (12882+21770); % cases: 12,882; controls: 21,770

% specify grid of hyper-parameters
h_rv      = 0.3;
theta0_rv = (-2.9:0.025:-2.85);

% specify stage of analysis
stage = 'step1';

% specify random start and algorithm 
myseed = 459;
method = 'squarem';
```

Note that here the grid of hyper-parameters is set minimum for illustration purpose.
In practice we are using much larger grid of hyper-parameters; please see
[Supplementary Tables 6-7](https://www.biorxiv.org/content/biorxiv/suppl/2018/07/16/160770.DC2/160770-2.pdf)
of [Zhu and Stephens (*bioRxiv*, 2017)][] for more details. 

**Step 3**. Fit baseline models across a grid of hyper-parameters in parallel.
We write a simple `sbatch` script, [`example5_null.sbatch`][], and submit it to a cluster
where [Slurm](https://slurm.schedmd.com/) has been installed.

```
sbatch example5_null.sbatch
```

The script [`example5_null.sbatch`][] makes three copies of [`example5_null.m`][]
(because there are three sets of `(h, theta0)` in the grid),
and then assigns one copy to each `(h, theta0)` setting.

After the `sbatch` submission, these three jobs run in three different nodes simultaneously.

```
xiangzhu@midway-login1: squeue -u xiangzhu
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
        29451148_3    sandyb  ibd2015 xiangzhu  R       1:29      1 midway249
        29451148_1    sandyb  ibd2015 xiangzhu  R       2:42      1 midway426
        29451148_2    sandyb  ibd2015 xiangzhu  R       2:42      1 midway427
```

Note that [`example5_null.sbatch`][] is designed for the computing clusters at
[University of Chicago Research Computing Center (RCC)](https://rcc.uchicago.edu/).
You may want to modify this script if you plan to run it on a different cluster environment.

**Step 4**. Aggregate the baseline results for enrichment analyses.
As soon as the three jobs in Step 3 are completed, the baseline results are saved as `mat` files
([version 7.3](https://www.mathworks.com/help/matlab/import_export/mat-file-versions.html)).

In case you cannot fit the baseline models in your own computing environment,
you can download the following baseline result files from this public
[page](https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/).

```
xiangzhu@midway-login1: ls ibd2015_null_*.mat
ibd2015_null_h_30_theta0_285_seed_459_squarem_step1.mat
ibd2015_null_h_30_theta0_288_seed_459_squarem_step1.mat
ibd2015_null_h_30_theta0_290_seed_459_squarem_step1.mat
```

### Fitting the enrichment model with [`example5_gsea.m`][]

Unlike fitting the baseline model,
fitting the enrichment model does not require any specialized toolbox.
With the baseline result files listed above in place,
everyone should be able to run this part by simply
typing the following line in a `Matlab` console:

```matlab
>> run example5_gsea.m
```

Since we use additional pathway information when fitting the enrichment model,
here we need to download and add a new data file [`ibd2015_sumstat_path2641.mat`][].
Let's look at the content of this file first.

```matlab
>> sumstat = matfile('ibd2015_sumstat_path2641.mat');
>> sumstat

sumstat = 
  matlab.io.MatFile
  Properties:
      Properties.Source: 'ibd2015_sumstat_path2641.mat'
    Properties.Writable: false                                                        
                  SiRiS: [2990x2990 double]                                           
                betahat: [2990x1    double]                                           
                    chr: [2990x1    double]                                           
                    pos: [2990x1    double]                                           
                     se: [2990x1    double]                                           
                   snps: [2990x1    double]
```

This new file contains GWAS summary statistics and LD estimates
of SNPs that are "inside" a biological pathway named
*IL23-mediated signaling events* (Pathway Commons 2, PID, 37 genes).
As in [Zhu and Stephens (*bioRxiv*, 2017)][], here we annotate a SNP
as "inside a pathway" if this SNP is within 100 kb of transcribed
region of any member gene in this pathway.

The results of fitting the enrichment model are saved as
[`ibd2015_gsea_seed_459_path2641_squarem.mat`][].

```matlab
>> load ibd2015_gsea_seed_459_path2641_squarem.mat
>> whos
  Name            Size                 Bytes  Class     Attributes

  alpha        2990x3x101            7247760  double              
  h               1x1                      8  double              
  log10bf         1x1                      8  double              
  logw0           3x1                     24  double              
  logw1           3x101                 2424  double              
  method          1x7                     14  char                
  mu           2990x3x101            7247760  double              
  myseed          1x1                      8  double              
  runtime         1x1                      8  double              
  s            2990x3x101            7247760  double              
  snps         2990x1                  23920  double              
  theta         101x1                    808  double              
  theta0          3x1                     24  double
```

There are three groups of output variables in the result file.

The first group consists of user-specified quantities in [`example5_gsea.m`][].

- `method`: model fitting algorithm.
- `myseed`: seed of random number generator.
- `h` & `theta0` & `theta`: grid values of hyper-parameters.
- `snps`: indices of "inside-pathway" SNPs.

The second group consists of estimated variational parameters
`alpha` & `mu` & `s` under the enrichment model,
where `alpha(:,i,j)` & `mu(:,i,j)` & `s(:,i,j)` correspond to
estimation under the hyper-parameter setting `[h, theta0(i), theta(j)]`.
Note that these variational parameter estimates can be further
used to **priortize genes within an enriched pathway**;
please see [Zhu and Stephens (*bioRxiv*, 2017)][] for more details.

The third group consists of estimated variational lower bounds and Bayes factor.

- `logw0`: variational lower bounds under baseline models.
- `logw1`: variational lower bounds under enrichment models.
- `log10bf`: log 10 enrichment Bayes factor.
- `runtime`: computational time (unit: seconds).

## More examples

The enrichment analyses of 31 complex traits and 4,026 gene sets
reported in [Zhu and Stephens (*bioRxiv*, 2017)][] are
essentially 124,806 “replications” of the example above.
Our full analysis results are publicly available
[online](http://xiangzhu.github.io/rss-gsea/).  
