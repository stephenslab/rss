## Example 5: Enrichment analysis of GWAS summary statistics

### Overview

This example illustrates how to perform enrichment analysis of GWAS summary statistics based on variational Bayes (VB) inference of RSS-BVSR models.

### Details

Here the enrichment analysis consists of fitting the following two models:

- **baseline model**: SNPs across the genome are equally likely to be associated with a target phenotype. 
- **enrichment model**: SNPs "near" a gene in a gene set are more likely (i.e. "enriched") to be associated with a target phenotype than remaining SNPs. 

For more details about these two models, please see [Carbonetto and Stephens (2013)](http://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1003770).

The key difference between our work and previous work [notably, Carbonetto and Stephens (2013)] is that we use **publicly available** summary data, rather than individual-level genetic data, in our enrichment analyses.

**To reproduce results of Example 5, please use scripts in the directory [example5](https://github.com/stephenslab/rss/tree/master/examples/example5), and follow the step-by-step instruction below.**

Before running scripts of this example, please make sure the [VB subroutines](https://github.com/stephenslab/rss/tree/master/src_vb) of RSS are installed. See instructions [here](RSS-via-VB).

### Step-by-step illustration

#### Fitting the baseline model ([`example5_null.m`](https://github.com/stephenslab/rss/blob/master/examples/example5/example5_null.m))

**Reminder**: Since we have to fit multiple regression analyses of 1.1 million common SNPs multiple times when fitting the baseline model, we need to use [rss_varbvsr_bigmem_squarem.m](https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem_squarem.m), which requires [Parallel Computing Toolbox](https://www.mathworks.com/help/distcomp/index.html). If you do not have this toolbox available, please skip this section, and use our result files for enrichment analyses instead (see Step 4). 

**Step 1**. [Download](http://projects.rcc.uchicago.edu/mstephens/data.txt) the input data file `ibd2015_sumstat.mat`, which contains the GWAS summary statistics and LD matrix estimates. Please contact us if you have trouble accessing this file.

Before proceeding to next step, let's look at the contents of `ibd2015_sumstat.mat`.
```matlab
>> sumstat = matfile('ibd2015_sumstat.mat');
>> sumstat

sumstat = 

  matlab.io.MatFile

  Properties:
      Properties.Source: '/project/mstephens/public_html/ibd2015_sumstat.mat'
    Properties.Writable: false                                               
                  SiRiS: [22x1 cell]                                         
                betahat: [22x1 cell]                                         
                    chr: [22x1 cell]                                         
                    pos: [22x1 cell]                                         
                     se: [22x1 cell]
```

Note that GWAS summary statistics and LD estimates are saved as [cell arrays](https://www.mathworks.com/help/matlab/cell-arrays.html). For each Chromosome `j`,

- `betahat{j, 1}` stores single-SNP effect size estimates of all SNPs in Chromosome `j`;
- `se{j, 1}` stores the corresponding standard errors of `betahat{j, 1}`;
- `chr{j, 1}` and `pos{j, 1}` store the physical positions of these SNPs;
- `SiRiS{j, 1}` stores a [sparse](https://www.mathworks.com/help/matlab/ref/sparse.html) matrix, defined as `repmat((1./se), 1, p) .* R .* repmat((1./se)', p, 1)`, where `R` is the estimated LD matrix of these `p` SNPs.

**Step 2**. Specify several dataset-specific variables that are required by [`null_template.m`](https://github.com/stephenslab/rss/blob/master/src_vb/null_template.m), the template script of fitting a baseline model.
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

% specify input and output paths
input_path  = '/project/mstephens/public_html/';                   
output_path = './';
```

Note that you may want to change `input_path` into other directories (e.g. `./`) where you save the input data file `ibd2015_sumstat.mat`.

**Step 3**. Fit baseline models across a grid of hyper-parameters in parallel. We write a simple `sbatch` script ([`example5_null.sbatch`](https://github.com/stephenslab/rss/blob/master/examples/example5/example5_null.sbatch)), and submit it to a cluster.
```
sbatch example5_null.sbatch
```

The script [`example5_null.sbatch`](https://github.com/stephenslab/rss/blob/master/examples/example5/example5_null.sbatch) will make three copies of `example5_null.m` (because there are three sets of `(h, theta0)` in the grid), and assigns one copy to each `(h, theta0)` setting.

After the submission, these three jobs will run in three nodes simultaneously.
```
xiangzhu@midway-login1: squeue -u xiangzhu
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
        29451148_3    sandyb  ibd2015 xiangzhu  R       1:29      1 midway249
        29451148_1    sandyb  ibd2015 xiangzhu  R       2:42      1 midway426
        29451148_2    sandyb  ibd2015 xiangzhu  R       2:42      1 midway427
```

Note that the script [`example5_null.sbatch`](https://github.com/stephenslab/rss/blob/master/examples/example5/example5_null.sbatch) is designed for the computing clusters at the [University of Chicago Research Computing Center (RCC)](https://rcc.uchicago.edu/). You may want to modify this script if you run it on other cluster environments.

**Step 4**. Aggregate the baseline results for enrichment analyses. Due to the "competitive" nature of our statistical method, the enrichment analyses are built on the basis of baseline modelling results. As soon as the three jobs in Step 3 are completed, the baseline results are saved as `mat` files ([version 7.3](https://www.mathworks.com/help/matlab/import_export/mat-file-versions.html)).

In case you cannot fit the baseline models in your own computing environment, we also share the baseline result files: `http://projects.rcc.uchicago.edu/mstephens/[file_name]` (listed below).

```
/project/mstephens/public_html/ibd2015_null_h_30_theta0_285_seed_459_squarem_step1.mat
/project/mstephens/public_html/ibd2015_null_h_30_theta0_288_seed_459_squarem_step1.mat
/project/mstephens/public_html/ibd2015_null_h_30_theta0_290_seed_459_squarem_step1.mat
```

#### Fitting the enrichment model ([`example5_gsea.m`](https://github.com/stephenslab/rss/blob/master/examples/example5/example5_gsea.m))

Unlike the baseline model, fitting the enrichment model does not require any specialized toolbox. With the baseline result files in place (see Step 4 above), everyone should be able to run this part by typing the following command in a Matlab console:
```matlab
>> run example5_gsea.m
```

Since we utilize additional pathway information in the enrichment model, we need to add a new data file `ibd2015_sumstat_path2641.mat` here. Let's look at the content of this file first.

```matlab
>> sumstat = matfile('ibd2015_sumstat_path2641.mat');
>> sumstat

sumstat = 

  matlab.io.MatFile

  Properties:
      Properties.Source: '/project/mstephens/public_html/ibd2015_sumstat_path2641.mat'
    Properties.Writable: false                                                        
                  SiRiS: [2990x2990 double]                                           
                betahat: [2990x1    double]                                           
                    chr: [2990x1    double]                                           
                    pos: [2990x1    double]                                           
                     se: [2990x1    double]                                           
                   snps: [2990x1    double]
```

This file contains GWAS summary statistics and LD estimates of SNPs that are "annotated" by a biological pathway (*IL23-mediated signaling events*, Pathway Commons 2, PID, 37 genes). Here we annotate a SNP as "inside pathway" if this SNP is within 100 kb of transcribed region of any member gene in the pathway.

As in fitting the baseline model, you may want to change `allstat_path` and `sumstat_path` into other directories (e.g. `./`) where you save the input data files `ibd2015_sumstat.mat` and `ibd2015_sumstat_path2641.mat`.

Finally, let's look at results of fitting the enrichment model (http://projects.rcc.uchicago.edu/mstephens/ibd2015_gsea_seed_459_path2641_squarem.mat).

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

There are four types of variables in the result file:

- user-specified quantities in `example5_gsea.m`: `method` (implementation type), `myseed` (seed of random number generator), `h` & `theta0` & `theta` (grid values of hyper-parameters), and `snps` (indices of pathway-annotated SNPs);
- estimated variational parameters: `alpha` & `mu` & `s`, where `alpha(:, i, j)` & `mu(:, i, j)` & `s(:, i, j)` correspond to the hyper-parameter setting `[h, theta0(i), theta(j)]`; **note that these parameter estimates can be directly used to priortize genes within an enriched pathway; see our [paper](http://www.biorxiv.org/content/early/2017/07/08/160770) for more details**.
- estimated variational lower bounds and Bayes factor: `logw0` (lower bounds under the baseline model), `logw1` (lower bounds under the enrichment model), and `log10bf` (log 10 Bayes factor);
- computational time: `runtime` (unit: seconds).

### More examples

Our enrichment analyses of 31 complex traits and 4,026 gene sets are essentially “replications” of the example above. Our full results are publicly available [online](http://xiangzhu.github.io/rss-gsea/results/).  
