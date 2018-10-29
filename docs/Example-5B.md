[Zhu and Stephens (2018)]: https://www.nature.com/articles/s41467-018-06805-x
[rss_varbvsr_bigmem_squarem.m]: https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem_squarem.m
[example5]: https://github.com/stephenslab/rss/tree/master/examples/example5
[example5_null.m]: https://github.com/stephenslab/rss/blob/master/examples/example5/example5_null.m
[example5_null.sbatch]: https://github.com/stephenslab/rss/blob/master/examples/example5/example5_null.sbatch
[example5_gsea.m]: https://github.com/stephenslab/rss/blob/master/examples/example5/example5_gsea.m
[example5_gene.m]: https://github.com/stephenslab/rss/blob/master/examples/example5/example5_gene.m
[ibd2015_sumstat.mat]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/ibd2015_sumstat.mat
[ibd2015_sumstat_path2641.mat]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/ibd2015_sumstat_path2641.mat
[ibd2015_gsea_seed_459_path2641_squarem.mat]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/ibd2015_gsea_seed_459_path2641_squarem.mat
[Supplementary Tables 6-7]: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-06805-x/MediaObjects/41467_2018_6805_MOESM1_ESM.pdf
[ibd2015_path2641_genes.mat]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/ibd2015_path2641_genes.mat 
[ibd2015_path2641_genes_results.mat]: https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/ibd2015_path2641_genes_results.mat

# Example 5: Enrichment and prioritization analysis of GWAS summary statistics using RSS (Part B).

## Overview

This is Part B of [Example 5](Example-5),
which illustrates how to perform enrichment analysis of
GWAS summary statistics based on variational Bayes (VB) inference of RSS-BVSR model.
This part describes an end-to-end enrichment analysis of inflammatory bowel disease (IBD)
GWAS summary statistics ([Liu et al, 2015](https://www.ncbi.nlm.nih.gov/pubmed/26192919)) and
a gene set named *IL23-mediated signaling events* (Pathway Commons 2, PID, 37 genes) using RSS.
This part illustrates the actual data analyses performed in [Zhu and Stephens (2018)][].

To reproduce results of Example 5 Part B,
please use scripts in the directory [example5][],
and follow the step-by-step guide below.
Before running any script in [example5][], please install the
[VB subroutines](https://github.com/stephenslab/rss/tree/master/src_vb) of RSS.
Please find installation instructions [here](RSS-via-VB).

Since a genome-wide enrichment analysis is conducted here,
this part is more complicated than [Example 5 Part A](Example-5A),
which is based on a relatively small simulated dataset.
It is advisable to go through [Example 5 Part A](Example-5A) before
diving into this real data example.

Note that the working directory here is assumed to be `rss/examples/example5`.
Please modify certain scripts accordingly if a different directory is used.

## Step-by-step illustration

### Fitting the baseline model with [example5_null.m][]

Since we have to perform Bayesian multiple regression analyses of
1.1 million common SNPs multiple times when fitting the baseline model,
we need to use the memory-efficient and parallel implementation of RSS
[rss_varbvsr_bigmem_squarem.m][],
which requires [MATLAB Parallel Computing Toolbox](https://www.mathworks.com/help/distcomp/index.html).
If you do not have this toolbox available, please skip this section,
and use our result files for the following enrichment analysis instead (see Step 4). 

**Step 1**. Download the input data file [ibd2015_sumstat.mat][],
which contains the GWAS summary statistics and LD matrix estimates.
This file is large (17G) because it has a LD matrix of 1.1 million common SNPs.
Please contact me (`xiangzhu[at]uchicago.edu`) if you have trouble accessing this file.

Before proceeding to the next step, let's look at the contents of [ibd2015_sumstat.mat][].

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

Here GWAS summary statistics and LD estimates are stored as
[cell arrays](https://www.mathworks.com/help/matlab/cell-arrays.html).
For each Chromosome `j`,

- `betahat{j,1}` stores single-SNP effect size estimates of all SNPs on Chromosome `j`;
- `se{j,1}` stores standard errors of `betahat{j, 1}`;
- `chr{j,1}` and `pos{j, 1}` store physical positions of these SNPs (GRCh37 build);
- `SiRiS{j,1}` stores a sparse matrix,
defined as `repmat((1./se),1,p) .* R .* repmat((1./se)',p,1)`,
where `R` is the estimated LD matrix of these `p` SNPs.

**Step 2**. Specify several dataset-specific variables that are required by
[null_template.m](https://github.com/stephenslab/rss/blob/master/src_vb/null_template.m),
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
In practice we are using much larger grid of hyper-parameters (because we seldomly
have a sensible guess for these parameters and we should learn them from data);
please see [Supplementary Tables 6-7][] of [Zhu and Stephens (2018)][] for more details. 

**Step 3**. Fit baseline models across a grid of hyper-parameters in parallel.
We write a simple `sbatch` script, [example5_null.sbatch][], and submit it to a cluster
where [Slurm](https://slurm.schedmd.com/) has been installed.

```bash
sbatch example5_null.sbatch
```

The script [example5_null.sbatch][] makes three copies of [example5_null.m][]
(because there are three sets of `(h, theta0)` in the grid specified in Step 2),
and then assigns one copy to each `(h, theta0)` setting.

After the `sbatch` submission, these three jobs run in three different nodes simultaneously.

```bash
xiangzhu@midway-login1: squeue -u xiangzhu
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
        29451148_3    sandyb  ibd2015 xiangzhu  R       1:29      1 midway249
        29451148_1    sandyb  ibd2015 xiangzhu  R       2:42      1 midway426
        29451148_2    sandyb  ibd2015 xiangzhu  R       2:42      1 midway427
```

Note that [example5_null.sbatch][] is designed for the computing clusters at
[University of Chicago Research Computing Center](https://rcc.uchicago.edu/).
You may need to modify this script if you plan to run it on a different environment.

**Step 4**. Aggregate the baseline results for enrichment analyses.
As soon as the three jobs in Step 3 are completed, the baseline results are saved as `mat` files
([version 7.3](https://www.mathworks.com/help/matlab/import_export/mat-file-versions.html)).

In case you cannot fit the baseline models in your own computing environment,
you can download the following baseline result files from
<https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example5/>.

```
xiangzhu@midway-login1: ls ibd2015_null_*.mat
ibd2015_null_h_30_theta0_285_seed_459_squarem_step1.mat
ibd2015_null_h_30_theta0_288_seed_459_squarem_step1.mat
ibd2015_null_h_30_theta0_290_seed_459_squarem_step1.mat
```

### Fitting the enrichment model with [example5_gsea.m][]

Unlike fitting the baseline model,
fitting the enrichment model does not require any specialized toolbox
or any high-performance computing cluster.
With the baseline result files listed above in place
(i.e. downloaded and saved in the working directory `rss/examples/example5`),
everyone should be able to run this part in a modern desktop
by simply typing the following line in a Matlab console:

```matlab
>> run example5_gsea.m
```

Since we use additional pathway information when fitting the enrichment model,
here we need to download and add a new data file [ibd2015_sumstat_path2641.mat][]
to the working directory `rss/examples/example5`.
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
As in [Zhu and Stephens (2018)][], here we annotate a SNP
as "inside a pathway" if this SNP is within 100 kb of transcribed
region of any member gene in this pathway.

The results of fitting the enrichment model are saved as
[ibd2015_gsea_seed_459_path2641_squarem.mat][].

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

The first group consists of user-specified quantities in [example5_gsea.m][]:
`method` (model fitting algorithm), `myseed` (seed of random number generator),
`h` & `theta0` & `theta` (grid of hyper-parameters),
and `snps` (indices of "inside-pathway" SNPs in the genome).

The second group consists of estimated variational parameters
`alpha` & `mu` & `s` under the enrichment model,
where `alpha(:,i,j)` & `mu(:,i,j)` & `s(:,i,j)` correspond to
estimation under the hyper-parameter setting `[h, theta0(i), theta(j)]`.
Note that these variational parameter estimates can be further
used to prioritize genes within an enriched pathway; please see the next section.

The third group consists of estimated variational lower bounds under the
baseline (`logw0`) and enrichment (`logw1`) model,
and the log 10 enrichment Bayes factor (`log10bf`).

The log 10 enrichment Bayes factor is 22.1355,
suggesting strong enrichment of genetic associations within *IL23-mediated signaling events* for IBD. 
This seems consistent with the important role of IL23 in immune response. 

### Prioritizing pathway genes with [example5_gene.m][]

Now we use the baseline and enrichment model fitting results above
to prioritize 37 member genes within the target pathway.
Similar to [Example 5 Part A](Example-5A), we need to download and
save a gene information file [ibd2015_path2641_genes.mat][] in the
working directory `rss/examples/example5` first.

```matlab
>> load ibd2015_path2641_genes.mat
>> whos
  Name                  Size                Bytes  Class     Attributes

  Aseg            1081481x37            160000304  double    sparse
  chr             1081481x1               4325924  int32
  gene_chr             37x1                   296  double
  gene_id              37x1                  4494  cell
  gene_start           37x1                   296  double
  gene_stop            37x1                   296  double
  pos             1081481x1               4325924  int32

>> whos Aseg
  Name            Size                Bytes  Class     Attributes

  Aseg      1081481x37            160000304  double    sparse

>> full(unique(Aseg(:)))'
ans =
     0     1
```

As shown above, there are 1,081,481 genome-wide SNPs and 37 pathway
member genes (both based on GRCh Build 37) in this simulated dataset.
Here we assign a SNP `j` to a gene `g` (i.e. `Aseg(j,g)==1`)
if and only if SNP `j` is within 100 kb of transcribed region
of gene `g` (i.e. `[gene_start-100e3 gene_stop+100e3]`).

With [ibd2015_path2641_genes.mat][] in place, we can then prioritize
pathway genes by simply typing the following line in a Matlab console:

```matlab
>> run example5_gene.m
```

The results of prioritizing pathway genes are saved as [ibd2015_path2641_genes_results.mat][].

```matlab
>> load ibd2015_path2641_genes_results.mat
>> whos
  Name             Size            Bytes  Class     Attributes

  b_p1            37x1               296  double
  b_p2            37x1               296  double
  e_p1            37x1               296  double
  e_p2            37x1               296  double
  gene_chr        37x1               296  double
  gene_id         37x1              4494  cell
  gene_start      37x1               296  double
  gene_stop       37x1               296  double
```

Similar to [Example 5 Part A](Example-5A),
`[b/e]_p_[1/2]` denotes the posterior probability of each locus
(gene with 100 kb window) containing at lease 1/2
trait-associated SNPs under the baseline/enrichment model.

We can easily load the gene-level results in R as follows:

```r
> suppressPackageStartupMessages(library(R.matlab))
> suppressPackageStartupMessages(library(dplyr))
> 
> res <- R.matlab::readMat("ibd2015_path2641_genes_results.mat")
>
> df <- data.frame(id=unlist(res$gene.id), chr=res$gene.chr, start=res$gene.start, stop=res$gene.stop, b_p1=res$b.p1, e_p1=res$e.p1, b_p2=res$b.p2, e_p2=res$e.p2)
> names(df) <- c("gene","chr","start","stop","b_p1","e_p1","b_p2","e_p2")
> df <- dplyr::arrange(df, -e_p1, -b_p1)
>
> sum(df$e_p1 >= df$b_p1)
[1] 35
> sum(df$e_p2 >= df$b_p2)
[1] 36
```

We can see that the enrichment model produces stronger gene-level signals
for almost all 37 pathway genes than the baseline model:

|gene    | chr|     start|      stop|   b_p1|   e_p1|   b_p2|   e_p2|
|:-------|---:|---------:|---------:|------:|------:|------:|------:|
|IL23R   |   1|  67632169|  67725662| 1.0000| 1.0000| 1.0000| 1.0000|
|IL12B   |   5| 158741791| 158757481| 1.0000| 1.0000| 0.9994| 1.0000|
|IL19    |   1| 206972215| 207016326| 1.0000| 1.0000| 0.4685| 0.9933|
|JAK2    |   9|   4985245|   5128183| 1.0000| 1.0000| 0.7649| 0.9995|
|IFNG    |  12|  68548550|  68553521| 1.0000| 1.0000| 0.2420| 0.9234|
|STAT3   |  17|  40465343|  40540513| 0.9999| 1.0000| 0.6887| 0.9944|
|CCL2    |  17|  32582296|  32584222| 0.9994| 1.0000| 0.0941| 0.8375|
|IL18RAP |   2| 103035254| 103069025| 0.9994| 1.0000| 0.3131| 0.9374|
|IL18R1  |   2| 102979097| 103015218| 0.9995| 1.0000| 0.3393| 0.9346|
|TYK2    |  19|  10461204|  10491248| 0.9629| 0.9999| 0.1699| 0.8993|
|STAT5A  |  17|  40439565|  40463961| 0.9998| 0.9997| 0.0618| 0.6792|
|IL2     |   4| 123372625| 123377650| 0.2659| 0.9989| 0.0370| 0.9074|
|STAT4   |   2| 191894306| 192015925| 0.7813| 0.9976| 0.1240| 0.8868|
|STAT1   |   2| 191833762| 191878976| 0.7642| 0.9953| 0.0735| 0.7793|
|IL6     |   7|  22766766|  22771621| 0.2801| 0.9776| 0.0364| 0.7868|
|NFKBIA  |  14|  35870716|  35873960| 0.2497| 0.9298| 0.0318| 0.7264|
|CXCL9   |   4|  76922623|  76928641| 0.4166| 0.9189| 0.0991| 0.7151|
|PIK3R1  |   5|  67511584|  67597649| 0.1645| 0.9125| 0.0136| 0.6976|
|PIK3CA  |   3| 178866311| 178952500| 0.1889| 0.9094| 0.0185| 0.6923|
|IL18    |  11| 112013976| 112034840| 0.1543| 0.8683| 0.0105| 0.5506|
|IL24    |   1| 207070788| 207077484| 0.1595| 0.8586| 0.0121| 0.5677|
|IL17F   |   6|  52101484|  52109298| 0.1160| 0.8515| 0.0069| 0.5720|
|IL17A   |   6|  52051185|  52055436| 0.1155| 0.8504| 0.0068| 0.5697|
|MPO     |  17|  56347217|  56358296| 0.1085| 0.8483| 0.0059| 0.5614|
|NFKB1   |   4| 103422486| 103538459| 0.1073| 0.8290| 0.0059| 0.5313|
|RELA    |  11|  65421067|  65430443| 0.1103| 0.8263| 0.0061| 0.5195|
|IL1B    |   2| 113587337| 113594356| 0.1193| 0.8263| 0.0071| 0.5224|
|IL12RB1 |  19|  18170371|  18197697| 0.1286| 0.8109| 0.0078| 0.4811|
|NOS2    |  17|  26083792|  26127555| 0.1178| 0.7861| 0.0069| 0.4544|
|ALOX12B |  17|   7975954|   7991021| 0.0773| 0.7716| 0.0030| 0.4345|
|CD3E    |  11| 118175295| 118186890| 0.0946| 0.7692| 0.0044| 0.4299|
|SOCS3   |  17|  76352859|  76356158| 0.1541| 0.7264| 0.0115| 0.3666|
|CXCL1   |   4|  74735109|  74736959| 0.0813| 0.7103| 0.0033| 0.3517|
|ITGA3   |  17|  48133340|  48167849| 0.0577| 0.6682| 0.0017| 0.3035|
|IL23A   |  12|  56732663|  56734194| 0.0370| 0.5185| 0.0007| 0.1672|
|CD4     |  12|   6898638|   6929976| 0.0337| 0.4747| 0.0006| 0.1366|
|TNF     |   6|  31543350|  31546112| 0.0000| 0.0000| 0.0000| 0.0000|

## More examples

The enrichment analyses of 31 complex traits and 4,026 gene sets
reported in [Zhu and Stephens (2018)][] are
essentially 124,806 “replications” of the example above
(with larger grids on hyper-parameters).
Our full analysis results are publicly available at
<https://xiangzhu.github.io/rss-gsea/>.  
