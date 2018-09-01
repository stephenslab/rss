[Zhu and Stephens (*Ann. Appl. Stat.*, 2017)]: https://projecteuclid.org/euclid.aoas/1507168840
[`example1.m`]: https://github.com/stephenslab/rss/blob/master/examples/example1.m

# Example 1: Fit RSS models via MCMC, and estimate SNP heritability (PVE).

## Overview

This example illustrates how to fit RSS models using MCMC
algorithms. Three types of prior distributions are considered:
BVSR in [Guan and Stephens (*Ann. Appl. Stat.*, 2011)](https://projecteuclid.org/euclid.aoas/1318514285),
BSLMM in [Zhou, Carbonetto and Stephens (*PLoS Genet.*, 2013)](https://doi.org/10.1371/journal.pgen.1003264),
and ASH in [Stephens (*Biostatistics*, 2017)](https://doi.org/10.1093/biostatistics/kxw041).
The MCMC output is further used to estimate the SNP heritability.
This example is closely related to Section 4.2 of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][]. 

## Details 
The single-SNP summary-level data are computed from a simulated GWAS dataset.
The GWAS data are simulated under the Scenario 2.1 in [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][].
Specifically, 100 "causal" SNPs are randomly drawn from 12758 SNPs on chromosome 16,
with effect sizes coming from standard normal N(0,1).
Effect sizes of remaining SNPs are zero. The true PVE (SNP heritability) is 0.2.

To reproduce results of Example 1, please read the step-by-step guide below and run [`example1.m`][].
Before running [`example1.m`][], please make sure the
[MCMC subroutines](https://github.com/stephenslab/rss/tree/master/src) are installed.
Please find installation instructions [here](RSS-via-MCMC).

## Step-by-step illustration

**Step 1**.
Download the simulated single-SNP summary-level data
[`example1.mat`](https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example1/).
Please contact me if you have trouble accessing this file.

The data file `example1.mat` contains the following elements.

- `betahat`: 12758 by 1 vector, single-SNP effect size estimate for each SNP
- `se`: 12758 by 1 vector, standard errors of the single-SNP effect size estimates
- `Nsnp`: 12758 by 1 vector, sample size of each SNP
- `R`: 12758 by 12758 matrix, LD matrix estimated from a reference panel
- `bwd`: integer, bandwidth of the banded matrix `R`
- `BR`: (`bwd`+1) by 12758 matrix, banded storage of matrix `R`   

Note that only `betahat`, `se`, `Nsnp` and `R` are needed for running MCMC.
The other two quantities, `bwd` and `BR`, are only used in SNP heritability calculation.

**Step 2**. Check the "small effects" model assumption.

Using single-SNP summary data, we can compute squared sample correlation between phenotype and each SNP,
and then check the "small effect" assumption by looking at these marginal squared correlation values
(please see Table 1 of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][] for more details).
The following code illustrates the "small effect" check in [`example1.m`][]<sup>1</sup>.

```matlab            
>> chatsqr = (betahat(:).^2) ./ (Nsnp(:).*(se(:).^2) + betahat(:).^2);
>> disp(prctile(log10(chatsqr), 0:25:100));
  -11.6029   -4.1154   -3.4721   -2.9962   -1.5982
```

Since our data are generated from genotypes of a single chromosome,
the simulated effect sizes per SNP are larger than would be expected in a typical GWAS
(Table 1 [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][]).

**Step 3**. Fit RSS-BVSR, RSS-BSLMM and RSS-ASH models via MCMC.

To fit RSS-BVSR and RSS-BSLMM models, we only need to specify the length of Markov chains for RSS software.

```matlab
Ndraw = 2e6;
Nburn = 2e5;
Nthin = 9e1;
[betasam, gammasam, hsam, logpisam, Naccept] = rss_bvsr(betahat, se, R, Nsnp, Ndraw, Nburn, Nthin);
[bsam, zsam, lpsam, hsam, rsam, Naccept] = rss_bslmm(betahat, se, R, Nsnp, Ndraw, Nburn, Nthin);
```

To fit RSS-ASH model, we need to specify the length of Markov chain
and a grid for the prior standard deviations of multiple-SNP effect sizes.

```matlab
Ndraw = 5e7;
Nburn = 1e7;
Nthin = 1e3;
sigma_beta = [0 0.001 0.003 0.01 0.03 0.1 0.3 1 3];
[bsam, zsam, wsam, lsam, Naccept] = rss_ash(betahat, se, R, Nsnp, sigma_beta, Ndraw, Nburn, Nthin);
```

**Step 4**. Estimate SNP heritability.

We use the posterior sample of multiple-SNP effect sizes (`bsam`)
to obtain the posterior sample of SNP heritability (`pvesam`).

```matlab
M = length(hsam); % the length of posterior simulations
pvesam = zeros(M,1); % preallocate the pve posterior estimates
for i = 1:M 
  pvesam(i) = compute_pve(bsam(i,:), betahat, se, Nsnp, bwd, BR, 1);
end
```

Recall that the SNP heritability estimator (Equation 3.8 of
[Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][]) involves vector-matrix-vector product.
To speed calculation, we exploit the banded structure of `R` and use the banded version
of vector-matrix-vector product implemented in `lapack`.
Hence, the banded storage `BR`, instead of the original form `R`, is used to calculate SNP heritability.

**Step 5**. Summarize results.

The dataset is simulated with the true SNP heritability (PVE) being 0.2.
The following table summarizes the posterior estimates (with 95% credible interval)
and the total computational time (including MCMC iterations and PVE calculations) for three models.

| Model     | PVE estimation       | Total time |
|-----------|----------------------|------------|
| RSS-BVSR  | 0.200 [0.125, 0.290] | 1.38 hours |
| RSS-BSLMM | 0.216 [0.136, 0.306] | 2.52 hours |
| RSS-ASH   | 0.197 [0.114, 0.286] | 6.69 hours |

The following histograms depict the posterior distributions of estimated SNP heritability under these three models.
![example1_pve](images/rss_example1_pve.png)

## More simulations

Simulations in Section 4.2 of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][]
are essentially "replications" of the example above.
To facilitate reproducible research, we make the simulated datasets
for Scenarios 2.1 and 2.2 in Section 4.2 of [Zhu and Stephens (*Ann. Appl. Stat.*, 2017)][]
available ([`rss_example1_simulations.tar.gz`](https://projects.rcc.uchicago.edu/mstephens/rss_wiki/example1/readme))<sup>2</sup>.

Each simulated dataset contains three files: `genotype.txt`, `phenotype.txt` and `simulated_data.mat`.
The files `genotype.txt` and `phenotype.txt` are the genotype and phenotype files for
[`GEMMA`](https://github.com/xiangzhou/GEMMA).
The file `simulated_data.mat` contains three cells.

```matlab
true_para = {pve, beta, gamma, sigma};
individual_data = {y, X};
summary_data = {betahat, se, Nsnp};
```

Only the `summary_data` cell above is used as the input for RSS methods.

The RSS methods also require an estimated LD matrix `R`.
This matrix `R` is provided in the file `genotype.mat`.

After applying RSS methods to these simulated data,
we obtain the following PVE estimation results.   

### Scenario 2.1 (sparse), True PVE = 0.2

<img src="images/pve2sparse_1.png" width="600">

### Scenario 2.1 (sparse), True PVE = 0.6

<img src="images/pve2sparse_2.png" width="600">

### Scenario 2.2 (polygenic), True PVE = 0.2

<img src="images/pve2polygenic_1.png" width="600">

### Scenario 2.2 (polygenic), True PVE = 0.6 

<img src="images/pve2polygenic_2.png" width="600">

--------

**Footnotes:**

1. The function [`prctile`](http://www.mathworks.com/help/stats/prctile.html)
used in [`example1.m`][] requires the [Statistics and Machine Learning Toolbox](http://www.mathworks.com/help/stats/index.html).
Please see this [commit](https://github.com/stephenslab/rss/pull/3/commits/566e149ed840a913bfef9c0d7bf82feb41d6735d)
(courtesy of Dr. [Dr. Peter Carbonetto](https://pcarbo.github.io/))
if this required toolbox is not available in your environment.

2. Currently these files are locked, since they contain individual-level genotypes
from Wellcome Trust Case Control Consortium (WTCCC, https://www.wtccc.org.uk/).
You need to get permission from WTCCC before we can share these files with you.  
