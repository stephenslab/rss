## Example 1: Fit the RSS model via MCMC

### Overview
This example illustrates how to fit an RSS model using MCMC simulation. Three types of prior distributions are considered: BVSR, BSLMM and ASH. The output of MCMC is further used to estimate SNP heritability.

### Details 
The summary-level data are computed from a simulated GWAS dataset. The GWAS data are simulated under the Scenario 2.1 in the RSS [paper](https://doi.org/10.1101/042457). Specifically, 100 "causal" SNPs are randomly drawn from 12758 SNPs on chromosome 16 (WTCCC UK Blood Service Control Group), with effect sizes coming from N(0,1). The true PVE (SNP heritability) is 0.2.

The population LD matrix is estimated from a reference panel (WTCCC 1958 British Birth Cohort), using the shrinkage estimator in [Wen and Stephens](http://stephenslab.uchicago.edu/assets/papers/Wen2010.pdf) (2010).

**To reproduce results of Example 1, please use [`example1.m`](https://github.com/stephenslab/rss/blob/master/examples/example1.m).** 

Before running [`example1.m`](https://github.com/stephenslab/rss/blob/master/examples/example1.m), please make sure the [MCMC subroutines](https://github.com/stephenslab/rss/tree/master/src) of RSS are installed. See instructions [here](RSS-via-MCMC).

### Step-by-step illustration

**Step 1**. Download the simulated summary-level data [`example1.mat`](https://uchicago.box.com/example1). Please contact us if you have trouble accessing this file.

The `example1.mat` contains the following data.
- `betahat`: 12758 by 1 vector, the single-SNP effect size estimate for each SNP
- `se`: 12758 by 1 vector, the standard errors of the single-SNP effect size estimates
- `Nsnp`: 12758 by 1 vector, the sample size of each SNP
- `R`: 12758 by 12758 matrix, the LD matrix estimated from a reference panel
- `bwd`: integer, the bandwidth of the matrix `R`
- `BR`: (`bwd`+1) by 12758 matrix, the banded storage of the matrix `R`   

Note that only `betahat`, `se`, `Nsnp` and `R` are needed for model fitting. The other two quantities, `bwd` and `BR`, are used in SNP heritability calculation.

**Step 2**. Check the "small effects" model assumption.

Using the summary data, we can compute the squared sample correlation between the phenotype and each SNP. We check the "small effects" assumption by looking at these marginal squared correlation values. (NB: The function [`prctile`](http://www.mathworks.com/help/stats/prctile.html) requires the [Statistics and Machine Learning Toolbox](http://www.mathworks.com/help/stats/index.html). Please see this [commit](https://github.com/stephenslab/rss/pull/3/commits/566e149ed840a913bfef9c0d7bf82feb41d6735d) (courtesy of Dr. [Dr. Peter Carbonetto](https://pcarbo.github.io/)) if the required toolbox is not available.)
```matlab            
>> chatsqr = (betahat(:).^2) ./ (Nsnp(:).*(se(:).^2) + betahat(:).^2);
>> disp(prctile(log10(chatsqr), 0:25:100));
  -11.6029   -4.1154   -3.4721   -2.9962   -1.5982
```
Since our data is simulated from genotypes on a single chromosome, the simulated effect sizes per SNP are larger than would be expected in a typical GWAS (Table 1 in the RSS [paper](https://doi.org/10.1101/042457)).

**Step 3**. Fit the RSS-BVSR, RSS-BSLMM and RSS-ASH model via MCMC.

To fit the RSS-BVSR and RSS-BSLMM model, only the length of MCMC simulation is needed.
```matlab
Ndraw = 2e6;
Nburn = 2e5;
Nthin = 9e1;
[betasam, gammasam, hsam, logpisam, Naccept] = rss_bvsr(betahat, se, R, Nsnp, Ndraw, Nburn, Nthin);
[bsam, zsam, lpsam, hsam, rsam, Naccept] = rss_bslmm(betahat, se, R, Nsnp, Ndraw, Nburn, Nthin);
```
The fitting of RSS-ASH model requires a grid for the prior standard deviations for the effect sizes.
```matlab
Ndraw = 5e7;
Nburn = 1e7;
Nthin = 1e3;
sigma_beta = [0 0.001 0.003 0.01 0.03 0.1 0.3 1 3];
[bsam, zsam, wsam, lsam, Naccept] = rss_ash(betahat, se, R, Nsnp, sigma_beta, Ndraw, Nburn, Nthin);
```

**Step 4**. Estimate the SNP heritability.

We now use the posterior of multiple-SNP effect sizes to estimate the SNP heritability.
```matlab
M = length(hsam); % the length of posterior simulations
pvesam = zeros(M,1); % preallocate the pve posterior estimates
for i = 1:M 
  pvesam(i) = compute_pve(bsam(i,:), betahat, se, Nsnp, bwd, BR, 1);
end
```
Recall that the SNP heritability estimator (Equation 3.8 in the RSS [paper](https://doi.org/10.1101/042457)) involves vector-matrix-vector product. To speed the calculation, we exploit the banded structure of `R` and use the banded version of vector-matrix-vector product implemented in `lapack`. Hence, the banded storage `BR`, instead of the original form `R`, is used to calculate the SNP heritability.

**Step 5**. Summarize the results.

The data is simulated with the true SNP heritability (PVE) being **0.2**. The following table summarizes the posterior estimates (with 95% credible interval) and the total computational time (including MCMC iterations and PVE calculations) for the three models.

| Model     | PVE estimation       | Total time |
|-----------|----------------------|------------|
| RSS-BVSR  | 0.200 [0.125, 0.290] | 1.38 hours |
| RSS-BSLMM | 0.216 [0.136, 0.306] | 2.52 hours |
| RSS-ASH   | 0.197 [0.114, 0.286] | 6.69 hours |

The following histograms depict the posterior distributions of estimated SNP heritability under the three models.
![example1_pve](images/rss_example1_pve.png)

### More simulations

The simulations in Section 4.2 of the RSS [paper](https://doi.org/10.1101/042457) are essentially "replications" of the example above. To facilitate reproducible research, we make the simulated datasets for Scenario 2.1 and 2.2 publicly available ([`rss_example1_simulations.tar.gz`](https://uchicago.box.com/example1)).

Each simulated dataset contains three files: `genotype.txt`, `phenotype.txt` and `simulated_data.mat`. The files `genotype.txt` and `phenotype.txt` are the genotype and phenotype files for [`GEMMA`](https://github.com/xiangzhou/GEMMA). The file `simulated_data.mat` contains three cells.
```matlab
true_para = {pve, beta, gamma, sigma};
individual_data = {y, X};
summary_data = {betahat, se, Nsnp};
```
Only the `summary_data` cell is used as the input for RSS methods.

The RSS methods also require an estimated LD matrix. This matrix `R` is provided in the file `genotype.mat`.

After applying RSS methods to these simulated data, we obtain the following results.   

**Scenario 2.1 (sparse), True PVE = 0.2**

<img src="images/pve2sparse_1.png" width="600">

**Scenario 2.1 (sparse), True PVE = 0.6**

<img src="images/pve2sparse_2.png" width="600">

**Scenario 2.2 (polygenic), True PVE = 0.2**

<img src="images/pve2polygenic_1.png" width="600">

**Scenario 2.2 (polygenic), True PVE = 0.6** 

<img src="images/pve2polygenic_2.png" width="600">
