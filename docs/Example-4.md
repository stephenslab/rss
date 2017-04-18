### Example 4: Fit the RSS model via VB

#### Overview

This example illustrates how to fit an RSS-BVSR model using variational Bayes (VB) approximation. 

Based on the theoretical derivations, we further set the input values of `se` and `R` such that the VB analysis of summary-level data is equivalent to the VB analysis of individual-level data ([Carbonetto and Stephens, 2012](https://projecteuclid.org/euclid.ba/1339616726)). Hence, this example provides an *in silico* sanity check for our theoretical work. 

#### Details

The summary-level data are computed from a simulated GWAS dataset. The GWAS data are simulated by [`enrich_datamaker.m
`](https://github.com/stephenslab/rss/blob/master/misc/enrich_datamaker.m), which makes "signal-enriched" data based on the gene set. Specifically, SNPs outside the gene set are selected to be causal ones with probability `10^theta0`, whereas SNPs inside the gene set are selected with a higher probability `10^(theta0+theta)`, where `theta>0`. For more details, see [Carbonetto and Stephens (2013)](http://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1003770). 

After generating the data, we feed the individual-level and summary-level data to [`varbvs`](https://github.com/pcarbo/varbvs) and [`rss-varbvsr`](https://github.com/stephenslab/rss/tree/master/src_vb) respectively, and then compare the VB output from these two methods.

**To reproduce results of Example 4, please use [example4.m](https://github.com/stephenslab/rss/blob/master/examples/example4.m).**

Before running [example4.m](https://github.com/stephenslab/rss/blob/master/examples/example4.m), please make sure the [VB subroutines](https://github.com/stephenslab/rss/tree/master/src_vb) of RSS are installed. See instructions [here](https://github.com/stephenslab/rss/wiki/RSS-via-VB).

#### Step-by-step illustration

**Step 1**. [Download](https://uchicago.box.com/v/example4) the input data `genotype.mat` and `AH_chr16.mat` for `enrich_datamaker.m`. Please contact us if you have trouble accessing this file.

**Step 2**. Install the `MATLAB` implementation of [`varbvs`](https://github.com/pcarbo/varbvs). Please follow the instruction [here](https://github.com/pcarbo/varbvs/tree/master/varbvs-MATLAB#large-scale-bayesian-variable-selection-for-matlab). After the installation, add `varbvs` to the search path (see the example below).
```matlab
addpath('/home/xiangzhu/varbvs-master/varbvs-MATLAB/');
```

**Step 3**. Extract the SNPs that inside the gene set. This step is where we need the input data [`AH_chr16.mat`](https://uchicago.box.com/v/example4). The index of SNPs inside the gene set is stored as `snps`.
```matlab
AH    = matfile('AH_chr16.mat');
H     = AH.H; % hypothesis matrix 3323x3160
A     = AH.A; % annotation matrix 12758x3323
paths = find(H(:,end)); % signal transduction (Biosystem, Reactome)
snps  = find(sum(A(:,paths),2)); % index of variants inside the pathway
```

**Step 4**. Simulate the enriched dataset. In addition to `genotype.mat` and `AH_chr16.mat`, four more input data are required in order to run `example4.m`. You can provide them through keyboard.

```matlab
% set the number of replications
prompt = 'What is the number of replications? ';
Nrep   = input(prompt);

% set the genome-wide log-odds 
prompt1 = 'What is the genome-wide log-odds? ';
theta0  = input(prompt1);

% set the log-fold enrichment
prompt2 = 'What is the log-fold enrichment? ';
theta   = input(prompt2);

% set the true pve
prompt3 = 'What is the pve (between 0 and 1)? ';
pve = input(prompt3);
``` 

With this in place, the data generation step is as follows.
```matlab
for i = 1:Nrep
  myseed = 617 + i;

  % generate data under enrichment hypothesis
  [true_para,individual_data,summary_data] = enrich_datamaker(C,thetatype,pve,myseed,snps);
  
  ... ...

end
```

**Step 5**. Ensure that `varbvs` and `rss-varbvsr` run in an almost identical environment. We fix all the hyper-parameters as their true values, and give the same parameter initialization.
```matlab
  % fix all hyper-parameters as their true values
  sigma   = true_para{4}^2;       % the true residual variance
  sa      = 1/sigma;              % sa*sigma being the true prior variance of causal effects
  sigb    = 1;                    % the true prior SD of causal effects 
  theta0  = thetatype(1);         % the true genome-wide log-odds (base 10)
  theta   = thetatype(2);         % the true log-fold enrichment (base 10)

  % initialize the variational parameters
  rng(myseed, 'twister');
  alpha0 = rand(p,1);
  alpha0 = alpha0 ./ sum(alpha0);

  rng(myseed+1, 'twister');
  mu0 = randn(p,1);
```

**Step 6**. Run `varbvs` on the individual-level genotypes and phenotypes. The VB analysis involves two steps: first fit the model assuming no enrichment; second fit the model using the enrichment information. We then calculate a Bayes factor based on the marginal likelihoods of these two models, and use it to test whether the gene set is enriched for genetic associations.
```matlab
fit_null = varbvs(X,[],y,[],'gaussian',options_n);
fit_gsea = varbvs(X,[],y,[],'gaussian',options_e);

bf_b = exp(fit_gsea.logw - fit_null.logw);
```

**Step 7**. Run `rss-varbvsr` on the single-SNP association summary statistics. Here we do not specify `se` and `R` as we actually do in our real data analyses. Instead, we set their values such that `rss-varbvsr` and `varbvs` are expected to produce the same output. Based on our derivation, we find that `rss-varbvsr` is equivalent to `varbvs` if and only if

1. the variance of phenotype is the same as residual variance `sigma.^2`;
2. the input LD matrix `R` is the same as the sample correlation matrix of cohort genotypes `X`

Interestingly, these two assumptions are **exactly the same assumptions in Proposition 2.1 of RSS [paper](http://biorxiv.org/content/early/2016/03/04/042457), which guarantees that the RSS likelihood is equivalent to the Gaussian likelihood of individual-level data.**

These two assumptions are implemented as follows.
```matlab
  % set the summary-level data for a perfect match b/t rss-varbvsr and varbvs
  betahat = summary_data{1};
  se      = sqrt(sigma ./ diag(X'*X));          % condition 1 for perfect matching
  Si      = 1 ./ se(:);
  R       = corrcoef(X);                        % condition 2 for perfect matching
  SiRiS   = sparse(repmat(Si, 1, p) .* R .* repmat(Si', p, 1));
```

Now we perform the VB analysis of summary data via `rss-varbvsr`.
```matlab
  [logw_nr,alpha_nr,mu_nr,s_nr] = rss_varbvsr(betahat,se,SiRiS,sigb,logodds_n,options);
  [logw_er,alpha_er,mu_er,s_er] = rss_varbvsr(betahat,se,SiRiS,sigb,logodds_e,options);

  bf_r = exp(logw_er-logw_nr);
```
**Step 8**. Compare VB output from `varbvs` and `rss-varbvsr`. Both `varbvs` and `rss-varbvsr` output an estimated posterior distribution of `beta` (the multiple regression coefficients, or, the multiple-SNP effects). Specifically, the estimated distributions have the same analytical form (Equations 6-7 in [Carbonetto and Stephens (2012)](https://projecteuclid.org/euclid.ba/1339616726)).

![](https://github.com/xiangzhu/pubfig/blob/master/wiki/varbvs_output.png)

Hence, it suffices to compare the estimated `[alpha, mu, s]`, and the estimated Bayes factors based on the VB output. Here is the result when I only run `example4.m` with one replication.
```matlab
>> run example4.m                           
What is the number of replications? 1
What is the genome-wide log-odds? -4
What is the log-fold enrichment? 2
What is the pve (between 0 and 1)? 0.3
```
The log10 maximum absolute differences of `[alpha, mu, s]` between `varbvs` and `rss-varbvsr`, under the null model.
```matlab
>> disp(log10([alpha_n_diff mu_n_diff s_n_diff]))
   -4.7084   -4.5409   -7.6421
``` 
The log10 maximum absolute differences of `[alpha, mu, s]` between `varbvs` and `rss-varbvsr`, under the enrichment model.
```matlab
>> disp(log10([alpha_e_diff mu_e_diff s_e_diff]))
   -3.9865   -4.0157   -7.6421
```

The log10 Bayes factors estimated from `varbvs` and `rss-varbvsr` are 4.0193 and 4.0193 respectively, and the log10 relative difference between them is -6.8345.

### More simulations

To see how `example4.m` behave on average, we need run more replication. Below are the results from 100 replications.
```matlab
>> run example4.m                           
What is the number of replications? 100
What is the genome-wide log-odds? -4
What is the log-fold enrichment? 2
What is the pve (between 0 and 1)? 0.3
``` 

![](https://github.com/xiangzhu/pubfig/blob/master/wiki/rss_example4_rep100.png) 
