function [logw1,alpha,mu,s] = gsea_single(method,file,se_all,Nsnp_all,h,theta0,theta,alpha,mu,logw0,alpha0,mu0)
% USAGE: perform the variational inference for the RSS-BVSR model under the enrichment hypothesis where
%	 the values of hyper-parameters are fixed (i.e. h, theta0 and theta); similar to null_single.m
% INPUT: 
%       method: the implementation of rss-varbvsr, character
%       file: the path of mat file that contains cell arrays of betahat, se and SiRiS, string
%       se_all: standard errors for all SNPs in the genome, numsnp by 1
%       Nsnp_all: sample sizes for all SNPs in the genome, numsnp by 1
%       h: the fixed proportion of phenotypic variance explained by available genotypes, scalar
%       theta0: the fixed (base 10) logarithm of the prior odds for inclusion, scalar
%	theta: the fixed log-fold (base 10) enrichment parameter, scalar
%       logw0: the log-importance weight for the given theta0 under null, scalar
%       alpha0: the variational posterior inclusion probabilities under null, p by 1
%       mu0: the variational posterior means of the additive effects (given snp included), p by 1
%       alpha: the initial values of the posterior inclusion probabilities under erichment, p by 1 
%       mu: the initial values of the expected additive effects (given snp included) under enrichment, p by 1
% OUTPUT:
%       logw1: the unnormalized log-importance weight for each combination of hyper-parameters
%       alpha: variational estimates of the posterior inclusion probabilities, p by 1
%       mu: posterior means of the additive effects (given snp included), p by 1
%       s: posterior variances of the additive effects (given snp included), p by 1

  % load the indices of assigned SNPs from the mat file
  sumstat = matfile(file);
  snps    = sumstat.snps;
  fprintf('There are %d SNPs assigned to a pre-defined gene set ...\n', length(snps));

  % get the number of SNPs in the whole genome (ng)
  ng = length(se_all);
  fprintf('There are %d SNPs on the whole genome ...\n', ng);

  % get the number of combinations of the hyperparameters (ns)
  % this function is only for the case where ns == 1
  ns = numel(h)*numel(theta0)*numel(theta);
  if ns ~= 1
    error('The hyperparameter grid exists; use gsea_wrapper.m or gsea_bigmem_wrapper.m instead ...');
  end

  % first compute the variational lower bound to the marginal log-likelihood
  % under the null hypothesis when there is no enrichment, only for SNPs in
  % the gene set; note that rss_varbvsr defines the log-odds ratio using the
  % natural logarithm, so we need to multiply theta0 by log(10)

  log10odds_null = theta0 * ones(ng,1); % genome-wide baseline  

  sigb0    = calc_beta_sd(se_all, Nsnp_all, sigmoid10(log10odds_null), h);
  logodds0 = log(10) * theta0;
  options  = struct('alpha',alpha0,'mu',mu0);

  % run rss-varbvsr under the null hypothesis
  % note: there is a duplicated calculation for the same (h, theta0)
  % but this issue does not exist for the wrapper versions (gsea_wrapper.m/gsea_bigmem_wrapper.m)
  F0 = rss_varbvsr_bigmem_wrapper(method,file,sigb0,logodds0,options);

  % next compute the marginal log-likelihood under the alternative
  % hypothesis only for SNPs in the pathway; note that rss_varbvsr
  % defines the log-odds ratio using the natural logarithm, so we
  % need to multiply (theta0+theta) by log(10)

  log10odds_gsea       = theta0 * ones(ng,1);  % genome-wide baseline
  log10odds_gsea(snps) = theta0 + theta;       % enriched in the gene set

  sigb1    = calc_beta_sd(se_all, Nsnp_all, sigmoid10(log10odds_gsea), h);
  logodds1 = log(10) * (theta0 + theta);
  options  = struct('alpha',alpha,'mu',mu);

  % run rss-varbvsr under the enrichment hypothesis
  [F1,alpha,mu,s] = rss_varbvsr_bigmem_wrapper(method,file,sigb1,logodds1,options);

  % compute the variational lower bound under the alternative; for an explanation
  % why we can decompose the variational lower bound in this way, see pp 5-6 of
  % Supplementary Text of Carbonetto and Stephens (PLoS Genetics, 2013)

  logw1 = logw0 + F1 - F0;

  fprintf('\n');

end

function sigb = calc_beta_sd(se, Nsnp, pival, h)
% USAGE: calculate sigma_beta using summary data, pi and pve under enrichment 

  pxy  = pival ./ (Nsnp .* (se.^2));
  psi  = h ./ sum(pxy);
  sigb = sqrt(psi);
end

function y = sigmoid10(x)
% USAGE: compute the inverse of logit10(x)

  y = 1./(1 + 10.^(-x));
end

