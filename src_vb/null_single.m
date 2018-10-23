function [logw, alpha, mu, s] = null_single(method, file, Nsnp, h, log10odds, alpha, mu)
% USAGE: perform the variational inference for the RSS-BVSR baseline model,
%        where the values of hyper-parameters (h and log10odds) are fixed
% INPUT: 
%       method: the implementation of rss-varbvsr, string
%       file: the mat file that contains cell arrays of betahat, se and SiRiS, string
%       Nsnp: total sample size of each genetic variant, p by 1
%       h: proportion of phenotypic variance explained by available genotypes, scalar
%       log10odds: logarithm (base 10) of the prior odds for inclusion, scalar
%       alpha: initial values of the posterior inclusion probabilities, p by 1
%       mu: initial values of the expected additive effects (if the SNP is included), p by 1
% OUTPUT:
%       logw: unnormalized log-importance weight for the given hyper-parameters
%       alpha: variational estimates of the posterior inclusion probabilities, p by 1
%       mu: posterior means of the additive effects (if the SNP is included), p by 1
%       s: posterior variances of the additive effects (if the SNP is included), p by 1

  % make sure that both hyper-parameters are scalar
  if ~isscalar(h)
    error('Hyper-parameter h must be scalar for null_single.m.');
  end
  if ~isscalar(log10odds)
    error('Hyper-parameter log10odds must be scalar for null_single.m.');
  end

  % load summary statistics from the mat file
  sumstat = matfile(file);
  se      = sumstat.se;

  % pre-compute a quantity that is used to induce sigma_beta from h 
  se_vec  = cell2mat(se);
  xxyysum = sum(1 ./ ( Nsnp(:) .* (se_vec(:).^2) ));
  clear se_vec;

  % get the prior variance of the additive effects (sigma_beta^2)
  sigb2 = calc_beta_variance(xxyysum, sigmoid10(log10odds), h);

  % fit the RSS-BVSR model via variational approximation
  sigb    = sqrt(sigb2);
  logodds = log(10)*log10odds;
  options = struct('alpha', alpha, 'mu', mu);

  fprintf('Hyperparameters: h = %0.3f, log10odds = %+0.2f (sd = %0.3f)\n',h,log10odds,sigb);

  % NB: Compute the unnormalized log-importance weight given values for the
  % hyperparameters (h and log10odds). Implicitly, the importance weight
  % includes these terms: the likelihood, the prior, and the proposal. The
  % proposal and prior cancel out from the expression for the importance
  % weight because both are assumed to be uniform for all the hyperparameters.
  [logw,alpha,mu,s] = rss_varbvsr_bigmem_wrapper(method,file,sigb,logodds,options);

  fprintf('\n');

end

function psi = calc_beta_variance(xxyysum, pival, h)
% USAGE: calculate sigma_beta^2 using summary data, pi and pve

  psi = h ./ (pival .* xxyysum);
end

function y = sigmoid10(x)
% USAGE: compute the inverse of logit10(x)

  y = 1./(1 + 10.^(-x));
end

