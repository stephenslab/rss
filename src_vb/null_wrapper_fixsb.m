function [logw, alpha, mu, s] = null_wrapper_fixsb(method, betahat, se, SiRiS, sigb, log10odds, alpha, mu)
% USAGE: perform the full variational inference for the RSS-BVSR baseline model
%        baseline model: the candidate gene set is not enriched for genotype-phenotype association
%        here the prior variance of causal effects is fixed as a scalar, sigb^2
% INPUT:
%       method: the implementation of rss-varbvsr, character 
%       betahat: effect size estimates under single-SNP model, p by 1 array or C by 1 cell array
%       se: standard errors of betahat, p by 1 array or C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), sparse matrix (CCS format), p by p array or C by 1 cell array
%       sigb: prior SD of the regression coefficients (if included), scalar
%       log10odds: logarithm (base 10) of the prior odds for inclusion, n0 by 1
%       alpha: initial values of the posterior inclusion probabilities, p by n0
%       mu: initial values of the expected additive effects (if the SNP is included), p by n0
% OUTPUT:
%       logw: the unnormalized log-importance weight for each combination of hyper-parameters, n0 by 1
%       alpha: the variational posterior inclusion probabilities, p by n0
%       mu: the variational posterior means of additive effects (if the SNP is included), p by n0
%       s: the variational posterior variances of additive effects (if the SNP is included), p by n0

% OVERVIEW:
% This inference procedure involves an inner loop and an outer loop. The
% inner loop consists of running a coordinate ascent algorithm to tighten
% the variational lower bound given a setting of the hyperparameters. The
% outer loop computes importance weights for all combinations of the
% hyperparameters.

  % make sure that prior SD sigb is scalar
  if ~isscalar(sigb)
    error('Hyper-parameter sigb must be scalar for null_wrapper_fixsb.m.');
  end

  % get the number of settings of the genome-wide log-odds (n0)
  n0 = numel(log10odds);

  % step 1: get the best initialization for the variational parameters
  fprintf('Finding best initialization for %d combinations ', n0);
  fprintf('of hyperparameters.\n');
  [logw,alpha,mu] = null_outerloop(method,betahat,se,SiRiS,sigb,alpha,mu,log10odds);

  % step 2: choose an initialization common to all the runs of the coordinate
  % ascent algorithm. This is chosen from the hyperparameters with the highest
  % variational estimate of the posterior probability
  [~, i] = max(logw(:));
  alpha  = repmat(alpha(:,i), [1 n0]);
  mu     = repmat(mu(:,i), [1 n0]);

  % step 3: compute the unnormalized log-importance weights
  fprintf('Computing importance weights for %d combinations ', n0);
  fprintf('of hyperparameters.\n');
  [logw,alpha,mu,s] = null_outerloop(method,betahat,se,SiRiS,sigb,alpha,mu,log10odds);

end

function [logw, alpha, mu, s] = null_outerloop(method, betahat, se, SiRiS, sigb, alpha, mu, log10odds)
% USAGE: the outer loop of the variational inference for the RSS-BVSR baseline model
% INPUT: 
%       method: the implementation of rss-varbvsr, character
%       betahat: effect size estimates under single-SNP model, p by 1 array or C by 1 cell array
%       se: standard errors of betahat, p by 1 array or C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), sparse matrix (CCS format), p by p array or C by 1 cell array
%       sigb: prior SD of the regression coefficients (if included), scalar
%       alpha: initial values of the posterior inclusion probabilities, p by n0
%       mu: initial values of the expected additive effects (if the SNP is included), p by n0
%       log10odds: logarithm (base 10) of the prior odds for inclusion, n0 by 1
% OUTPUT:
%       logw: the unnormalized log-importance weight for each combination of hyper-parameters, n0 by 1
%       alpha: variational estimates of the posterior inclusion probabilities, p by n0
%       mu: posterior means of the additive effects (if the SNP is included), p by n0
%       s: posterior variances of the additive effects (if the SNP is included), p by n0

  % get the total number of SNPs analyzed (p)
  if iscell(betahat)
    p = length(cell2mat(betahat));
  else
    p = length(betahat);
  end
 
  % get the number of settings of the genome-wide log-odds (n0)
  n0 = numel(log10odds);

  % initialize storage for the unnormalized log-importance weights, and
  % the variational variances of the additive effects
  logw = zeros(n0,1);
  s    = zeros(p,n0);

  % repeat for each combination of the hyperparameters
  for i = 1:n0

    logodds = log(10)*log10odds(i);
    options = struct('alpha',alpha(:,i),'mu',mu(:,i));

    % display the current status of the inference procedure
    fprintf('(%03d) log10odds = %+0.2f (sd = %0.3f)\n',i,log10odds(i),sigb);
  
    % NB: Compute the unnormalized log-importance weight given values for the
    % hyperparameters (h and log10odds). Implicitly, the importance weight
    % includes these terms: the likelihood, the prior, and the proposal. The
    % proposal and prior cancel out from the expression for the importance
    % weight because both are assumed to be uniform for all the hyperparameters.
    [logw(i),alpha(:,i),mu(:,i),s(:,i)] = rss_varbvsr_wrapper(method,betahat,se,SiRiS,sigb,logodds,options);

    fprintf('\n');

  end

end



