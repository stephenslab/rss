function [log10bf,logw1,alpha,mu,s] = gsea_wrapper_fixsb(method,betahat,se,SiRiS,snps,sigb,theta0,theta,logw0,alpha0,mu0)
% USAGE: perform the full variational inference for the RSS-BVSR model under enrichment hypothesis
%        enrichment hypothesis: the gene set is enriched for genotype-phenotype association
%        here the prior variance of causal effects is fixed as sigb^2
% INPUT:
%       method: the implementation of rss-varbvsr, character 
%       betahat: the effect size estimates under single-SNP model, p by 1 array or C by 1 cell array
%       se: standard errors of betahat, p by 1 array or C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), double precision sparse matrix (ccs format), p by p array or C by 1 cell array
%	snps: the genomic indices of SNPs that are assigned to the pathway
%       sigb: the prior SD of the regression coefficients (if included), scalar
%       theta0: the grid of the genome-wide log-odds (base 10), n0 by 1
%       theta: the grid of the enrichment (base 10), n1 by 1
%	logw0: the log-importance weights for the grid of theta0 under null, n0 by 1
%       alpha0: the variational posterior inclusion probabilities under null, p by n0
%       mu0: the variational posterior means of the additive effects (given snp included), p by n0
% OUTPUT:
%	log10bf: log 10 of gene set enrichment BF based on summary-level data, scalar
%       logw1: the unnormalized log-importance weight for each combination of hyper-parameters, n0 by n1
%       alpha: the variational posterior inclusion probabilities, p by n0 by n1
%       mu: the variational posterior means of the additive effects (given snp included), p by n0 by n1
%       s: the variational posterior variances of the additive effects (given snp included), p by n0 by n1

% OVERVIEW:
% This inference procedure involves an inner loop and an outer loop. The
% inner loop consists of running a coordinate ascent algorithm to tighten
% the variational lower bound given a setting of the hyperparameters. The
% outer loop computes importance weights for all combinations of the
% hyperparameters.

  % Get the number of SNPs assigned to the enriched pathway (p), 
  % the number of settings of the genome-wide log-odds (n0),
  % and the number of settings of the enrichment parameter (n1).
  n0 = numel(theta0);
  n1 = numel(theta);

  % Set a initialization of the variational parameters for each
  % combination of the hyperparameters based on null analysis.
  alpha = repmat(alpha0, [1 1 n1]);
  mu    = repmat(mu0, [1 1 n1]);
 
  % First get the best initialization for the variational parameters.
  fprintf('Finding best initialization for %d combinations ',n0*n1);
  fprintf('of hyperparameters.\n');
  [logw1,alpha,mu,~] = gsea_outerloop(method,betahat,se,SiRiS,snps,sigb,theta0,theta,alpha,mu);

  % Choose an initialization common to all the runs of the coordinate
  % ascent algorithm. This is chosen from the hyperparameters with the
  % highest variational estimate of the importance weight.
  [~,i] = max(logw1(:));
  alpha = repmat(alpha(:,i),[1 n0 n1]);
  mu    = repmat(mu(:,i),[1 n0 n1]);

  % Compute the unnormalized log-importance weights.
  fprintf('Computing importance weights for %d combinations ',n0*n1);
  fprintf('of hyperparameters.\n');
  [logw1,alpha,mu,s] = gsea_outerloop(method,betahat,se,SiRiS,snps,sigb,theta0,theta,alpha,mu);

  % Compute the marginal log-likelihood under the null hypothesis using
  % importance sampling.
  c     = max(logw0(:));
  logZ0 = c + log(mean(exp(logw0(:) - c)));

  % Compute the marginal log-likelihood under the enrichment hypothesis using
  % importance sampling.
  c     = max(logw1(:));
  logZ1 = c + log(mean(exp(logw1(:) - c)));

  % Get the numerical estimate of the Bayes factor.
  log10bf = (logZ1 - logZ0) / log(10);

  % print the BF value.
  fprintf('Log 10 BF = %0.2e\n',log10bf);

end

function [logw1,alpha,mu,s] = gsea_outerloop(method,betahat,se,SiRiS,snps,sigb,theta0,theta,alpha,mu)
% USAGE: the outer loop of the variational inference for the RSS-BVSR model under enrichment
% INPUT:
%       method: the implementation of rss-varbvsr, character 
%       betahat: the effect size estimates under single-SNP model, p by 1 array or C by 1 cell array
%       se: standard errors of betahat, p by 1 array or C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), double precision sparse matrix (ccs format), p by p array or C by 1 cell array
%       snps: the genomic indices of SNPs that are assigned to the pathway
%       sigb: the prior SD of the regression coefficients (if included), scalar
%       theta0: the grid of the genome-wide log-odds (base 10), n0 by 1
%       theta: the grid of the enrichment (base 10), n1 by 1
%	alpha: initial values of the posterior inclusion probabilities, p by n0 by n1
%	mu: initial values of the expected additive effects (given snp included), p by n0 by n1
% OUTPUT:
%       logw1: the unnormalized log-importance weight for each combination of hyper-parameters, n0 by n1
%       alpha: the variational posterior inclusion probabilities, p by n0 by n1
%       mu: the variational posterior means of the additive effects (given snp included), p by n0 by n1
%       s: the variational posterior variances of the additive effects (given snp included), p by n0 by n1

  % Get the number of SNPs analyzed (p),
  % the number of settings of the genome-wide log-odds (n0),
  % and the number of settings of the enrichment parameter (n1).
  if iscell(betahat)
    p = length(cell2mat(betahat));
  else
    p = length(betahat);
  end
  n0 = numel(theta0);
  n1 = numel(theta);

  % Initialize storage for log-importance weights under the alternative
  % hypothesis, and for the variances of the additive effects.
  logw1 = zeros(n0,n1);
  s     = zeros(p,n0,n1);

  % Repeat for each setting of the genome-wide log-odds (theta0).
  iter = 0;
  for i = 1:n0
    
    logodds = (log(10)*theta0(i)) * ones(p,1);

    % Repeat for each setting of the enrichment parameter (theta).
    for j = 1:n1

      iter = iter + 1;
      fprintf('(%03d) theta0 = %+0.2f, theta = %0.2f\n',iter,theta0(i),theta(j));

      options = struct('alpha',alpha(:,i,j),'mu',mu(:,i,j));
      logodds(snps) = log(10) * (theta0(i) + theta(j));
      [logw1(i,j),alpha(:,i,j),mu(:,i,j),s(:,i,j)] = rss_varbvsr_wrapper(method,betahat,se,SiRiS,sigb,logodds,options);
      fprintf('\n');

    end

  end

end

