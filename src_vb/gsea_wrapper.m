function [log10bf,logw1,alpha,mu,s] = gsea_wrapper(method,betahat,se,SiRiS,se_all,Nsnp_all,snps,h,theta0,theta,logw0,alpha0,mu0,alpha,mu)
% USAGE: perform the full variational inference for the RSS-BVSR model under alternative hypothesis
%        alternative hypothesis: the gene set is enriched for genotype-phenotype association
% INPUT:
%       method: the implementation of rss-varbvsr, character 
%       betahat: the effect size estimates under single-SNP model, p by 1 array or C by 1 cell array
%       se: standard errors of betahat, p by 1 array or C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), double precision sparse matrix (ccs format), p by p array or C by 1 cell array
%	se_all: standard errors for all SNPs in the genome, numsnp by 1
%	Nsnp_all: sample sizes for all SNPs in the genome, numsnp by 1
%	snps: the genomic indices of SNPs that are assigned to the pathway, p by 1 
%       h: the fixed proportion of phenotypic variance explained by available genotypes, nh by 1
%       theta0: the grid of the genome-wide log-odds (base 10), n0 by 1
%       theta: the grid of the enrichment (base 10), n1 by 1
%       logw0: the log-importance weights for the grid of theta0 under null, n0 by nh
%       alpha0: the variational posterior inclusion probabilities under null, p by n0 by nh
%	mu0: the variational posterior means of the additive effects (given snp included), p by n0 by nh
%       alpha: the initial values of the posterior inclusion probabilities under erichment, p by n0 by n1 by nh 
%       mu: the initial values of the expected additive effects (given snp included) under enrichment, p by n0 by n1 by nh
% OUTPUT:
%	log10bf: log 10 of gene set enrichment BF based on summary-level data, scalar
%       logw1: the log-importance weights for the grid of (theta0,theta,h) under enrichment, n0 by n1 by nh
%       alpha: the variational posterior inclusion probabilities under enrichment, p by n0 by n1 by nh
%       mu: the variational posterior means of the additive effects under enrichment, p by n0 by n1 by nh
%       s: the variational posterior variances of the additive effects under enrichment, p by n0 by n1 by nh

% NOTE: 
% To calculate the BF, we integrate over the hyperparameters using importance sampling under the
% assumption that SNPs outside the enriched pathways are unaffected by pathways (a posteriori).
% See pp 5-6 of Supplementary text of Carbonetto and Stephens (PLoS Genetics, 2013).

  % get the number of settings of the genome-wide log-odds (n0)
  n0 = numel(theta0);
  % get the number of settings of the enrichment parameter (n1)
  n1 = numel(theta);
  % get the number of settings of the heritability parameter (nh)
  nh = numel(h);

  % sanity check on subsetting the se data
  if iscell(se)
    se_vec = cell2mat(se);
  else
    se_vec = se;
  end
  if all(se_vec==se_all(snps)) ~= 1
    error('Inconsistent data input ...');
  end
  clear se_vec;

  % first get the best initialization for the variational parameters
  fprintf('Finding best initialization for %d combinations \n',n0*n1*nh);
  fprintf('of hyperparameters.\n');
  [logw1,alpha,mu] = gsea_outerloop(method,betahat,se,SiRiS,se_all,Nsnp_all,snps,h,theta0,theta,alpha,mu,logw0,alpha0,mu0);

  % next choose an initialization common to all the runs of the coordinate ascent algorithm
  % this is chosen from the hyperparameters with the highest variational estimate of the log-importance weight
  [~, i] = max(logw1(:));
  alpha  = repmat(alpha(:,i),[1 n0 n1 nh]);
  mu     = repmat(mu(:,i),[1 n0 n1 nh]);

  % compute the unnormalized log-importance weights
  fprintf('Computing importance weights for %d combinations \n',n0*n1*nh);
  fprintf('of hyperparameters.\n');
  [logw1,alpha,mu,s] = gsea_outerloop(method,betahat,se,SiRiS,se_all,Nsnp_all,snps,h,theta0,theta,alpha,mu,logw0,alpha0,mu0);

  % compute the marginal log-likelihood under null using importance sampling
  c     = max(logw0(:));
  logZ0 = c + log(mean(exp(logw0(:) - c)));

  % compute the marginal log-likelihood under enrichment using importance sampling
  c     = max(logw1(:));
  logZ1 = c + log(mean(exp(logw1(:) - c)));

  % get the numerical estimate of the Bayes Factor (BF)
  log10bf = (logZ1 - logZ0) / log(10);

  % print the BF value.
  fprintf('Log 10 BF = %0.2e\n',log10bf);

end

function [logw1,alpha,mu,s] = gsea_outerloop(method,betahat,se,SiRiS,se_all,Nsnp_all,snps,h,theta0,theta,alpha,mu,logw0,alpha0,mu0)
% USAGE: the outer loop of the variational inference for the RSS-BVSR model under enrichment
% INPUT:
%       method: the implementation of rss-varbvsr, character 
%       betahat: the effect size estimates under single-SNP model, p by 1 array or C by 1 cell array
%       se: standard errors of betahat, p by 1 array or C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), double precision sparse matrix (ccs format), p by p array or C by 1 cell array
%       se_all: standard errors for all SNPs in the genome, numsnp by 1
%       Nsnp_all: sample sizes for all SNPs in the genome, numsnp by 1
%       snps: the genomic indices of SNPs that are assigned to the pathway, p by 1 
%       h: the fixed proportion of phenotypic variance explained by available genotypes, nh by 1
%       theta0: the grid of the genome-wide log-odds, n0 by 1
%       theta: the grid of the enrichment, n1 by 1
%       logw0: the log-importance weights for the grid of theta0 under null, n0 by nh
%       alpha0: the variational posterior inclusion probabilities under null, p by n0 by nh
%       mu0: the variational posterior means of the additive effects (given snp included), p by n0 by nh
%       alpha: the initial values of the posterior inclusion probabilities under erichment, p by n0 by n1 by nh 
%       mu: the initial values of the expected additive effects (given snp included) under enrichment, p by n0 by n1 by nh
% OUTPUT:
%       logw1: the log-importance weights for the grid of (theta0,theta,h) under enrichment, n0 by n1 by nh
%       alpha: the variational posterior inclusion probabilities under enrichment, p by n0 by n1 by nh
%       mu: the variational posterior means of the additive effects under enrichment, p by n0 by n1 by nh
%       s: the variational posterior variances of the additive effects under enrichment, p by n0 by n1 by nh

  % get the number of SNPs assigned to the enriched pathway (p)
  if iscell(betahat)
    p = length(cell2mat(betahat));
  else
    p = length(betahat);
  end
  % get the number of settings of the genome-wide log-odds (n0)
  n0 = numel(theta0);
  % get the number of settings of the enrichment parameter (n1)
  n1 = numel(theta);
  % get the number of settings of the heritability parameter (nh)
  nh = numel(h);
  % get the number of SNPs in the whole genome (ng)
  ng = length(se_all);

  % initialize storage for log-importance weights under the alternative
  % hypothesis, and for the variances of the additive effects
  logw1 = zeros(n0,n1,nh);
  s     = zeros(p,n0,n1,nh);

  iter = 0;
  % repeat for each setting of the heritability (h)
  for k = 1:nh
    % repeat for each setting of the genome-wide log-odds (theta0)
    for i = 1:n0
      % NOTE:
      % First compute the variational lower bound to the marginal log-likelihood
      % under the null hypothesis when there is no enrichment, only for SNPs in
      % the gene set. Note that rss_varbvsr defines the log-odds ratio using the
      % natural logarithm, so we need to multiply theta0 by log(10).
      log10odds_null = theta0(i) * ones(ng,1); % genome-wide baseline  

      sigb0    = calc_beta_sd(se_all, Nsnp_all, sigmoid10(log10odds_null), h(k));
      logodds0 = log(10) * theta0(i);
      options  = struct('alpha',alpha0(:,i,k),'mu',mu0(:,i,k),'verbose',true);

      % run rss-varbvsr under null
      F0 = rss_varbvsr_wrapper(method,betahat,se,SiRiS,sigb0,logodds0,options);

      % repeat for each setting of the enrichment parameter (theta)
      for j = 1:n1

        % display the current status of the inference procedure
        iter = iter + 1;
        fprintf('(%04d) theta0 = %+0.2f, theta = %0.2f, h = %0.2f \n',iter,theta0(i),theta(j),h(k));

        % NOTE:
        % Next compute the marginal log-likelihood under the alternative
        % hypothesis only for SNPs in the pathway. Note that rss_varbvsr
        % defines the log-odds ratio using the natural logarithm, so we
        % need to multiply (theta0+theta) by log(10).
	log10odds_gsea       = theta0(i) * ones(ng,1);  % genome-wide baseline
      	log10odds_gsea(snps) = theta0(i) + theta(j);    % enriched in the pathway

	sigb1    = calc_beta_sd(se_all, Nsnp_all, sigmoid10(log10odds_gsea), h(k));
        logodds1 = log(10) * (theta0(i) + theta(j));
        options  = struct('alpha',alpha(:,i,j,k),'mu',mu(:,i,j,k),'verbose',true);

        % run rss-varbvsr under enrichment
        [F1,alpha(:,i,j,k),mu(:,i,j,k),s(:,i,j,k)] = rss_varbvsr_wrapper(method,betahat,se,SiRiS,sigb1,logodds1,options);

        % NOTE:
        % Compute the importance weight under the alternative. The prior and
        % proposal for the genome-wide log-odds and enrichment parameters are
        % assumed to be uniform, so they cancel out from the importance
        % weights. For an explanation why we can decompose the variational
        % lower bound in this way, see pp 5-6 of Supplementary text of 
        % Carbonetto and Stephens (PLoS Genetics, 2013).
        logw1(i,j,k) = logw0(i,k) + F1 - F0;
      end
    end
  end
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

