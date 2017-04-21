function [lnZ, alpha, mu, s, info] = rss_varbvsr_wrapper(method, betahat, se, SiRiS, sigb, logodds, options)
% USAGE: a wrapper for the mean-field variational approximation of the RSS-BVSR model given the hyperparameters
% INPUT:
%       method: the implementation of rss-varbvsr, character
%       betahat: the effect size estimates under single-SNP model, p by 1 array or C by 1 cell array
%       se: standard errors of betahat, p by 1 array or C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), double precision sparse matrix (ccs format), p by p array or C by 1 cell array
%       sigb: the prior SD of the regression coefficients (if included), scalar
%       logodds: the prior log-odds (i.e. log(prior PIP/(1-prior PIP))) of inclusion for each SNP, p by 1
%       options: user-specified behaviour of the algorithm, structure
% OUTPUT:
%       lnZ: scalar, the variational lower bound of the marginal log likelihood (up to some constant)
%       alpha: p by 1, variational estimates of the posterior inclusion probabilities 
%       mu: p by 1, posterior means of the additive effects (given snp included)
%       s: p by 1, posterior variances of the additive effects (given snp included)
%       info: structure with following fields 
%               - iter: integer, number of iterations
%               - maxerr: the maximum relative difference between the parameters at the last two iterations
%               - sigb: scalar, the maximum likelihood estimate of sigma_beta
%               - loglik: iter by 1, the variational lower bound at each iteration

% IMPLEMENTATION METHODS
% 'original': the summary data must be p-by-1 vectors and p-by-p matrix
% 'squarem': the summary data must be p-by-1 vectors and p-by-p matrix; use SQUAREM add-on
% 'parallel': the summary data must be C-by-1 cell arrays; use parallel for loop
% 'pasquarem': the summary data must be C-by-1 cell arrays; use parallel for loop and SQUAREM

  % make sure summary-level data are stored as C by 1 cell arrays
  cell_check = prod([iscell(betahat),iscell(se),iscell(SiRiS)]);
  if cell_check == 1
    fprintf('Summary-level data are stored as cell arrays. \n');
  end

  C = length(betahat);
  if prod([length(se)==C,length(SiRiS)==C]) == 0
    error('Summary-level data must have the same dimension.');
  end

  % decide the implementation type of rss-varbvsr
  switch method

    case 'original'
      if C ~= 1 && cell_check == 1
        error('Use parallel implementations.');
      end
      if cell_check == 1
	betahat = cell2mat(betahat);
	se 	= cell2mat(se);
	SiRiS 	= cell2mat(SiRiS);
      end

      fprintf('The original implementation of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr(betahat,se,SiRiS,sigb,logodds,options);

    case 'squarem'
      if C ~= 1 && cell_check == 1
	error('Use parallel implementations.');
      end
      if cell_check == 1
        betahat = cell2mat(betahat);
        se      = cell2mat(se);
        SiRiS   = cell2mat(SiRiS);
      end

      fprintf('The serial implementation with SQUAREM add-on of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_squarem(betahat,se,SiRiS,sigb,logodds,options);

    case 'parallel'
      if cell_check ~= 1
        error('Use serial implementations.');
      end

      fprintf('The parallel implementation of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_parallel(betahat,se,SiRiS,sigb,logodds,options);

    case 'pasquarem'
      if cell_check ~= 1
        error('Use serial implementations.');
      end

      fprintf('The parallel implementation with SQUAREM add-on of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_pasquarem(betahat,se,SiRiS,sigb,logodds,options);

    otherwise
      error('Unexpected implementation method.');

  end

end
