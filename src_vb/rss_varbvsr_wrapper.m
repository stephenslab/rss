function [lnZ,alpha,mu,s,info] = rss_varbvsr_wrapper(method,betahat,se,SiRiS,sigb,logodds,options)
% USAGE: mean-field variational inference of RSS-BVSR model for a given set of hyperparameters
% INPUT:
%       method: the implementation of rss-varbvsr, string
%       betahat: effect size estimates under single-SNP model, p by 1 array or C by 1 cell array
%       se: standard errors of betahat, p by 1 array or C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), sparse matrix (CCS format), p by p array or C by 1 cell array
%       sigb: prior SDs of the regression coefficients (if included), p by 1 or scalar
%       logodds: log(prior PIP/(1-prior PIP)) of inclusion for each SNP, p by 1 or scalar
%       options: user-specified behaviour for each implementation, structure
% OUTPUT:
%       lnZ: scalar, variational lower bound of the marginal log likelihood (up to some constant)
%       alpha: p by 1, variational estimates of the posterior inclusion probabilities 
%       mu: p by 1, posterior means of the additive effects (if the SNP is included)
%       s: p by 1, posterior variances of the additive effects (if the SNP is included)
%       info: structure with following fields 
%               - iter: integer, number of iterations till convergence
%               - maxerr: maximum relative difference between the parameters at the last 2 iterations
%               - loglik: iter by 1, variational lower bound at each iteration

% IMPLEMENTATION METHODS
% 'original': input data must be p-by-1 vectors and p-by-p matrix
% 'squarem': input data must be p-by-1 vectors and p-by-p matrix; only use SQUAREM
% 'parallel': input data must be C-by-1 cell arrays; only use parallel for loop
% 'pasquarem': input data must be C-by-1 cell arrays; use both parallel for loop and SQUAREM

  % check whether input data are stored as C by 1 cell arrays
  cell_check = prod([iscell(betahat),iscell(se),iscell(SiRiS)]);
  if cell_check == 1
    fprintf('Input data are stored as cell arrays. \n');
  end

  C = length(betahat);
  if prod([length(se)==C,length(SiRiS)==C]) == 0
    error('Input data must have consistent dimensions.');
  end

  % decide the implementation type of rss-varbvsr
  switch method

    case 'original'
      if C ~= 1 && cell_check == 1
        error('Use parallel implementations instead.');
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
	error('Use parallel implementations instead.');
      end
      if cell_check == 1
        betahat = cell2mat(betahat);
        se      = cell2mat(se);
        SiRiS   = cell2mat(SiRiS);
      end

      fprintf('The serial SQUAREM implementation of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_squarem(betahat,se,SiRiS,sigb,logodds,options);

    case 'parallel'
      if cell_check ~= 1
        error('Use serial implementations instead.');
      end

      fprintf('The parallel implementation of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_parallel(betahat,se,SiRiS,sigb,logodds,options);

    case 'pasquarem'
      if cell_check ~= 1
        error('Use serial implementations instead.');
      end

      fprintf('The parallel SQUAREM implementation of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_pasquarem(betahat,se,SiRiS,sigb,logodds,options);

    otherwise
      error('Unexpected implementation method.');

  end

end
