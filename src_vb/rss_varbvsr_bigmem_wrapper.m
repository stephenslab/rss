function [lnZ, alpha, mu, s, info] = rss_varbvsr_bigmem_wrapper(method, file, sigb, logodds, options)
% USAGE: mean-field variational inference of RSS-BVSR model for a given set of hyperparameters
%        this function includes all implementations that deal with datasets requiring large memory
% INPUT:
%       method: the implementation of rss-varbvsr, string
%       file: the mat file that contains cell arrays of betahat, se and SiRiS, string
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
% 'original': use the original parallel implementation of rss-varbvsr
% 'squarem': add SQAURE to the original bigmem implementation of rss-varbvsr 

  % decide the implementation type of rss-varbvsr
  switch method

    case 'original'
      fprintf('The parallel bigmem implementation of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_bigmem(file,sigb,logodds,options);

    case 'squarem'
      fprintf('The parallel bigmem SQUAREM implementation of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_bigmem_squarem(file,sigb,logodds,options);

    otherwise
      error('Unexpected implementation method.');

  end

end
