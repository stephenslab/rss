function [lnZ, alpha, mu, s, info] = rss_varbvsr_bigmem_wrapper(method, file, sigb, logodds, options)
% USAGE: a wrapper for the mean-field variational approximation of the RSS-BVSR model given the hyperparameters
%	 this function includes all the implementations that deal with the datasets requiring large memory
% INPUT:
%       method: the implementation of rss-varbvsr, character
%       file: the path of mat file that contains cell arrays of betahat, se and SiRiS, string
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
% 'original': use the original parallel implementation of rss-varbvsr
% 'squarem': add SQAURE to the original bigmem implementation of rss-varbvsr 

  % decide the implementation type of rss-varbvsr
  switch method

    case 'original'
      fprintf('The parallel bigmem implementation of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_bigmem(file,sigb,logodds,options);

    case 'squarem'
      fprintf('The parallel bigmem implementation with SQUAREM add-on of rss-varbvsr is used.\n');
      [lnZ,alpha,mu,s,info] = rss_varbvsr_bigmem_squarem(file,sigb,logodds,options);

    otherwise
      error('Unexpected implementation method.');

  end

end
