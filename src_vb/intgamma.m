% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/intgamma.m
% INTGAMMA(LOGODDS,ALPHA) computes an integral that appears in the
% variational lower bound of the marginal log-likelihood. This integral is
% the expectation on the prior inclusion probabilities taken with respect to
% the variational approximation. Input LOGODDS may either be a scalar, or a
% vector, specifying the prior log-odds. Input ALPHA specifies the mixture
% weights for the variational approximation. See function VARBVS for details
% on these inputs.
function I = intgamma (logodds, alpha)

  % This is the same as 
  %
  %    SUM(ALPHA.*LOG(Q) + (1-ALPHA).*LOG(1-Q)).
  %  
  I = sum((alpha-1) .* logodds + logsigmoid(logodds));
