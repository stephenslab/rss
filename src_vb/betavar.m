% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/betavar.m
% BETAVAR(P,MU,S) returns the variance of X, in which X is drawn from the
% normal distribution with probability P, and X is zero with probability
% 1-P. Inputs MU and S specify the mean and variance of the normal
% density. Inputs P, MU and S must be scalars, or arrays of the same
% dimension. This function is useful for calculating the variance of the
% coefficients under the fully-factorized variational approximation.
function v = betavar (p, mu, s)

  % Note that this is the same as 
  % 
  %    V = P*(S + MU^2) - (P*MU)^2.
  %
  v = p.*(s + (1 - p).*mu.^2);
