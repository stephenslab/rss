% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/logpexp.m
% LOGPEXP(X) returns LOG(1 + EXP(X)). For large X, LOGPEXP should be
% approximately X. The computation is performed in a numerically stable
% manner.
function y = logpexp (x)

  % For large entries of X, LOG(1 + EXP(X)) is effectively the same as X.
  y = x;

  % Find entries of X that are not large. For these entries, compute
  % LOG(1 + EXP(X)).
  I    = find(x < 16);
  y(I) = log(1 + exp(x(I)));
