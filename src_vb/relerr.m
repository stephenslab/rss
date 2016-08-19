% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/relerr.m
% RELERR(X,Y) returns the absolute relative error between X and Y. Note that
% RELERR(X,Y) = RELERR(Y,X), and we define RELERR([],[]) = 0.
function y = relerr (x1, x2)
  if isempty(x1) | isempty(x2)
    y = 0;
  else
    y = abs(x1 - x2)./(abs(x1) + abs(x2) + eps);
  end
  
