function [loo_mean, loo_var] = mvnloo(x_vec, m_vec, v_mat)
% USAGE: compute the univariate complete conditional of MVN(m_vec, v_mat)
% INPUT:
%	x_vec: p by 1 vector, the observation from the multivariate normal
%	m_vec: p by 1 vector, the mean of the multivariate normal
%	v_mat: p by p matrix, the covariance matrix 
% OUTPUT:
%	loo_mean: p by 1 vector, the complete conditional mean
%	loo_var: p by 1 vector, the complete conditional variance
% SOURCE:
% https://www2.stat.duke.edu/~km68/materials/214.3%20(MVN%20Theory).pdf

  % compute the precision matrix
  p     = length(m_vec);
  k_mat = v_mat \ speye(p);
  disp('precision matrix is ready ...');

  % get the complete conditional variance
  loo_var = 1 ./ diag(k_mat);

  % compute the complete conditional mean
  loo_mean = x_vec - (k_mat * (x_vec-m_vec)) .* loo_var;  
 
end
