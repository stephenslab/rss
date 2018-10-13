function I = intklbeta_rssbvsr(alpha, mu, s, sigb_square)
% USAGE: compute an integral in the variational lower bound of the marginal log-likelihood
% INPUT:
%	alpha: p by 1, variational estimates of the posterior inclusion probabilities
%	mu: p by 1, posterior means of the additive effects (given snp included)
%	s: p by 1, posterior variances of the additive effects (given snp included)
%	sigb_square: p by 1 or scalar, prior variance of the regression coefficients (if included)
% OUTPUT:
%	I: negative KL divergence between the approximation and the prior of regression coefficients
 
  I = sum(alpha) + alpha'*log(s./sigb_square);
  I = (I - alpha'*((s+mu.^2)./sigb_square)) * 0.5;
  I = I - alpha'*log(alpha+eps) - (1-alpha)'*log(1-alpha+eps);

end  
