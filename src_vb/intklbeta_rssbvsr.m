function I = intklbeta_rssbvsr(alpha, mu, sigma_square, sigma_beta_square)
% USAGE: compute an integral in the variational lower bound of the marginal log-likelihood.
% INPUT: see rss_varbvsr for details on the inputs to this function.
% OUTPUT: the negative Kullback-Leibler divergence between the approximation and the prior of the coefficients.
  
  I = (sum(alpha) + alpha'*log(sigma_square/sigma_beta_square) ... 
      - alpha'*(sigma_square + mu.^2)/sigma_beta_square) * 0.5 ...
      - alpha'*log(alpha + eps) - (1 - alpha)'*log(1 - alpha + eps);

end  
