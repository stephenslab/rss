% parent function: 	rss_ash
% main function: 	update_bz
% local functions: 	bz_gibbs, logexpsum

function [beta_new, bs_new, zlabel_new, labelcounts_new, loglik_new, logpos_new] = update_bz(zs, se, R, beta, bs, zlabel, labelcounts, loglik, logpos, omega, sig_b, matrix_type, p_gamma)
% USAGE: random-scan Gibbs sampler to update (beta, zlabel) under rss-ash model
% INPUT: 
%	zs: 1 by p, betahat ./ se
%	se: 1 by p
%	R: p by p
%	beta: 1 by p, current value of beta
%	bs: 1 by p, current value of beta ./ se
%	zlable: 1 by p, current value of labels
%	labelcounts: 1 by K, label counts of current zlabel
%	loglik: scalar, current value of log p(betahat|beta)
%	logpos: scalar, current value of log p(beta, zlabel|betahat, theta)
%	omega: 1 by K, current value of allocation probabilities
%	sig_b: 1 by K, ordered and non-negative
%	matrix_type: 0 if R is identity; 1 otherwise
%	p_gamma: the re-ordered pmf of 0.3*Unif + 0.7*Geometric based on single-SNP p-values, 1 by p;
% OUTPUT:
%	beta_new: 1 by p, updated beta
%	bs_new: 1 by p, updated beta ./ se
%	zlabel_new: 1 by p, updated labels
%	labelcounts_new: 1 by K, updated label counts
%	loglik_new: scalar, updated value of log p(betahat|beta)
%	logpos_new: scalar, updated value of log p(beta, zlabel|betahat, theta)

	% updated beta is different with current beta in ONLY one entry
	% updated zlabel is different with current zlabel in AT MOST one entry
	% updated labelcounts has either zero or two elements different with current one
	% such LOCAL properties faciliate fast update of loglik and logpos 
	
	% initialize beta, zlabel and labelcounts
	beta_new 	= beta;
	bs_new 		= bs;
	zlabel_new 	= zlabel;
	labelcounts_new = labelcounts;

	% randomly select ONE coordinate to update based on the marginal association
	j = sample(p_gamma);
	
	% gibbs update of the selected coordinate
	switch matrix_type
		case 0 
			linear_tmp = zs(j);
		case 1 
			linear_tmp = zs(j) - ( bs * R(:, j) - bs(j) );
	end

	linear_part = linear_tmp / se(j);
       	[beta_tmp, zlabel_tmp] = bz_gibbs(linear_part, se(j), omega, sig_b);

	% update beta and zlabel
        beta_new(j) 	= beta_tmp;
	bs_new(j) 	= beta_new(j) / se(j);
        zlabel_new(j) 	= zlabel_tmp;

	% update the log likelihood: log p(betahat|beta)
	% p(betahat|beta) = Normal(SRS^{-1}beta, SRS)
	bs_diff 	= bs_new(j) - bs(j);
        bs_square_diff 	= (bs_new(j) + bs(j)) * bs_diff;
	loglik_new 	= loglik + linear_tmp * bs_diff - 0.5 * bs_square_diff;

	% update the log posterior: log p(beta|z, theta)
	% p(beta, z|betahat, theta) is proportional to
        % p(betahat|beta) * p(beta|z, theta) * p(z|theta)
	logprior 	= logpos - loglik;
	logprior_new 	= logprior - 0.5*(beta_new(j)/sig_b(zlabel_new(j)))^2; 
	logprior_new 	= logprior_new + 0.5*(beta(j)/sig_b(zlabel(j)))^2;

	% the proposed label is different with the current one
        if zlabel_new(j) ~= zlabel(j)
                label_rmv = zlabel(j);
                label_add = zlabel_new(j);
		% update label counts
		labelcounts_new(label_rmv) = labelcounts(label_rmv) - 1;
                labelcounts_new(label_add) = labelcounts(label_add) + 1;
		% update log prior
		logprior_new = logprior_new - (log(omega(label_rmv)) - log(sig_b(label_rmv))); 
		logprior_new = logprior_new + (log(omega(label_add)) - log(sig_b(label_add)));  
        end
	logpos_new = loglik_new + logprior_new;
end

function [beta_j, z_j] = bz_gibbs(linear_j, se_j, omega, sig_b)
% USAGE: sample (beta_j, z_j) as a block from its full conditional 
% INPUT:
%       linear_j: scalar, related to the mean term in the full conditional
%       se_j: scalar, the SE of univariate effect size estimate for the jth SNP
%       omega: 1 by K, positive and sum to 1
%       sig_b: 1 by K, ordered and non-negative
% OUTPUT:
%       beta_j: scalar, the gibbs updated value of beta_j
%       z_j: scalar, the gibbs updated value of z_j 

        % STEP 1: sample z_j from categorical distribution
        varsuminv = 1 ./ (se_j^2 + sig_b.^2);
        % use log transform to avoid underflow
        logpmf 	= log(omega) + 0.5*log(varsuminv) - 0.5*(se_j^4)*(linear_j^2) .* varsuminv;
        pmf 	= exp(logpmf - logexpsum(logpmf));
        z_j 	= sample(pmf, 1);

        % STEP 2: given z_j, sample beta_j from normal distribution
        sig = sig_b(z_j);
        if sig == 0
                beta_j = 0;
        else
                sd_beta_j = (sig * se_j) / sqrt(sig^2 + se_j^2);
                mu_beta_j = sd_beta_j^2 * linear_j;
                beta_j 	  = mu_beta_j + sd_beta_j * randn;
        end
end

function V = logexpsum(a)
% USAGE: find the log of sum of exp(a(1)), exp(a(2)),... when a(i) are very LARGE
% INPUT:
%       a: vector, usually very large
% OUTPUT:
%       V: log(sum(exp(a)))

        K = max(a);
        V = log(sum(exp(a-K))) + K;
end
