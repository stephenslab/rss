function [bsam, zsam, wsam, lsam, Naccept] = rss_ash(betahat, se, R, Nsnp, sigma_beta, Ndraw, Nburn, Nthin)
% USAGE: simulate posterior draws under the RSS model with ASH prior
% INPUT:
%	betahat: 1 by p, effect size estimates under single-SNP model
%	se: 1 by p, standard errors of betahat
%	R: p by p, population LD, symmetric, positive definite matrix
%	Nsnp: 1 by p, the number of individuals for each SNP
%	sigma_beta: 1 by K, ordered and non-negative, prior SD on beta
%	Ndraw: scalar, the total number of MCMC draws
%	Nburn: scalar, the length of burn-in
%	Nthin: scalar, the length of thinning interval
% OUTPUT:
%	bsam: Nsam by p, the output chain of beta
%	zsam: Nsam by p, the output chain of zlabel
%	wsam: Nsam by K, the output chain of omega
%	lsam: Nsam by 1, the output chain of lambda
% 	Naccept: scalar, the total number of acceptance in MH step

	% make sure the input summary-level data are row vectors
	if (rows(betahat) ~= 1)
		betahat = betahat';
	end
	if (rows(se) ~= 1)
		se = se';
	end
	if (rows(Nsnp) ~= 1)
		Nsnp = Nsnp';
	end

	p = length(betahat);
	fprintf('total number of variants analyzed: %d \n', p);
	if length(se) ~= p
		error('length of betahat and se must be the same!');
	end

	K = length(sigma_beta);
	fprintf('number of mixture components in ASH prior: %d \n', K);

        Nstart 	  = Nburn + Nthin;
        Nsam 	  = length(Nstart:Nthin:Ndraw);
	tsam 	  = 0;
	Naccept   = 0;

	% preallocate memory for mcmc output
	zsam = uint8(zeros(Nsam, p));
	bsam = zeros(Nsam, p);
	wsam = zeros(Nsam, K);
	lsam = zeros(Nsam, 1);
	
	% pre-compute some quantities used in mcmc iterations
	q  = betahat ./ (se.^2);
	zs = betahat ./ se;

	% lower/upper bound on lambda and step size of proposal
	lambda_LB    = 1e-12;
	lambda_UB    = 1e1;
	sigma_lambda = 0.01;

	% inititate latent label based on marginal z-scores
	abz 			= abs(zs);
	[gamma_start, snp_rank] = initiate_model(p, abz);
	z 			= gamma_start + 1;
	% set starting labels for the large effect
	z(z==2) = round(2*K/3) + 1;
	% initiate the label counts given z
        zc 		   = zeros(1, K);
        zc(1) 		   = sum(z==1);
        zc(round(2*K/3)+1) = p - zc(1);
	% initiate beta given the initial z
        b  = sigma_beta(z) .* randn(1, p);
        bs = b ./ se;
        % initiate lambda by 1 
        l  = 1;
        lt = log((l-lambda_LB)/(lambda_UB-l));
        % initiate omega by Dir(1,...,1)
        w = dirichlet_sample(ones(1, K));
	fprintf('model initiation done \n')

	% rank-based proposal distribution
        p_gamma = pmf_ugmix(p, 2e3);
        p_gamma = p_gamma(snp_rank);

	% if R is an identity matrix, then we have a faster implementation
	matrix_type = 1;
	if R == eye(p)
		matrix_type = 0;
		fprintf('faster computation avaialbe if R is identity matrix \n')
	end

	% viz the mcmc progress
	progress_bar = progress('init','start MCMC iteration');
	
	% the FIRST mcmc step is run outside FOR loop

	% update (beta, zlabel) given data and (omega, lambda)
	[b_loglik, b_logpos] = compute_likpos(q, se, R, sigma_beta, b, z, zc, w, matrix_type);
	[b, bs, z, zc, b_loglik, b_logpos] = update_bz(zs, se, R, b, bs, z, zc, b_loglik, b_logpos, w, sigma_beta, matrix_type, p_gamma);

	% update omega given data and (beta, zlabel, lambda)
        w = update_omega(zc, l, K);

	% update lambda given data and (beta, zlabel, omega)
	[lt_loglik, l] = loglik_theta(lt, w, K, lambda_LB, lambda_UB);
	[l, lt, lt_loglik, Naccept] = update_lambda(w, l, lt, lt_loglik, sigma_lambda, K, lambda_LB, lambda_UB, Naccept);

        if Nburn == 0 && Nthin == 1
		tsam = tsam + 1;

                bsam(tsam, :) = b;
                zsam(tsam, :) = uint8(z);
                wsam(tsam, :) = w;
                lsam(tsam)    = l;
        end
	progress_bar = progress(progress_bar, 1/Ndraw);
	
	%% mcmc rountine to sample (beta, z, omega, lambda)
	for t=2:Ndraw
		
		% update (beta, zlabel) given data and (omega, lambda)
		[b, bs, z, zc, b_loglik, b_logpos] = update_bz(zs, se, R, b, bs, z, zc, b_loglik, b_logpos, w, sigma_beta, matrix_type, p_gamma);
	
		% update omega given data and (beta, zlabel, lambda)
		w = update_omega(zc, l, K);

		% update lambda given data and (beta, zlabel, omega)
		[l, lt, lt_loglik, Naccept] = update_lambda(w, l, lt, lt_loglik, sigma_lambda, K, lambda_LB, lambda_UB, Naccept);

		% record the MCMC simulated draws of parameters
		if t > Nburn && mod(t-Nburn, Nthin) == 0
			tsam = tsam + 1;

			bsam(tsam, :) = b;
			zsam(tsam, :) = uint8(z);
			wsam(tsam, :) = w;
			lsam(tsam)    = l;
		end
	
		progress_bar = progress(progress_bar, t/Ndraw);
	
	end
	
end

function [gamma_start, snp_rank] = initiate_model(p, abz)
% USAGE: assign initial value to gamma for MCMC run by marginal association p-value/z-score
% INPUT:
%       p: the number of snps analyzed
%	abz: the absolute z-score obtained from single-marker model
% OUTPUT:
%       gamma_start: the initial state of the inclusion vector, 0/1, 1 by p
%	snp_rank: the rank of each snp based on its marginal p-value, integer, 1 by p

        gamma_start 	= zeros(1, p);
        q_genome 	= invnormcdf( 1 - (0.025/p) ); % or use matlab stat tool box: norminv 
        in_loci 	= abz > q_genome;
	ngamma_start 	= sum(in_loci);
	[~, snp_rank] 	= sort(abz, 'descend');

	if ngamma_start < 10
		gamma_start(snp_rank(1:10)) = 1;
	else
		gamma_start(in_loci) = 1;
	end
end

function p_gamma = pmf_ugmix(p, geomean)
% USAGE: find pmf of mixture of uniform (0.3) and truncated geometric (0.7) rvs
% INPUT:
% 	p: scalar, specifying the support for the mixture rv on 1, 2, ... p
% 	geomean: the mean parameter of the truncated geometrical rv, scalar
% OUTPUT:
% 	p_gamma: pmf vector, 1 by p

	gp 	  = 1 - exp(-1/geomean); % the same as piMASS
	unif_term = (0.3/p) * ones(1, p);
	geom_seq  = (1-gp) .^ (0:p-1);
	geom_term = 0.7*gp/(1-(1-gp)^p) * geom_seq;
	p_gamma   = unif_term + geom_term;
	p_gamma   = p_gamma / sum(p_gamma);
end

function [loglik, logpos] = compute_likpos(q, se, R, sigmabeta, beta, zlabel, labelcounts, omega, matrix_type)
% USAGE: compute log p(betahat|beta) and log p(beta, zlabel|betahat, theta)
% INPUT:
%       q: 1 by p, betahat ./ (se.^2)
%       se: 1 by p
%       R: p by p correlation matrix
%       sigmabeta: 1 by K, prior SD on beta
%       beta: 1 by p
%       zlabel: 1 by p
%       labelcounts: 1 by K
%       omega: 1 by K, postive and sum to 1
%       matrix_type: 0 if R is identity; 1 otherwise
% OUTPUT:
%       loglik: log p(betahat|beta) up to some constant
%       logpos: log p(beta, z|betahat, theta) up to some constant

        % p(beta, z|betahat, theta) is proportional to
        % p(betahat|beta) * p(beta|z, theta) * p(z|theta)
        logpos = 0;

        % log p(z|theta)
        logpos = logpos + dot(labelcounts, log(omega));

        % log p(beta|z, theta)
        logpos  = logpos - dot(labelcounts, log(sigmabeta));
        betastd = sigmabeta(zlabel);
        logpos  = logpos - 0.5* sum( (beta ./ betastd).^2 );

        % log p(betahat|beta)
        sb = beta ./ se;
        switch matrix_type
                case 0
                        loglik = dot(q, beta) - 0.5*dot(sb, sb);
                case 1
                        loglik = dot(q, beta) - 0.5*(sb * R * sb');
        end
        logpos = logpos + loglik;
end

function omega = update_omega(labelcounts, lambda, K)
% USAGE: Gibbs update of omega
% INPUT:
%       labelcounts: 1 by K
%       lambda: current value of lambda
%       K: the number of label types
% OUTPUT:
%       omega: updated value of omega, 1 by K, positive and sum to 1

        % update the parameter of dirichlet distribution
        para = labelcounts + lambda * ones(1, K);

        % sample omega from the new dirichlet
        omega = dirichlet_sample(para);
        omega = omega ./ sum(omega);
end

function r = dirichlet_sample(a,n)
% Author: Thomas Minka
% Source: http://research.microsoft.com/en-us/um/people/minka/software/fastfit/
% DIRICHLET_SAMPLE   Sample from Dirichlet distribution.
%
% DIRICHLET_SAMPLE(a) returns a probability vector sampled from a 
% Dirichlet distribution with parameter vector A.
% DIRICHLET_SAMPLE(a,n) returns N samples, collected into a matrix, each 
% vector having the same orientation as A.
%
%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986

% This is essentially a generalization of the method for Beta rv's.
% Theorem 4.1, p.594

	if nargin < 2
  		n = 1;
	end

	row = (rows(a) == 1);

	a = a(:);
	% y = gamrnd(repmat(a, 1, n),1);
	% randgamma is faster
	y = randgamma(repmat(a, 1, n));
	r = col_sum(y);
    	% r(find(r == 0)) = 1; % Minka's version
	r(r == 0) = 1;
	r = y./repmat(r, rows(y), 1);
	if row
  		r = r';
	end
end

function [loglik, lambda] = loglik_theta(theta, omega, K, a, b)
% USAGE: transform theta to lambda; compute log p(theta|omega)
% INPUT:
%       theta: scalar, logit transformation of lambda
%       omega: 1 by K, all positive and sum 1
%       K: scalar
%       a & b: lower and upper bound of lambda
% OUTPUT:
%       loglik: scalar, log p(theta|omega)
%       lambda: scalar

        lambda = (a + b*exp(theta)) / (1 + exp(theta));

        loglik = gammaln(K*lambda) - K * gammaln(lambda);
        loglik = loglik + lambda * sum(log(omega));
        loglik = loglik + theta - 2 * log(1 + exp(theta));
end

function [lambda, theta, theta_loglik, Naccept] = update_lambda(omega, lambda_old, theta_old, loglik_old, sigma_theta, K, a, b, Naccept)
% USAGE: update theta:=logit(lambda); use Metropolis-Hastings (MH) normal random walk
% INPUT:
%       omega: 1 by K, positive and sum to 1
%       lambda_old: scalar, the current value of lambda
%       theta_old: scalar, logit transform of lambda_old
%       loglik_old: scalar, the log likelihood at theta_old
%       sigma_theta: scalar, the SD of the normal random walk
%       K: number of mixture components
%       a & b: the lower and upper bounds of lambda
%       Naccept: the total number of acceptance up til now
% OUTPUT:
%       lambda: the inverse logit transform of theta
%       theta: the MH updated value of theta
%       theta_loglik: the log likelihood at updated value of theta
%       Naccept: the total number of acceptance after this MH update 

        if rand < 0.33
                extrastep = randperm(19);
                repeat = 1 + extrastep(1);
        else
                repeat = 1;
        end

        theta_new = propose_theta(theta_old, sigma_theta, repeat);
        [loglik_new, lambda_new] = loglik_theta(theta_new, omega, K, a, b);
        theta_ratio = loglik_new - loglik_old;

        if theta_ratio > 0 || log(rand) < theta_ratio
                Naccept = Naccept + 1;
                lambda = lambda_new;
                theta = theta_new;
                theta_loglik = loglik_new;
        else
                lambda = lambda_old;
                theta = theta_old;
                theta_loglik = loglik_old;
        end
end

function theta_new = propose_theta(theta_old, sigma_theta, repeat)
% USAGE: propose the gaussian random walk update for theta in the MH step
% INPUT:
%       theta_old: the current value of theta, scalar
%       sigma_theta: the sd of the normal move added to theta_old, scalar, positive
%       repeat: the number of local proposals to compound, scalar
% OUTPUT:
%       theta_new: the newly proposed value of theta, scalar
% NB: this is a symmetric jump, so the contribution to log MH ratio is ZERO

        theta = theta_old;
        for i=1:repeat
                theta = theta + randn * sigma_theta; % add a normal error and reflect about the boundary if outside
        end

        theta_new = theta;
end


