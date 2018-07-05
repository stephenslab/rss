function [betasam, gammasam, hsam, logpisam, msam, ngsam, Naccept] = rss_bvsr_mdp(betahat, se, R, Nsnp, Ndraw, Nburn, Nthin)
% USAGE: simulate posterior draws under RSS model with BVSR prior (allowing for multiplicative dispersion)
% INPUT:
%	betahat: the marginal effects estimates, p by 1
%	se: the standard errors of betahat, p by 1
%	q: betahat ./ se.^2, p by 1
%	R: symmetric, positive definite matrix, p by p
%	Nsnp: sample size of SNPs, p by 1
%	Ndraw: the total number of MCMC samples to draw, integer
%	Nburn: the number of initial samples to discard ("burn-in"), integer
%	Nthin: keep every Nthin-th simulation draw from the sequence ("thin"), integer
% OUTPUT:
%	loglik: the log likelihood, Ndraw by 1
%	Naccept: the number of accepted moves in MH step, scalar
%	hsam: the MCMC sample of h, Nsam by 1
%	gammasam: the MCMC sample of z, Nsam by p
%	logpisam: the MCMC sample of log(pi), Nsam by 1
%	betasam: the MCMC sample of beta, Nsam by p
%	msam: the MCMC sample of m, Nsam by 1
%	ngsam: the MCMC sample of |gamma|, Nsam by 1

	% make sure the input summary-level data are column vectors
	betahat = betahat(:);
	se 	= se(:);
	Nsnp 	= Nsnp(:);
	p 	= length(betahat);
	fprintf('total number of variants analyzed: %d \n', p);
	if length(se) ~= p
		error('length of betahat and se must be the same!')
	end

	Naccept = 0;
	Nstart 	= Nburn + Nthin;
        Nsam 	= length(Nstart:Nthin:Ndraw);
	tsam 	= 0;

	% pre-compute some quantities used in mcmc iterations
	xxyy 	= 1 ./ ( Nsnp .* (se.^2) );
        q 	= betahat ./ (se.^2);
        Si 	= 1 ./ se;
	tic;
        SiRiS 	= repmat(Si, 1, p) .* R .* repmat(Si', p, 1);
	clear R;
	L 	= chol(SiRiS, 'lower');
	Liq 	= L \ q;
	clear L;
	bAb 	= dot(Liq, Liq);
	clear Liq;
	mtime = toc;
	fprintf('heavy matrix computation completed after %d hours \n', round(mtime/3600));
 
	% preallocate memory for the output matrices
	gammasam = zeros(Nsam, p);
    	hsam 	 = zeros(Nsam, 1);
    	betasam  = zeros(Nsam, p);
    	logpisam = zeros(Nsam, 1);
    	msam 	 = zeros(Nsam, 1);
	ngsam 	 = zeros(Nsam, 1);

	% initiate model
	abz = abs(betahat ./ se);
	[gamma_start, ngamma_start, snp_rank] = initiate_model(p, abz);
	gamma_old = gamma_start; 					% 0/1, 1 by p 
	rank_old  = find(gamma_old); 					% integer, 1 by length(gamma_old == 1)
	logpi_old = log(ngamma_start / p);
	h_old 	  = rand;
	m_old	  = 1+rand; 						% m is close to 1 for most of autosomes 
	fprintf('model initiation done \n')

	% rank-based proposal distribution
	p_gamma = pmf_ugmix(p, 2e3);
	p_gamma = p_gamma(snp_rank);

	% viz the mcmc progress
	progress_bar = progress('init','start MCMC iteration');

	% calculate the posterior p(h gamma pi|betahat) for the FIRST iteration
	psi_old = calc_beta_variance(xxyy, logpi_old, h_old);
	[betapost_old, ~, logpost_old] = calc_posterior_bvsr_mdp(q, SiRiS, bAb, rank_old, psi_old, logpi_old, m_old);
	beta_old = zeros(1, p); beta_old(rank_old) = betapost_old;
 
for t = 1:Ndraw
	
	% MH sampler for [gamma, h, log(pi)] with small world proposal
	if rand < 0.33
		extrastep = randperm(19);
		repeat 	  = 1 + extrastep(1);
	else
		repeat = 1; 			
	end	
	
	logMHratio = 0;
	
	% symmetric uniform random walk proposal for h
	h_new 	   = propose_h(h_old, repeat);
	logMHratio = logMHratio + 0;

	% symmetric normal random walk proposal for m
	m_new 	   = m_old + 1e-2*randn;
	logMHratio = logMHratio + 0;  

	% 'rank based proposal' for gamma (Guan & Stephens 2011)
	[rank_new, gamma_logratio] = propose_gamma(rank_old, p_gamma, repeat);
	logMHratio 		   = logMHratio + gamma_logratio;
	
	% symmetric uniform random walk proposal for log(pi)
	[logpi_new, pi_logratio] = propose_logpi(logpi_old, repeat, p);
	logMHratio 		 = logMHratio + pi_logratio;

	% calculate the posterior p(h gamma pi|betahat) for the proposed value
	psi_new = calc_beta_variance(xxyy, logpi_new, h_new);
        [betapost_new, ~, logpost_new] = calc_posterior_bvsr_mdp(q, SiRiS, bAb, rank_new, psi_new, logpi_new, m_new);
	beta_new = zeros(1, p); beta_new(rank_new) = betapost_new;

	logMHratio = logMHratio + logpost_new - logpost_old; % from the posterior
	
	% make a choice for MH step
	if logMHratio>0 || log(rand)<logMHratio
		Naccept 	= Naccept + 1;
		logpost_old 	= logpost_new;
		rank_old 	= rank_new;
		beta_old 	= beta_new;
		h_old 		= h_new;
		logpi_old 	= logpi_new;
		m_old 		= m_new;	
	end
	rank_new = []; beta_new = []; betapsm_new = [];  %#ok<NASGU>
	
	gamma_old = zeros(1,p); gamma_old(rank_old) = 1;

	% save the MCMC sample of [beta gamma h log(pi)] 
	if t > Nburn && mod(t-Nburn, Nthin) == 0
		tsam 			= tsam + 1;
		gammasam(tsam, :) 	= gamma_old; 
		hsam(tsam) 		= h_old; 
		betasam(tsam, :) 	= beta_old; 
		logpisam(tsam) 		= logpi_old; 
		ngsam(tsam) 		= length(rank_old); 
		msam(tsam) 		= m_old;
	end
	
	progress_bar = progress(progress_bar, t/Ndraw);

end

end

function [gamma_start, ngamma_start, snp_rank] = initiate_model(p, abz)
% USAGE: initiate latent labels by marginal association absolute z-score
% INPUT:
%	p: the number of snps analyzed 
%	abz: the absolute z-score obtained from single-marker model 
% OUTPUT:
%	gamma_start: 1 by p, the initial latent labels
%	ngamma_start: integer, the initial number of snps included in the model
%	snp_rank: integer, 1 by p, the rank of snps based their marginal p-values

	gamma_start = zeros(1, p);
	q_genome = invnormcdf( 1 - (0.025/p) ); % matlab stat tool box: norminv 
	in_loci = (abz > q_genome);
	ngamma_start = sum(in_loci);
	[~, snp_rank] = sort(abz, 'descend');
	baseline = min(10, p-1); % avoid that initial pi == 1 (-Inf log likelihood)
	if ngamma_start < baseline
		ngamma_start = baseline;
		gamma_start(snp_rank(1:baseline)) = 1;
	else
		gamma_start(in_loci) = 1;
	end
end

function p_gamma = pmf_ugmix(p, geomean)
% USAGE: find pmf of mixture of uniform (0.3) and truncated geometric (0.7) rvs
% INPUT:
%	p: scalar, specifying the support for the mixture rv on 1, 2, ... p
%	geomean: the mean parameter of the truncated geometrical rv, scalar
% OUTPUT:
%	p_gamma: pmf vector, 1 by p

	gp = 1 - exp(-1/geomean); % the same as piMASS 
	unif_term = (0.3/p) * ones(1, p);
	geom_seq = (1-gp) .^ (0:p-1); 
	geom_term = 0.7*gp/(1-(1-gp)^p) * geom_seq;
	p_gamma = unif_term + geom_term;
	p_gamma = p_gamma / sum(p_gamma);
end

function psi = calc_beta_variance(xxyy, logpi, h)
% USAGE: calculate sigma_beta^2 in the same way as BSLMM
% INPUT:
%	xxyy: scalar
%	logpi: scalar
%	h: scalar
% OUTPUT: psi: scalar, variance of beta, sigma_beta^2
	
	xxyysum = sum(xxyy);
	pival = exp(logpi);
	psi = h / (pival * xxyysum);
end

function h_new = propose_h(h_old, repeat)
% USAGE: symmetric uniform random walk proposal for h
%               h_new = h_old + Unif(-0.1, 0.1) in the range [0,1]
% INPUT:
%       h_old: the current value of h, scalar
%       repeat: small-world proposal, scalar
% OUTPUT:
%       h_new: the proposed value of h, scalar

        h = h_old;
        for i=1:repeat
                h = h + (rand-0.5) * 0.2;
                while true
                        if h < 0
                                h = 2*0 - h;
                        end
                        if h > 1
                                h = 2*1 - h;
                        end
                        if (h >= 0 && h <= 1); break; end % until the proposed value is inside [0, 1]
                end
        end
        h_new = h;
% NB: this is a symmetric jump, so its contribution to log MH ratio is ZERO
end

function [logpi_new, pi_logratio] = propose_logpi(logpi_old, repeat, p)
% USAGE: symmetric uniform random walk proposal for log(pi)
%               log(pi_new) = log(pi_old) + Unif(-0.05, 0.05) range log(1/p, 1)
% INPUT:
%       logpi_old: the current value of log(pi), scalar
%       repeat: small-world proposal, integer
%       p: the total number of snps, integer
% OUTPUT:
%       logpi_new: the proposed value of log(pi), scalar
%       pi_logratio: the contribution of proposal for log(pi) to log MH acceptance ratio
% NB: replace 1 with 1-1e-8 to avoid drawing pi == 1, which leads to -Inf log likelihood

        pi_logratio = 0;
        for i=1:repeat
                logpi_new = logpi_old + (rand-0.5)*0.1;
		while true
                	if logpi_new < log(1/p)
                        	logpi_new = 2*log(1/p) - logpi_new;
                	end
                	if logpi_new > log(1-1e-8)
                        	logpi_new = 2*log(1) - logpi_new;
                	end
			if (logpi_new >= log(1/p) && logpi_new < log(1-1e-8)); break; end
		end
                pi_logratio = pi_logratio + (logpi_new-logpi_old);
                logpi_old = logpi_new;
        end
end

function [betapost, loglik, logpost] = calc_posterior_bvsr_mdp(q, SiRiS, bAb, rank, psi, logpi, m)
% USAGE: calculate p(h gamma pi m|betahat) up to some constant 
% INPUT:
%       q: betahat ./ se^2, p by 1
%       SiRiS: S^{-1}RS^{-1}, p by p, sparse matrix
%	bAb: betahat' * (SRS)^{-1} * betahat, scalar
%       rank: the marginal rank of snps included in the current model, 1 by n_gamma
%       logpi: the current value of log(pi), scalar with range log(1/p,1)
%       psi: the current value of prior variance of beta, scalar
%	m: the multiplicative dispersion parameter, scalar
% OUTPUT:
%       logpost: the log posterior p(h gamma pi|betahat), scalar
%       betapost: sample positive components of beta from its conditional posterior, 1 by n_gamma
%       loglik: the log likelihood p(betahat|h gamma), scalar

        logpost = 0;
        n_gamma = length(rank); 		% the number of snps in the current model, integer
        p 	= length(q); 			% the total number of snps, integer
        Psiz 	= psi * ones(n_gamma, 1); 	% the vector of current value of psi(gamma, h), n_gamma by 1

        if n_gamma ~= 0

                RSBSMz 	  = q(rank); 			% n_gamma by 1
                OmegaInvz = SiRiS(rank, rank) / m; 	% n_gamma by n_gamma

		% add a diagonal matrix Psiz to OmegaInvz
                OmegaInvz(1:(n_gamma+1):n_gamma^2) = OmegaInvz(1:(n_gamma+1):n_gamma^2) + 1 ./ Psiz';

                L   = chol(OmegaInvz, 'lower'); 
                U   = L';
                T   = solve_tril(L, RSBSMz) / m;
                muz = solve_triu(U, T);

                % sample positive components of beta from Normal(muz, Omega)
                betapost_err = solve_triu(U, randn(n_gamma,1));
                betapost     = muz' + betapost_err'; % convert the column to row vector

                % calculate log posterior for MH
                % p(h gamma pi m|betahat) = Const. * p(betahat|h gamma) * p(h) * p(gamma|pi) * p(pi) * p(m) 
                
		% first compute likelihood p(betahat| h gamma)
                logdet = - sum(log(diag(L))) - 0.5 * sum(log(Psiz)) - 0.5*p*log(m);
                logquad = (dot(RSBSMz, muz) - bAb) / (2*m);
                logpost = logpost + logdet + logquad;

        else % rare, extreme case: no snps included in the model
                betapost = [];
                logpost = 0;
        end

        loglik = logpost;

        % add the part from p(gamma|pi) * p(pi)
        logpost = logpost + (p-n_gamma) * log(1-exp(logpi)) + (n_gamma-1) * logpi;

	% add the part from p(m): normal(1, 0.1^2)
	logpost = logpost - (m-1)^2 / (2*0.01);  
end

