function [bsam, zsam, lpsam, hsam, rsam, Naccept] = rss_bslmm(betahat, se, R, Nsnp, Ndraw, Nburn, Nthin)
% USAGE: simulate posterior draws under the RSS model with BSLMM prior
% INPUT:
%	betahat: 1 by p, effect size estimates under single-SNP model
%	se: 1 by p, standard errors of betahat
%	R: p by p, population LD, symmetric, positive definite matrix
%	Nsnp: 1 by p, the number of individuals for each SNP
%	Ndraw: scalar, the total number of MCMC draws
%	Nburn: scalar, the length of burn-in
%	Nthin: scalar, the length of thinning interval
% OUTPUT:
%	bsam: Nsam by p, the output chain of beta
%	zsam: Nsam by p, the output chain of zlabel (0/1) mapping to (1/2)
%	lpsam: Nsam by 1, the output chain of log(pi)
%	hsam: Nsam by 1, the output chain of h
%	rsam: Nsam by 1, the output chain of rho
% 	Naccept: 4 by 1, the total number of acceptance in four Metropolis-Hastings (MH) steps

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

	Nstart 	= Nburn + Nthin;
        Nsam 	= length(Nstart:Nthin:Ndraw);
	tsam 	= 0;

	z_Naccept = 0;
	p_Naccept = 0;
	h_Naccept = 0;
	r_Naccept = 0;
	stepsize  = 1; % sd of gaussian random walk

	% preallocate memory for mcmc output
	zsam 	= zeros(Nsam, p);
	bsam 	= zeros(Nsam, p);
	lpsam 	= zeros(Nsam, 1);
	hsam 	= zeros(Nsam, 1);
	rsam 	= zeros(Nsam, 1);

	% pre-compute some quantities used in mcmc iterations
	q 	= betahat ./ (se.^2);
	zs 	= betahat ./ se;
	xxyysum = sum( 1 ./ ( Nsnp .* (se.^2) ) );  % used in induced prior SDs on beta

	% save SiRiS in banded storage to exploit benefit of lapack
	Si 	= 1 ./ se(:);
	SiRiS 	= repmat(Si, 1, p) .* R .* repmat(Si', p, 1);
	% bandwidth of SiRiS
	bw = findbandwidth(SiRiS);
	% band storage of SiRiS
	BM = bandstorage(SiRiS, bw);
	% band storage of R
	BR = bandstorage(R, bw);
	fprintf('banded storage done \n')

	% initiate the latent label based on marginal z-scores
	abz = abs(zs);
	[gamma_start, ngamma_start, snp_rank] = initiate_model(p, abz);
	z = gamma_start + 1;
	% initiate the label counts given z
	zc 	= zeros(1, 2); 
	zc(1) 	= sum(z==1); 
	zc(2) 	= p - zc(1);
	% record the binary status of two label types
	bin = [(z==1); (z==2)];
	% initiate betatilde by standard normal
	bt = randn(1, p);
	% initiate pi based on marginal p-values
	lp = log(ngamma_start / p);
	% initiate (h, rho) by unif(0,1)
	h = rand;
	r = rand;
	% initiate beta given betatilde, zlabel, pi, h and rho
	sig_b 	= compute_sd_bslmm(xxyysum, lp, h, r);
	b 	= sig_b(z) .* bt;	
	% compute the initial value of beta ./ se
	bs = b ./ se;
	fprintf('model initilization done \n')

	% rank-based proposal distribution
	p_gamma = pmf_ugmix(p, 2e3);
	p_gamma = p_gamma(snp_rank);

	% if R is an identity matrix, then we have a faster implementation
	matrix_type = 1;
	if R == eye(p)
		matrix_type = 0;
		fprintf('faster computation avaialbe if R is identity matrix \n')
	end

	% run the FIRST mcmc step outside for loop
	[loglik, lpart, qpart] = compute_loglik_rssbslmm(q, bw, BR, b, bs, bin, matrix_type);
	
	% viz the mcmc progress
	progress_bar = progress('init','start MCMC iteration');
	
	% mcmc rountine to sample (betatilde, zlabel, pi, h, rho)
	for t=1:Ndraw

		% update zlabel given betatilde, pi, h, rho
		[b, bs, z, bin, zc, ~, lpart, qpart, z_Naccept] = update_zlabel(q, zs, se, R, b, bs, z, bin, zc, bt, lp, sig_b, loglik, lpart, qpart, z_Naccept, matrix_type, p_gamma);

		% update betatilde given zlabel, pi, h, rho
		[b, bs, bt, loglik, lpart, qpart] = update_betatilde(q, zs, se, R, bw, BR, BM, b, bs, z, bin, bt, sig_b, lpart, qpart, matrix_type, p_gamma);

		% update pi given betatilde, zlabel, h, rho
		[b, bs, lp, sig_b, loglik, lpart, qpart, p_Naccept] = update_pi(xxyysum, b, bs, z, zc, lp, h, r, sig_b, loglik, lpart, qpart, p, p_Naccept, stepsize);
		
		% update h given betatilde, zlabel, pi, rho
		[b, bs, h, sig_b, loglik, lpart, qpart, h_Naccept] = update_h(xxyysum, b, bs, z, lp, h, r, sig_b, loglik, lpart, qpart, h_Naccept, stepsize);
		
		% update rho given betatilde, zlabel, pi, h
		[b, bs, r, sig_b, loglik, lpart, qpart, r_Naccept] = update_rho(xxyysum, b, bs, z, lp, h, r, sig_b, loglik, lpart, qpart, r_Naccept, stepsize);

		% record the MCMC simulated draws of parameters
		if t > Nburn && mod(t-Nburn, Nthin) == 0
			tsam = tsam + 1;

			bsam(tsam, :) 	= b;
			zsam(tsam, :) 	= z;
			lpsam(tsam) 	= lp;
			hsam(tsam) 	= h;
			rsam(tsam) 	= r;
		end
		
		progress_bar = progress(progress_bar, t/Ndraw);
	end
	
	Naccept = [z_Naccept p_Naccept h_Naccept r_Naccept];
end

function bandwidth = findbandwidth(A)
% USAGE: find the bandwidth of a symmetric band matrix
% SOURCE: Numerical Computing with MATLAB: Revised Reprint pp 72
        [i, j] = find(A);
        bandwidth = max(abs(i-j));
end

function B = bandstorage(A, p)
% USAGE: create the banded storage of a symmetric banded matrix A
% SOURCE: http://www.mathworks.com/matlabcentral/fileexchange
%         /31131-efficient-cholesky-decomposition-of-symmetric-banded-matrix 
% INPUT: upper or lower bandwidth p and a symmetric matrix A
% OUTPUT: B, compressed form of A

        dim=size(A);
        if ~(dim(1)==dim(2))
                error('A must be square')
        end
        if (all((all(A)~=all(A'))))
                error('A must be symmetric')
        end

        n = dim(1);
        B = zeros(p+1, n);
        for i=1:n % column
                if i<=n-p
                        for j=i:p+i
                                B(j-i+1, i) = A(j, i);
                        end
                else
                        for j=i:n
                                B(j-i+1, i) = A(j, i);
                        end
                end
        end
end

function [gamma_start, ngamma_start, snp_rank] = initiate_model(p, abz)
% USAGE: assign initial value to gamma for MCMC run by marginal association p-value/z-score
% INPUT:
%       p: the number of snps analyzed
%	abz: the absolute z-score obtained from single-marker model
% OUTPUT:
%       gamma_start: the initial state of the inclusion vector, 0/1, 1 by p
% 	ngamma_start: integer, the initial number of snps included in the model
%	snp_rank: the rank of each snp based on its marginal p-value, integer, 1 by p

        gamma_start 	= zeros(1, p);
        q_genome 	= invnormcdf( 1 - (0.025/p) ); % or use matlab stat tool box: norminv 
        in_loci 	= (abz > q_genome);
	ngamma_start 	= sum(in_loci);
	[~, snp_rank] 	= sort(abz, 'descend');
	baseline 	= min(10, p-1); % avoid that initial pi == 1 (-Inf log likelihood)

	if ngamma_start < baseline
		ngamma_start = baseline;
		gamma_start(snp_rank(1:baseline)) = 1;
	else
		gamma_start(in_loci) = 1;
	end
end 

function sig_b = compute_sd_bslmm(xxyysum, lp, h, r)
% USAGE: compute the induced SDs in BSLMM prior
% INPUT:
%       xxyysum: scalar
%       lp: scalar, current value of log(pi)
%       h: scalar, current value of h
%       r: scalar, current value of rho
% OUTPUT:
%       sig_b: 1 by 2

        sig2_beta = (h*r) / (exp(lp) * xxyysum);
        sig2_poly = (h*(1-r)) / xxyysum;

        sig2_b = [sig2_poly sig2_poly+sig2_beta];
        sig_b = sqrt(sig2_b);
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

function [loglik, lpart, qpart] = compute_loglik_rssbslmm(q, bwd, BR, b, bs, bin, matrix_type)
% USAGE: compute log p(betahat|se, R, beta) under bslmm prior
% INPUT:
%       q: 1 by p, beta ./ (se.^2)
%	bwd: bandwidth of R matrix
%       BR: banded storage of estimated LD matrix
%       b: 1 by p, beta
%       bs: 1 by p, beta ./ se
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       matrix_type: 0 if R is identity; 1 otherwise
% OUTPUT:
%       loglik: log p(betahat|se, R, beta)
%       lpart: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]

        lpart = zeros(1, 2);
        qpart = zeros(1, 3);

        lpart(1) = q(bin(1, :)) * b(bin(1, :))';
        lpart(2) = q(bin(2, :)) * b(bin(2, :))';

        switch matrix_type
                case 0
                        qpart(1) = bs(bin(1, :)) * bs(bin(1, :))';
                        qpart(2) = bs(bin(2, :)) * bs(bin(2, :))';
                case 1
			% direct implementation
                        %qpart(1) = bs(bin(1, :)) * R(bin(1, :), bin(1, :)) * bs(bin(1, :))';
                        %qpart(2) = bs(bin(2, :)) * R(bin(2, :), bin(2, :)) * bs(bin(2, :))';
                        %qpart(3) = bs * R * bs' - qpart(1) - qpart(2);
			
			p = length(q);
			y = lapack('DSBMV(h, i, i, d, d, i, d, i, d, D, i)', 'L', p, bwd, 1, BR, bwd+1, bs', 1, 0, zeros(p, 1), 1);
			qpart_all = bs * y;
			
			bs_1 = bs .* double(bin(1, :));
			y_1 = lapack('DSBMV(h, i, i, d, d, i, d, i, d, D, i)', 'L', p, bwd, 1, BR, bwd+1, bs_1', 1, 0, zeros(p, 1), 1);
			qpart(1) = bs_1 * y_1;

			bs_2 = bs .* double(bin(2, :));
                        y_2 = lapack('DSBMV(h, i, i, d, d, i, d, i, d, D, i)', 'L', p, bwd, 1, BR, bwd+1, bs_2', 1, 0, zeros(p, 1), 1);
                        qpart(2) = bs_2 * y_2;

			qpart(3) = qpart_all - qpart(1) - qpart(2);
        end

        loglik = sum(lpart) - 0.5*sum(qpart);
end

function [beta, bs, logpi, sigmabeta, loglik, lpart, qpart, Naccept] = update_pi(xxyysum, beta_old, bs_old, zlabel, labelcounts, logpi_old, h_old, rho_old, sigmabeta_old, loglik_old, lpart_old, qpart_old, p, Naccept, stepsize)
% USAGE: Metropolis-Hastings update of pi given zlabel and betatilde
% INPUT:
%       xxyysum: scalar, used in inducing sigma_beta and sigma_polygenic
%       beta_old: 1 by p
%       bs_old: 1 by p, beta_old ./ se
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       labelcounts: 1 by 2, [sum(zlabel==1) sum(zlabel==2)]
%       logpi_old: scalar, current value of log(pi)
%       h_old: scalar
%       rho_old: scalar
%       sigmabeta_old: 1 by 2, current value of [sigma_polygenic, sqrt(sigma_beta^2+sigma_polygenic)]
%       loglik_old: scalar, log p(betahat|betatilde, zlabel, theta_old)
%       lpart_old: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart_old: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%       p: integer, the total number of variants
%       Naccept: integer
%       stepsize: scalar, the sd of normal random walk
% OUTPUT:
%       beta: 1 by p
%       bs: 1 by p, beta ./ se
%       logpi: scalar
%       sigmabeta: 1 by 2
%       loglik: scalar
%       lpart: 1 by 2
%       qpart: 1 by 3
%       Naccept: integer

        %% STEP 0: compute log p(pi|betahat, betatilde, zlabel, h, rho)
        logpos_old = loglik_old + log(1-exp(logpi_old))*labelcounts(1) + logpi_old*labelcounts(2) - logpi_old;
        logMHratio = 0;

        %% STEP 1: propose new values for pi
        if rand < 0.33
                extrastep = randperm(19);
                repeat = 1 + extrastep(1);
        else
                repeat = 1;
        end

        % symmetric gaussian random walk proposal for logit transformed log(pi)
        [logpi_new, pi_logratio] = propose_pi(logpi_old, stepsize, repeat, p);
        logMHratio = logMHratio + pi_logratio;

        %% STEP 2: calculate proposed sigmabeta and mh acceptance ratio
        sigmabeta_new = compute_sd_bslmm(xxyysum, logpi_new, h_old, rho_old);
        [loglik_new, lpart_new, qpart_new] = update_loglik_sigma(sigmabeta_old, sigmabeta_new, lpart_old, qpart_old);
        logpos_new = loglik_new + dot([log(1-exp(logpi_new)) logpi_new], labelcounts) - logpi_new;
        logMHratio = logMHratio + logpos_new - logpos_old;

        %% STEP 3: Metropolis-Hastings update
        if logMHratio>0 || log(rand)<logMHratio
                Naccept 	= Naccept + 1;
                loglik 		= loglik_new;
                lpart 		= lpart_new;
                qpart 		= qpart_new;
                logpi 		= logpi_new;
                sigmabeta 	= sigmabeta_new;
                linear_scalar 	= sigmabeta_new ./ sigmabeta_old;
                beta 		= linear_scalar(zlabel) .* beta_old;
                bs 		= linear_scalar(zlabel) .* bs_old;
        else
                loglik 		= loglik_old;
                lpart 		= lpart_old;
                qpart 		= qpart_old;
                logpi 		= logpi_old;
                sigmabeta 	= sigmabeta_old;
                beta 		= beta_old;
                bs 		= bs_old;
        end
end

function [beta, bs, h, sigmabeta, loglik, lpart, qpart, Naccept] = update_h(xxyysum, beta_old, bs_old, zlabel, logpi_old, h_old, rho_old, sigmabeta_old, loglik_old, lpart_old, qpart_old, Naccept, stepsize)
% USAGE: Metropolis-Hastings update of h given zlabel and betatilde
% INPUT:
%       xxyysum: scalar, used in inducing sigma_beta and sigma_polygenic
%       beta_old: 1 by p
%       bs_old: 1 by p, beta_old ./ se
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       logpi_old: scalar, current value of log(pi)
%       h_old: scalar
%       rho_old: scalar
%       sigmabeta_old: 1 by 2, current value of [sigma_polygenic, sqrt(sigma_beta^2+sigma_polygenic)]
%       loglik_old: scalar, log p(betahat|betatilde, zlabel, theta_old)
%       lpart_old: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart_old: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%       Naccept: integer
%       stepsize: scalar, the sd of normal random walk
% OUTPUT:
%       beta: 1 by p
%       bs: 1 by p, beta ./ se
%       h: scalar
%       sigmabeta: 1 by 2
%       loglik: scalar
%       lpart: 1 by 2
%       qpart: 1 by 3
%       Naccept: integer

        %% STEP 0: compute log p(h|betahat, betatilde, zlabel, pi, rho)
        logpos_old = loglik_old;
        logMHratio = 0;

        %% STEP 1: propose new values for h
        if rand < 0.33
                extrastep = randperm(19);
                repeat = 1 + extrastep(1);
        else
                repeat = 1;
        end

        % symmetric gaussian random walk proposal for log(h/(1-h))
        [h_new, h_logratio] = propose_h(h_old, stepsize, repeat);
        logMHratio = logMHratio + h_logratio;

        %% STEP 2: calculate proposed sigmabeta and mh acceptance ratio
        sigmabeta_new = compute_sd_bslmm(xxyysum, logpi_old, h_new, rho_old);
        [loglik_new, lpart_new, qpart_new] = update_loglik_sigma(sigmabeta_old, sigmabeta_new, lpart_old, qpart_old);
        logpos_new = loglik_new;
        logMHratio = logMHratio + logpos_new - logpos_old;

        %% STEP 3: metropolis-hastings update
        if logMHratio>0 || log(rand)<logMHratio
                Naccept 	= Naccept + 1;
                loglik 		= loglik_new;
                lpart 		= lpart_new;
                qpart 		= qpart_new;
                h 		= h_new;
                sigmabeta 	= sigmabeta_new;
                linear_scalar 	= sigmabeta_new ./ sigmabeta_old;
                beta 		= linear_scalar(zlabel) .* beta_old;
                bs 		= linear_scalar(zlabel) .* bs_old;
        else
                loglik 	  = loglik_old;
                lpart 	  = lpart_old;
                qpart 	  = qpart_old;
                h 	  = h_old;
                sigmabeta = sigmabeta_old;
                beta 	  = beta_old;
                bs 	  = bs_old;
        end
end

function [beta, bs, rho, sigmabeta, loglik, lpart, qpart, Naccept] = update_rho(xxyysum, beta_old, bs_old, zlabel, logpi_old, h_old, rho_old, sigmabeta_old, loglik_old, lpart_old, qpart_old, Naccept, stepsize)
% USAGE: Metropolis-Hastings update of rho given zlabel and betatilde
% INPUT:
%       xxyysum: scalar, used in inducing sigma_beta and sigma_polygenic
%       beta_old: 1 by p
%       bs_old: 1 by p, beta_old ./ se
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       logpi_old: scalar, current value of log(pi)
%       h_old: scalar
%       rho_old: scalar
%       sigmabeta_old: 1 by 2, current value of [sigma_polygenic, sqrt(sigma_beta^2+sigma_polygenic)]
%       loglik_old: scalar, log p(betahat|betatilde, zlabel, theta_old)
%       lpart_old: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart_old: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%       Naccept: integer
%       stepsize: scalar, the sd of normal random walk
% OUTPUT:
%       beta: 1 by p
%       bs: 1 by p, beta ./ se
%       rho: scalar
%       sigmabeta: 1 by 2
%       loglik: scalar
%       lpart: 1 by 2
%       qpart: 1 by 3
%       Naccept: integer

        %% STEP 0: compute log p(rho|betahat, betatilde, zlabel, pi, h)
        logpos_old = loglik_old;

        logMHratio = 0;

        %% STEP 1: propose new values for rho
        if rand < 0.33
                extrastep = randperm(19);
                repeat = 1 + extrastep(1);
        else
                repeat = 1;
        end

        % symmetric gaussian random walk proposal for log(rho/(1-rho))
        [rho_new, rho_logratio] = propose_h(rho_old, stepsize, repeat);

        logMHratio = logMHratio + rho_logratio;

        %% STEP 2: calculate proposed sigmabeta and mh acceptance ratio
        sigmabeta_new = compute_sd_bslmm(xxyysum, logpi_old, h_old, rho_new);
        [loglik_new, lpart_new, qpart_new] = update_loglik_sigma(sigmabeta_old, sigmabeta_new, lpart_old, qpart_old);
        logpos_new = loglik_new;
        logMHratio = logMHratio + logpos_new - logpos_old;

        %% STEP 3: metropolis-hastings update
        if logMHratio>0 || log(rand)<logMHratio
                Naccept 	= Naccept + 1;
                loglik 		= loglik_new;
                lpart 		= lpart_new;
                qpart 		= qpart_new;
                rho 		= rho_new;
                sigmabeta 	= sigmabeta_new;
                linear_scalar 	= sigmabeta_new ./ sigmabeta_old;
                beta 		= linear_scalar(zlabel) .* beta_old;
                bs 		= linear_scalar(zlabel) .* bs_old;
        else
                loglik 		= loglik_old;
                lpart 		= lpart_old;
                qpart 		= qpart_old;
                rho 		= rho_old;
                sigmabeta 	= sigmabeta_old;
                beta 		= beta_old;
                bs 		= bs_old;
        end
end

function [logpi_new, pi_logratio] = propose_pi(logpi_old, sigma_gauss, repeat, p)
% USAGE: symmetric gaussian random walk proposal for logit transformed log(pi)
% INPUT:
%       logpi_old: the current value of log(pi), scalar
%       sigma_gauss: SD of the gaussian proposal
%       repeat: small world proposal
%       p: the total number of SNPs
% OUTPUT:
%       logpi_new: the proposed value of log(pi), scalar
%       pi_logratio: the contribution of proposal for log(pi) to log MH acceptance ratio

        pi_logratio = 0;

        theta_old = (log(p)+logpi_old) / log(p);
        kappa_old = log(theta_old) - log(1-theta_old);

        for i=1:repeat
                kappa_new = kappa_old + sigma_gauss * randn;
                theta_new = exp(kappa_new) / ( exp(kappa_new) + 1 );
                logpi_new = (theta_new - 1) * log(p);

                pi_logratio = pi_logratio + log(theta_new) + log(1-theta_new);
                pi_logratio = pi_logratio - log(theta_old) - log(1-theta_old);
                pi_logratio = pi_logratio + logpi_new;
                pi_logratio = pi_logratio - logpi_old;

                logpi_old = logpi_new;
                theta_old = theta_new;
                kappa_old = kappa_new;
        end
end

function [h_new, h_logratio] = propose_h(h_old, sigma_gauss, repeat)
% USAGE: symmetric gaussian random walk proposal for logit(h)
%               logit(h_new) = logit(h_old) + Normal(0, sigma_gauss^2)
% INPUT:
%       h_old: the current value of h, scalar
%       sigma_gauss: SD of the gaussian proposal
%       repeat: small-world proposal, scalar
% OUTPUT:
%       h_new: the proposed value of h, scalar
%       h_logratio: - log J(h_new|h_old) + log J(h_old|h_new)

        h_logratio = 0;
        lh_old 	   = log(h_old / (1-h_old));

        for i=1:repeat
                lh_new = lh_old + randn * sigma_gauss;
                h_new  = exp(lh_new) / (exp(lh_new)+1);

                h_logratio = h_logratio + log(h_new) + log(1-h_new);
                h_logratio = h_logratio - log(h_old) - log(1-h_old);

                lh_old = lh_new;
                h_old  = h_new;
        end
end

function [loglik_new, linear_part_new, quadra_part_new] = update_loglik_sigma(sigmabeta, sigmabeta_new, linear_part, quadra_part)
% USAGE: compute log p(betahat|betatilde, zlabel, theta) and friends when hyper theta is updated by metropolis-hastings
% INPUT:
%       sigmabeta: 1 by 2
%       sigmabeta_new: 1 by 2
%       linear_part: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       quadra_part: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
% OUTPUT:
%       loglik_new: log p(betahat|betatilde, zlabel, theta_new)
%       linear_part_new: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       quadra_part_new: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]

        % update the two components of the linear term: q' * beta
        linear_scalar   = sigmabeta_new ./ sigmabeta;
        linear_part_new = linear_part .* linear_scalar;

        % update the three components of the quadratic term: bs * R * bs'
        quadra_scalar   = [linear_scalar(1)^2 linear_scalar(2)^2 linear_scalar(1)*linear_scalar(2)];
        quadra_part_new = quadra_part .* quadra_scalar;

        % update the log p(betahat|betatilde, zlabel, theta)
        loglik_new = sum(linear_part_new) - 0.5 * sum(quadra_part_new);
end


