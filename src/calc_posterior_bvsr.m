% parent function: rss_bvsr
% main function:   calc_posterior_bvsr
% local functions: calc_posterior_bvsr_identity, calc_posterior_bvsr_general

function [betapost, loglik, logpost] = calc_posterior_bvsr(q, se, SiRiS, rank, psi, logpi, matrix_type)
% USAGE: compute p(h gamma pi|betahat) under the RSS-BVSR model

        switch matrix_type
                case 0
                        [betapost, loglik, logpost] = calc_posterior_bvsr_identity(q, se, rank, psi, logpi);
                case 1
                        [betapost, loglik, logpost] = calc_posterior_bvsr_general(q, SiRiS, rank, psi, logpi);
        end
end

function [betapost, loglik, logpost] = calc_posterior_bvsr_identity(q, se, rank, psi, logpi)
% USAGE: calculate p(h gamma pi|betahat) up to some constant when R = identity
% INPUT:
%       q: betahat ./ se^2, p by 1
%       se: standard errors of the marginal effects estimates, p by 1
%       rank: the positions of snps included in the current model, 1 by n_gamma
%       psi: the current value of the prior variance of beta, scalar
%       logpi: the current value of log(pi), scalar with range log(1/p,1)
% OUTPUT:
%       betapost: sample positive components of beta from its conditional posterior, 1 by n_gamma
%       loglik: the log likelihood p(betahat|h gamma), scalar
%       logpost: the log posterior p(h gamma pi|betahat), scalar

        logpost = 0;
        n_gamma = length(rank); 		% the number of snps included in the current model, integer
        p 	= length(q); 			% the total number of snps analyzed, integer
        Psiz 	= psi * ones(n_gamma, 1); 	% the vector of current value of psi(gamma, h), n_gamma by 1

        if n_gamma ~= 0 

                RSBSMz 	 = q(rank); 		% n_gamma by 1
                S2iz 	 = 1 ./ (se(rank).^2); 	% n_gamma by 1
                OmegaInv = S2iz + 1 ./ Psiz; 	% n_gamma by 1
                muz 	 = RSBSMz ./ OmegaInv; 	% n_gamma by 1

                % sample positive components of beta from Normal(muz, Omega)
                betapost_err = randn(n_gamma,1) ./ sqrt(OmegaInv); 		% n_gamma by 1
                betapost     = muz' + betapost_err'; 				% convert the column to row vector

                % calculate log posterior for MH acceptance ratio
                % p(h gamma pi|betahat) = Const. * p(betahat|h gamma) * p(h) * p(gamma|pi) * p(pi) 

                % first compute likelihood p(betahat|h gamma) 
                logdet  = - 0.5 * sum(log(OmegaInv)) - 0.5 * sum(log(Psiz));
                logquad = 0.5 * dot(RSBSMz, muz);
                logpost = logpost + logdet + logquad;

        else % rare, extreme case: no snps included in the current model
                betapost = [];
                logpost  = 0;
        end

        loglik = logpost;
	
	% add the part from p(gamma|pi) * p(pi)
        logpost = logpost + (p-n_gamma) * log(1-exp(logpi)) + (n_gamma-1) * logpi;
end

function [betapost, loglik, logpost] = calc_posterior_bvsr_general(q, SiRiS, rank, psi, logpi)
% USAGE: calculate p(h gamma pi|betahat) up to some constant -- general case
% INPUT:
%       q: betahat ./ se^2, p by 1
%       SiRiS: S^{-1}RS^{-1}, p by p matrix
%       rank: the marginal rank of snps included in the current model, 1 by n_gamma
%       psi: the current value of prior variance of beta, scalar
%       logpi: the current value of log(pi), scalar with range log(1/p,1)
% OUTPUT:
%       betapost: sample positive components of beta from its conditional posterior, 1 by n_gamma
%       loglik: the log likelihood p(betahat|h gamma), scalar
%       logpost: the log posterior p(h gamma pi|betahat), scalar

        logpost = 0;
        n_gamma = length(rank); 		% the number of snps in the current model, integer
        p 	= length(q); 			% the total number of snps, integer
        Psiz 	= psi * ones(n_gamma, 1); 	% the vector of current value of psi(gamma, h), n_gamma by 1

        if n_gamma ~= 0

                RSBSMz 	  = q(rank); 		% n_gamma by 1
                OmegaInvz = SiRiS(rank, rank); 	% n_gamma by n_gamma

		% add a diagonal matrix Psiz to OmegaInvz
                OmegaInvz(1:(n_gamma+1):n_gamma^2) = OmegaInvz(1:(n_gamma+1):n_gamma^2) + 1 ./ Psiz';

                L   = chol(OmegaInvz, 'lower'); 
                U   = L';
                T   = solve_tril(L, RSBSMz);
                muz = solve_triu(U, T);

                % sample positive components of beta from Normal(muz, Omega)
                betapost_err = solve_triu(U, randn(n_gamma,1));
                betapost     = muz' + betapost_err'; % convert the column to row vector

                % calculate log posterior for MH
                % p(h gamma pi|betahat) = Const. * p(betahat|h gamma) * p(h) * p(gamma|pi) * p(pi) 
                
		% first compute likelihood p(betahat| h gamma)
                logdet 	= - sum(log(diag(L))) - 0.5 * sum(log(Psiz));
                logquad = 0.5 * dot(RSBSMz, muz);
                logpost = logpost + logdet + logquad;

        else % rare, extreme case: no snps included in the model
                betapost = [];
                logpost  = 0;
        end

        loglik = logpost;

        % add the part from p(gamma|pi) * p(pi)
        logpost = logpost + (p-n_gamma) * log(1-exp(logpi)) + (n_gamma-1) * logpi;
end

