% parent function: 	rss_bslmm
% main function: 	update_zlabel
% local functions: 	propose_z, propose_gamma, update_loglik_onez, update_loglik_twoz 
 
function [b, bs, z, bin, zc, loglik, lpart, qpart, Naccept] = update_zlabel(q, zs, se, R, b, bs, z, bin, zc, bt, logpi, sig_b, loglik, lpart, qpart, Naccept, matrix_type, p_gamma)
% USAGE: Metropolis-Hastings update of zlabel (rank-based proposal)
% SOURCE: Guan & Stephens (2011) Appendix A and Zhou et al (2013) Text S2
% INPUT:
%	q: 1 by p, betahat ./ (se.^2)
%       zs: 1 by p, betahat ./ se
%       se: 1 by p
%       R: p by p
%       b: 1 by p
%       bs: 1 by p, beta ./ se
%       z: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%	zc: 1 by 2, [sum(zlabel==1), sum(zlabel==2)]
%       bt: 1 by p
%       sig_b: 1 by 2, [sigma_poly, sqrt(sigma_poly^2+sigma_beta^2)]
%       loglik: log p(betahat|betatilde, zlabel, theta)
%       lpart: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%	Naccept: scalar
%	matrix_type: 0 if R is identity; 1 otherwise
%	p_gamma: rank-based proposal
% OUTPUT:
%       b: 1 by p
%       bs: 1 by p, beta_new ./ se
%       z: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       bin: 2 by p, [(zlabel_new == 1) ; (zlabel_new == 2)]
%       zc: 1 by 2, [sum(zlabel==1), sum(zlabel==2)]
%       loglik: log p(betahat|betatilde, zlabel_new, theta)
%       lpart: 1 by 2, indexed by [zlabel_new == 1, zlabel_new == 2]
%       qpart: 1 by 2, indexed by [both zlabel_new == 1, both zlabel_new == 2, cross term]
%	Naccept: scalar

	% compute log p(z|betahat, betatilde, pi, h, rho)
        logpos = loglik + log(1-exp(logpi))*zc(1) + logpi*zc(2);

	% small-world proposal
        if rand < 0.33
                extrastep = randperm(19);
                repeat 	  = 1 + extrastep(1);
        else
                repeat = 1;
        end

        logMHratio = 0;
	
	% compound the local proposals
	b_tmp 		= b;
        bs_tmp 		= bs;
        z_tmp 		= z;
        bin_tmp 	= bin;
	zc_tmp 		= zc;
        lpart_tmp 	= lpart;
        qpart_tmp 	= qpart;
        loglik_tmp 	= loglik;
        logpos_tmp 	= logpos;

	z_logratio = 0;
	for k=1:repeat
		[b_new, bs_new, z_new, bin_new, zc_new, loglik_new, logpos_new, lpart_new, qpart_new, z_logratio_new] ...
		= propose_z(q, zs, se, R, b_tmp, bs_tmp, z_tmp, bin_tmp, zc_tmp, bt, logpi, sig_b, loglik_tmp, logpos_tmp, lpart_tmp, qpart_tmp, matrix_type, p_gamma);

		b_tmp   	= b_new;
        	bs_tmp     	= bs_new;
        	z_tmp 		= z_new;
        	bin_tmp    	= bin_new;
		zc_tmp 		= zc_new;
        	lpart_tmp  	= lpart_new;
        	qpart_tmp  	= qpart_new;
        	loglik_tmp 	= loglik_new;
        	logpos_tmp 	= logpos_new;
		z_logratio 	= z_logratio + z_logratio_new;
	end
	
	logMHratio = logMHratio + z_logratio;
	logMHratio = logMHratio + logpos_new - logpos;

	if logMHratio>0 || log(rand)<logMHratio
                Naccept 	= Naccept + 1;
		b 		= b_new;
                bs 		= bs_new;
                z 		= z_new;
                bin 		= bin_new;
		zc 		= zc_new;
                lpart 		= lpart_new;
                qpart 		= qpart_new;
                loglik 		= loglik_new;
	end

end

function [beta_new, bs_new, zlabel_new, bin_new, labelcounts_new, loglik_new, logpos_new, lpart_new, qpart_new, z_logratio] = propose_z(q, zs, se, R, beta, bs, zlabel, bin, labelcounts, betatilde, logpi, sigmabeta, loglik, logpos, lpart, qpart, matrix_type, p_gamma)
% USAGE: Metropolis-Hastings update of zlabel using rank-based proposal (Guan & Stephens, 2011)
% INPUT:
%	q: 1 by p, betahat ./ (se.^2)
%       zs: 1 by p, betahat ./ se
%       se: 1 by p
%       R: p by p
%       beta: 1 by p
%       bs: 1 by p, beta ./ se
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       labelcounts: 1 by 2, [sum(zlabel==1), sum(zlabel==2)]
%       betatilde: 1 by p
%	logpi: scalar, log(pi)
%       sigmabeta: 1 by 2, [sigma_poly, sqrt(sigma_poly^2+sigma_beta^2)]
%       loglik: log p(betahat|betatilde, zlabel, theta)
%       logpos: log p(zlabel|betahat, betatilde, theta)
%       lpart: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%       matrix_type: 0 if R is identity; 1 otherwise
%	p_gamma: rank-based proposal 
% OUTPUT:
%       beta_new: 1 by p
%       bs_new: 1 by p, beta_new ./ se
%       zlabel_new: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       bin_new: 2 by p, [(zlabel_new == 1) ; (zlabel_new == 2)]
%       labelcounts_new: 1 by 2, [sum(zlabel_new==1), sum(zlabel_new==2)]
%       loglik_new: log p(betahat|betatilde, zlabel_new, theta)
%       logpos_new: log p(zlabel_new|betahat, betatilde, theta)
%       lpart_new: 1 by 2, indexed by [zlabel_new == 1, zlabel_new == 2]
%       qpart_new: 1 by 2, indexed by [both zlabel_new == 1, both zlabel_new == 2, cross term]
%       z_logratio: log ratio contribution of zlabel to the MH log ratio, scalar

        % initialize quantities related to zlabel
        beta_new        = beta;
        bs_new          = bs;
        zlabel_new      = zlabel;
        bin_new         = bin;
        labelcounts_new = labelcounts;

        lpart_new       = lpart;
        qpart_new       = qpart;
        loglik_new      = loglik;
        logpos_new      = logpos;

        % rank-based proposal of zlabel
        [snp_change, z_logratio] = propose_gamma(bin, p_gamma);

        % case 1: swap the zlabel of snp a and b; a and b have DIFFERENT zlabels
        if length(snp_change) == 2

                a = snp_change(1);
                b = snp_change(2);

                zlabel_new(a) = zlabel(b);
                zlabel_new(b) = zlabel(a);
                bin_new(zlabel(a), a) = 0; bin_new(zlabel(a), b) = 1;
                bin_new(zlabel(b), a) = 1; bin_new(zlabel(b), b) = 0;
                beta_new(a) = betatilde(a) * sigmabeta(zlabel_new(a));
                beta_new(b) = betatilde(b) * sigmabeta(zlabel_new(b));
                bs_new(a) = beta_new(a) / se(a);
                bs_new(b) = beta_new(b) / se(b);

                [loglik_new, lpart_new, qpart_new] = update_loglik_twoz(a, b, zs, se, R, bs, betatilde, zlabel, bin, lpart, qpart, matrix_type);

                logprior = logpos - loglik;
                % p(zlabel_new|theta) == p(zlabel|theta) b/c swapping does not change label counts
                logpos_new = loglik_new + logprior;

        end

        % case 2: flip the zlabel of snp j
        if length(snp_change) == 1

                j = snp_change;

                zlabel_new(j) = 3 - zlabel(j); % b/c 1+2=3
                bin_new(zlabel(j), j) 	  = 0;
                bin_new(zlabel_new(j), j) = 1;
                labelcounts_new(zlabel(j)) 	= labelcounts(zlabel(j)) - 1;
                labelcounts_new(zlabel_new(j)) 	= labelcounts(zlabel_new(j)) + 1;
                beta_new(j) = sigmabeta(zlabel_new(j)) * betatilde(j);
                bs_new(j)   = beta_new(j) / se(j);

                add_comp = zlabel_new(j);
                rmv_comp = zlabel(j);

                [loglik_new, lpart_new, qpart_new] = update_loglik_onez(j, add_comp, rmv_comp, q, R, beta, beta_new, bs, bs_new, bin, lpart, qpart, matrix_type);

                logprior 	= logpos - loglik;
                logprior_new 	= logprior + ( zlabel_new(j) - zlabel(j) ) * ( logpi - log(1-exp(logpi)) );
                logpos_new 	= loglik_new + logprior_new;

        end
end

function [snp_change, z_logratio] = propose_gamma(bin, p_gamma)
% USAGE: rank-based proposal of updating zlabel (Guan & Stephens, 2011)
% INPUT:
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       p_gamma: 1 by p, the pmf of the mixture of discrete uniform and geometric rv
% OUTPUT:
%       snp_change: position(s) of updated snp(s); flip or swap
%       z_logratio: scalar, log ratio contribution of zlabel to the MH log ratio

        p 	   = length(p_gamma);           % total number of snps
        z_logratio = 0;
        ngamma 	   = sum(bin(2, :));        	% total number of snps with type 2 zlabel 
        rank 	   = find(bin(2, :));         	% locations of snps with type 2 zlabel 

        % rare case 1: no snp in the current model; must add one snp
        if ngamma == 0
                r_add 	   = sample(p_gamma);                           % randomly pick one snp based on p-values
                z_logratio = z_logratio - log(p_gamma(r_add));          % compute the contribution to log MH ratio
                snp_change = r_add;
        end

        % rare case 2: all snps in the current model; must remove one snp
        if ngamma == p
                r_rmv 	   = randint(1,ngamma);                   % uniformly pick one snp to remove
                z_logratio = z_logratio + log(ngamma);          % compute the contribution to log MH ratio
                snp_change = r_rmv;
        end

        if (ngamma >= 1) && (ngamma <= p-1)
                gamma_flag = rand; % decide the type of move

                if gamma_flag < 0.4 % add a snp
                        % propose a snp NOT in the current model
                        % randomly pick a snp until the picked snp is NOT in the current model
                        while true
                                r_add = sample(p_gamma);
                                if (bin(2, r_add) == 0); break; end
                        end
                        % calculate Pr(picking a snp NOT in the current model)
                        prob_total = 1;
                        prob_total = prob_total - sum(p_gamma(rank));
                        % calculate the contribution to log MH ratio
                        z_logratio = z_logratio - log(p_gamma(r_add) / prob_total) - log(ngamma+1);
                        snp_change = r_add;

                elseif gamma_flag >=0.4 && gamma_flag < 0.8 % remove a snp
                        % uniformly draw a snp from the current model
                        col_id = randint(1,ngamma);
                        r_rmv  = rank(col_id); % return the location of the chosen snp
                        % calculate Pr(picking a snp NOT in the proposed model)
                        prob_total = 1;
                        prob_total = prob_total - sum(p_gamma(rank));
                        prob_total = prob_total + p_gamma(r_rmv); % b/c this snp is NOT included in the proposed model
                        % calculate the contribution to log MH ratio
                        z_logratio = z_logratio + log(p_gamma(r_rmv) / prob_total) + log(ngamma);
                        snp_change = r_rmv;

                else % swap inclusion/exclusion status of two snps
                        col_id = randint(1,ngamma); % pick a snp in the model to remove
                        r_rmv  = rank(col_id);
                        while true
                                r_add = sample(p_gamma); % randomly pick a snp based on marginal association rank
                                if (bin(2, r_add) == 0); break; end % until the picked snp is NOT in the current model
                        end
                        % calculate the probability of picking a snp NOT in the current model
                        prob_total = 1;
                        prob_total = prob_total - sum(p_gamma(rank));
                        z_logratio = z_logratio + log( p_gamma(r_rmv) / (prob_total + p_gamma(r_rmv) - p_gamma(r_add) ) );
                        z_logratio = z_logratio - log( p_gamma(r_add) / prob_total );
                        snp_change = [r_rmv r_add];
                end
        end
end

function [loglik_new, linear_part_new, quadra_part_new] = update_loglik_onez(j, add_comp, rmv_comp, q, R, beta, beta_new, bs, bs_new, bin, linear_part, quadra_part, matrix_type)
% USAGE: compute log p(betahat|betatilde, zlabel, theta) and friends when zlable of ONE snp is updated by gibbs 
% INPUT:
%       j: scalar, the jth SNP
%       add_comp: zlabel_new(j)
%       rmv_comp: zlabel(j)
%       q: 1 by p, betahat ./ (se^2)
%       R: p by p
%       beta & beta_new: 1 by p; beta(i) == beta_new(i) for all i ~= j
%       bs & bs_new: 1 by p; bs(i) == bs_new(i) for all i ~= j
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       linear_part: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       quadra_part: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%       matrix_type: 0 if R is identity; 1 otherwise
% OUTPUT:
%       loglik_new: log p(betahat|betatilde, zlabel_new, theta)
%       linear_part_new: 1 by 2, indexed by [zlabel_new == 1, zlabel_new == 2]
%       quadra_part_new: 1 by 2, indexed by [both zlabel_new == 1, both zlabel_new == 2, cross term]


        one_idx = bin(1, :);
        two_idx = bin(2, :);

        if add_comp == 1
                add_idx = one_idx;
                rmv_idx = two_idx;
        else
                rmv_idx = one_idx;
                add_idx = two_idx;
        end

        % update the two components of the linear term: q' * beta
        linear_part_new = linear_part;
        linear_part_new(rmv_comp) = linear_part_new(rmv_comp) - q(j) * beta(j);
        linear_part_new(add_comp) = linear_part_new(add_comp) + q(j) * beta_new(j);

        % update the three components of the quadratic term: bs * R * bs'
        quadra_part_new = quadra_part;
        switch matrix_type
                case 0
                        add_sum = 0;
                        rmv_sum = bs(j);
                case 1
                        add_sum = bs(add_idx) * R(add_idx, j);
                        rmv_sum = bs(rmv_idx) * R(rmv_idx, j);
        end
        quadra_part_new(rmv_comp) = quadra_part_new(rmv_comp) - 2 * bs(j) * rmv_sum + bs(j)^2;
        quadra_part_new(add_comp) = quadra_part_new(add_comp) + 2 * bs_new(j) * add_sum + bs_new(j)^2;
        quadra_part_new(end) = quadra_part_new(end) - 2 * bs(j) * add_sum;
        quadra_part_new(end) = quadra_part_new(end) + 2 * bs_new(j) * rmv_sum;
        quadra_part_new(end) = quadra_part_new(end) - 2 * bs(j) * bs_new(j);

        % update the log p(betahat|betatilde, zlabel, theta)
        loglik_new = sum(linear_part_new) - 0.5 * sum(quadra_part_new);
end

function [loglik_new, linear_part_new, quadra_part_new] = update_loglik_twoz(a, b, zs, se, R, bs, betatilde, zlabel, bin, linear_part, quadra_part, matrix_type)
% USAGE: compute log p(betahat|betatilde, zlabel, theta) and friends when zlable of TWO snps are swapped
% INPUT:
%       a & b : the ath and bth snp
%       zs: 1 by p, betahat ./ se
%       se: 1 by p
%       R: p by p
%       bs: 1 by p, beta ./ se
%       betatilde: 1 by p
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       linear_part: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       quadra_part: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%       matrix_type: 0 if R is identity; 1 otherwise
% OUTPUT:
%       loglik_new: log p(betahat|betatilde, zlabel_new, theta)
%       linear_part_new: 1 by 2, indexed by [zlabel_new == 1, zlabel_new == 2]
%       quadra_part_new: 1 by 2, indexed by [both zlabel_new == 1, both zlabel_new == 2, cross term]


        one_idx = bin(1, :);
        two_idx = bin(2, :);

        ba_ratio = (betatilde(b)/betatilde(a)) * (se(a)/se(b));
        ab_ratio = 1 / ba_ratio;

        % update the two components of the linear term: q' * beta
        linear_part_new = linear_part;
        linear_part_new(zlabel(a)) = linear_part_new(zlabel(a)) + bs(a) * (zs(b)*ba_ratio-zs(a));
        linear_part_new(zlabel(b)) = linear_part_new(zlabel(b)) + bs(b) * (zs(a)*ab_ratio-zs(b));

        % update the three components of the quadratic term: bs * R * bs'
        quadra_part_new = quadra_part;

        if zlabel(a) == 1
                other_a = one_idx; other_a(a)=0;
                other_b = two_idx; other_b(b)=0;
        else
                other_a = two_idx; other_a(a)=0;
                other_b = one_idx; other_b(b)=0;
        end

        switch matrix_type
                case 0
                        sum_aa = 0;
                        sum_ab = 0;
                        sum_ba = 0;
                        sum_bb = 0;
                case 1
                        sum_aa = bs(other_a) * R(other_a, a);
                        sum_ab = bs(other_a) * R(other_a, b);
                        sum_ba = bs(other_b) * R(other_b, a);
                        sum_bb = bs(other_b) * R(other_b, b);
        end

        quadra_part_new(zlabel(a)) = quadra_part_new(zlabel(a)) + bs(a) * 2 * (ba_ratio*sum_ab - sum_aa) + (ba_ratio^2-1) * bs(a)^2;
        quadra_part_new(zlabel(b)) = quadra_part_new(zlabel(b)) + bs(b) * 2 * (ab_ratio*sum_ba - sum_bb) + (ab_ratio^2-1) * bs(b)^2;
        quadra_part_new(end) = quadra_part_new(end) + 2 * bs(a) * (ba_ratio*sum_bb - sum_ba);
        quadra_part_new(end) = quadra_part_new(end) + 2 * bs(b) * (ab_ratio*sum_aa - sum_ab);

        % update log p(zlabel|betahat, betatilde, theta)
        loglik_new = sum(linear_part_new) - 0.5 * sum(quadra_part_new);
end

