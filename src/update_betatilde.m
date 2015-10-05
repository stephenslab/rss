% parent function: 	rss_bslmm
% main function: 	update_betatilde
% local functions: 	update_bt_global, compute_loglik_rssbslmm, update_bt_local, propose_bt_gibbs, update_loglik_onebt 

function [beta, bs, betatilde, loglik, lpart, qpart] = update_betatilde(q, zs, se, R, bw, BR, BM, beta, bs, zlabel, bin, betatilde, sigmabeta, lpart, qpart, matrix_type, p_gamma)
% USAGE: update betatilde given summary-level data and current value of zlabel and theta
% INPUT:
%       q: 1 by p, betahat ./ (se^2)
%       zs: 1 by p, betahat ./ se
%       se: 1 by p
%       R: p by p
%       bw: bandwidth of R
%	BR: banded storage of R
%       BM: banded storage of SiRiS
%       beta: 1 by p
%       bs: 1 by p, beta ./ se
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       betatilde: 1 by p
%       sigmabeta: 1 by 2, [sigma_poly, sqrt(sigma_poly^2+sigma_beta^2)]
%       lpart: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%	matrix_type: 0 if R is identity; 1 otherwise
%	p_gamma: 1 by p, the re-ordered pmf of 0.3*unif + 0.7*geometric based on single-snp p-value
% OUTPUT:
%       beta: 1 by p
%       bs: 1 by p, beta_new ./ se
%	betatilde: 1 by p
%       loglik: log p(betahat|betatilde_new, zlabel, theta)
%       lpart: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]

	if rand < 0.001 % restrict the number of expensive global moves
		[beta, bs, betatilde, loglik, lpart, qpart] = update_bt_global(q, se, BR, bw, BM, zlabel, bin, sigmabeta, matrix_type);
	else
		[beta, bs, betatilde, loglik, lpart, qpart] = update_bt_local(q, zs, se, R, beta, bs, zlabel, bin, betatilde, sigmabeta, lpart, qpart, matrix_type, p_gamma);
	end
end

function [beta, bs, betatilde, loglik, lpart, qpart] = update_bt_global(q, se, BR, bw, BM, zlabel, bin, sigmabeta, matrix_type)
% USAGE: update betatilde by its full conditional
% INPUT:
%       q: 1 by p, betahat ./ (se.^2)
%       se: 1 by p
%       BR: banded storage of R
%       bw: bandwidth of R
%       BM: banded storage of SiRiS
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       sigmabeta: 1 by 2, [sigma_poly, sqrt(sigma_poly^2+sigma_beta^2)]
%       matrix_type: 0 if R is identity; 1 otherwise
% OUTPUT:
%       beta: 1 by p
%       bs: 1 by p, beta ./ se
%       betatilde: 1 by p
%       loglik: log p(betahat|beta)
%       lpart: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]

        p = length(q);
        D = sigmabeta(zlabel);

        switch matrix_type
                case 0
                        var_vec    = 1 ./ ((D ./ se).^2 + 1);
                        mean_term  = var_vec(:) .* D(:) .* q(:);
                        error_term = sqrt(var_vec(:)) .* randn(p, 1);
                        betatilde  = mean_term + error_term;
                case 1
                        BOmega = BM;
                        BOmega(1, :) = BOmega(1, :) + 1 ./ (D.^2);

                        BL = lapack('DPBTRF(h, i, i, D, i, i)', 'L', p, bw, BOmega, bw+1, 0);

                        q_col = q(:);

                        Uiq = lapack('DTBTRS(h, h, h, i, i, i, d, i, D, i, i)', 'L', 'N', 'N', p, bw, 1, BL, bw+1, q_col, p, 0);
                        Udb = Uiq + randn(p, 1);
                        Dbt = lapack('DTBTRS(h, h, h, i, i, i, d, i, D, i, i)', 'L', 'T', 'N', p, bw, 1, BL, bw+1, Udb, p, 0);

                        betatilde = Dbt ./ D(:);
        end

        betatilde = betatilde'; % convert col to row vector

        beta = betatilde .* D;
        bs = beta ./ se;

	[loglik, lpart, qpart] = compute_loglik_rssbslmm(q, bw, BR, beta, bs, bin, matrix_type);
end

function [loglik, lpart, qpart] = compute_loglik_rssbslmm(q, bwd, BR, b, bs, bin, matrix_type)
% USAGE: compute log p(betahat|se, R, beta) under RSS-BSLMM model
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

function [beta_new, bs_new, betatilde_new, loglik_new, lpart_new, qpart_new] = update_bt_local(q, zs, se, R, beta, bs, zlabel, bin, betatilde, sigmabeta, lpart, qpart, matrix_type, p_gamma)
% USAGE: random scan gibbs update of zlabel, i.e., only update one snp's betatilde
% INPUT:
%       q: 1 by p, betahat ./ (se^2)
%       zs: 1 by p, betahat ./ se
%       se: 1 by p
%       R: p by p
%       beta: 1 by p
%       bs: 1 by p, beta ./ se
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       betatilde: 1 by p
%       sigmabeta: 1 by 2, [sigma_poly, sqrt(sigma_poly^2+sigma_beta^2)]
%       lpart: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       qpart: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%       matrix_type: 0 if R is identity; 1 otherwise
%       p_gamma: 1 by p, the re-ordered pmf of 0.3*unif + 0.7*geometric based on single-snp p-value
% OUTPUT:
%       beta_new: 1 by p
%       bs_new: 1 by p, beta_new ./ se
%       betatilde_new: 1 by p
%       loglik_new: log p(betahat|betatilde, zlabel_new, theta)
%       lpart_new: 1 by 2, indexed by [zlabel_new == 1, zlabel_new == 2]
%       qpart_new: 1 by 2, indexed by [both zlabel_new == 1, both zlabel_new == 2, cross term]

        % initialize quantities related to betatilde
        betatilde_new 	= betatilde;
        beta_new 	= beta;
        bs_new 		= bs;

        % randomly select ONE snp to update based on the marginal association
        j = sample(p_gamma);

        % update betatilde for the selected snp; might remain unchanged
        betatilde_new(j) = propose_bt_gibbs(j, zs, se, R, bs, zlabel, sigmabeta, matrix_type);

        % update quantities related to betatilde
        beta_new(j) = sigmabeta(zlabel(j)) * betatilde_new(j);
        bs_new(j)   = beta_new(j) / se(j);

        chg_comp = zlabel(j);
        [loglik_new, lpart_new, qpart_new] = update_loglik_onebt(j, chg_comp, q, R, bin, beta, beta_new, bs, bs_new, lpart, qpart, matrix_type);
end

function betatilde_proposed = propose_bt_gibbs(j, zs, se, R, bs, zlabel, sigmabeta, matrix_type)
% USAGE: propose a betatilde based on the full conditional: p(betatilde(j)|betahat, betatilde(-j), zlabel, theta)
% INPUT:
%       j: the jth snp to update
%       zs: 1 by p, betahat ./ se
%       se: 1 by p
%       R: p by p
%       bs: 1 by p, beta ./ se
%       zlabel: 1 by p, 1 <- polygenic and 2 <- polygenic+sparse
%       sigmabeta: 1 by 2, [sigma_poly, sqrt(sigma_poly^2+sigma_beta^2)]
%       matrix_type: 0 if R is identity; otherwise 1
% OUTPUT:
%       betatilde_proposed: the proposed value for betatilde(j)

        sig = sigmabeta(zlabel(j));

        switch matrix_type
                case 0
                        linear_tmp = zs(j);
                case 1
                        linear_tmp = zs(j) - ( bs * R(:, j) - bs(j) );
        end

        varsum = sig^2 + se(j)^2;

        sd_beta_j 	   = se(j) / sqrt(varsum);
        mu_beta_j 	   = ( sig * se(j) * linear_tmp ) / varsum;
        betatilde_proposed = mu_beta_j + sd_beta_j * randn;
end

function [loglik_new, linear_part_new, quadra_part_new] = update_loglik_onebt(j, chg_comp, q, R, bin, beta, beta_new, bs, bs_new, linear_part, quadra_part, matrix_type)
% USAGE: compute log p(betahat|betatilde, zlabel, theta) and friends when betatilde of ONE snp is updated by gibbs
% INPUT:
%       j: scalar, the jth SNP
%       chg_comp: zlabel(j)
%       q: 1 by p, betahat ./ (se^2)
%       R: p by p
%       bin: 2 by p, [(zlabel == 1) ; (zlabel == 2)]
%       beta & beta_new: 1 by p; beta(i) == beta_new(i) for all i ~= j
%       bs & bs_new: 1 by p; bs(i) == bs_new(i) for all i ~= j
%       linear_part: 1 by 2, indexed by [zlabel == 1, zlabel == 2]
%       quadra_part: 1 by 2, indexed by [both zlabel == 1, both zlabel == 2, cross term]
%       matrix_type: 0 if R is identity; 1 otherwise
% OUTPUT:
%       loglik_new: log p(betahat|betatilde, zlabel_new, theta)
%       linear_part_new: 1 by 2, indexed by [zlabel_new == 1, zlabel_new == 2]
%       quadra_part_new: 1 by 2, indexed by [both zlabel_new == 1, both zlabel_new == 2, cross term]

        % update the two components of the linear term: q' * beta
        linear_part_new = linear_part;
        linear_part_new(chg_comp) = linear_part_new(chg_comp) + q(j) * ( beta_new(j) - beta(j) );

        % update the three components of the quadratic term: bs * R * bs'
        quadra_part_new = quadra_part;
        switch matrix_type
                case 0
                        chg_sum = 0;
                        nch_sum = 0;
                case 1
                        chg_idx = bin(chg_comp, :);
                        nch_idx = bin(3-chg_comp, :);
                        chg_sum = bs(chg_idx) * R(chg_idx, j) - bs(j);
                        nch_sum = bs(nch_idx) * R(nch_idx, j);
        end
        bs_diff = bs_new(j) - bs(j);
        bs_square_diff = (bs_new(j) + bs(j)) * bs_diff;
        quadra_part_new(chg_comp) = quadra_part_new(chg_comp) + 2 * bs_diff * chg_sum + bs_square_diff;
        quadra_part_new(end) = quadra_part_new(end) + 2 * bs_diff * nch_sum;

        % update the log p(betahat|betatilde, zlabel, theta)
        loglik_new = sum(linear_part_new) - 0.5 * sum(quadra_part_new);
end
