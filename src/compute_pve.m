function pve = compute_pve(beta, betahat, se, samplesize, bwd, BR, matrix_type)
% USAGE: compute PVE (approximately) using summary-level data
% INPUT:
%	beta: p by 1
%       betahat: p by 1
%       se: p by 1
%       samplesize: p by 1
%	bwd: bandwidth of R
%       BR: banded storage of R
%       matrix_type: 0 if R is identity; 1 otherwise
% OUTPUT:
%       pve: scalar

        beta = beta(:);
	betahat = betahat(:);
        se = se(:);
        samplesize = samplesize(:);
	
	p = length(beta);
        yyxx_sqrt = sqrt(samplesize .* (se.^2) + betahat.^2); 
        comp1 = beta ./ yyxx_sqrt;

        switch matrix_type
                case 1
                        %comp2 = R * comp1; % a general way of computing matrix-vector product
			comp2 = lapack('DSBMV(h, i, i, d, d, i, d, i, d, D, i)', 'L', p, bwd, 1, BR, bwd+1, comp1, 1, 0, zeros(p, 1), 1);
                case 0
                        comp2 = comp1;
        end
        pve = sum(comp1 .* comp2);
end
