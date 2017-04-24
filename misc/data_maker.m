% main function: data_maker
% local functions: XY_maker, effectsize_maker, usepve, single_linreg 

function [true_para, individual_data, summary_data] = data_maker(X, betatype, pve, myseed)
% USAGE: simulate individual-level gwas data and compute the single-snp summary-level data
% INPUT:
%	X: n by p genotype matrix
%	betatype: number of causal effects, [num_large, num_small]
%	pve: scalar, user-defined pve
%	myseed: integer, random seed used in data generation	
% OUTPUT:
%	true_para: true pve, beta, gamma and sigma 
%	individual_data: centered phenotype y and genotype X
%	summary_data: single-snp analysis betahat and se 

	%% simulate the individual-level data of gwas
	[y, X, beta, gamma, Nsnp, sigma] = XY_maker(X, betatype, pve, myseed);

	%% compute the single-snp summary stats
	[betahat, se] = single_linreg(y, X);

	%% output
	true_para 	= {pve, beta, gamma, sigma};
	individual_data = {y, X};
	summary_data 	= {betahat, se, Nsnp};
end

function [y, X, beta, gamma, Nsnp, sigma] = XY_maker(X, betatype, pve, myseed)
% USAGE: generate continuous phenotype under additive model
% INPUT:
%	X: n by p genotype matrix
%	betatype: [num_large, num_small]
%	pve: scalar, user-defined pve
%	myseed: integer, random seed used in data generation
% OUTPUT:
%	y: n by 1 centered trait vector
%	X: n by p column-centered genotype matrix
%	beta: p by 1, true multi-snp genetic effect
%	gamma: p by 1, causal indicator for each snp
%	Nsnp: p by 1, sample size for each snp
%	sigma: residual SD to obtain the given pve under model Y=XB+E

	rng(myseed, 'twister');
	[n, p] 		= size(X);
	Nsnp 		= n * ones(p, 1);	
	[beta, gamma] 	= effectsize_maker(p, betatype); 
        X 		= X - repmat(mean(X),n,1); 		% center columns of genotypes
	sigma 		= usepve(X, beta, n, pve);		% decide residual sd based on pve
	y 		= X * beta + sigma * randn(n,1);		
	y 		= y - mean(y); 				% center trait
end

function [beta, gamma] = effectsize_maker(p, betatype)
% USAGE: generate effect sizes for various scenarios
% INPUT:
%	p: scalar, the total number of snps analyzed
%	betatype: [num_large, num_small]
% OUTPUT:
%       beta: p by 1, true multi-snp genetic effect
%       gamma: p by 1, causal indicator for each snp (0.01 for small effect)

	num_large = betatype(1);
	num_small = betatype(2);

	beta  = zeros(p, 1);
	gamma = zeros(p, 1);

	II 	 = randperm(p); 
	I 	 = II(1:num_large);
	beta(I)  = randn(num_large, 1);
	gamma(I) = 1;
	
	if num_small ~= 0
		LL 	  = II((num_large+1):end);
		JJ 	  = LL(randperm(p-num_large)); 
		num_small = min(p-num_large, num_small);
		J 	  = JJ(1:num_small);
		beta(J)   = 0.001 * randn(num_small, 1);
		gamma(J)  = 0.001;
	end

	if (sum(gamma == 1) ~= num_large) || (sum(gamma == 0.001) ~= num_small)
		error('check the generation of effect sizes')
	end
end

function sigma = usepve(X, beta, n, pve)
% USAGE: find residual SD to obtain given pve under model Y=XB+E
% INPUT:
%	X: n by p genotype matrix
%	beta: p by 1 genetic effect
%	n: sample size
%	pve: prespecified pve value
% OUTPUT:
%	sigma: residual SD

	part_1 = (1-pve) / pve;
	xb     = X * beta;  
	part_2 = dot(xb, xb) / n;
	sigma  = sqrt(part_1 * part_2);
end

function [betahat, se] = single_linreg(y, X)
% USEAGE: variable-by-variable regression; n obs and p vars;
% INPUT:	
%	y: n*1 vector; response; 
%	X: n*p matrix; design matrix; 
% OUTPUT: 	
%	betahat: p*1 vector; marginal effect estimate;
%	se: p*1 vector; mle std error of marginal effect;

% NB: y must be centered and X must be column centered; THIS IS IMPORTANT
	[n, p] 	 = size(X);
	SX2 	 = sum(X .* X); 
	betahat  = (X'* y)./ SX2';
	yrep 	 = repmat(y, 1, p);
	brep 	 = repmat(betahat', n, 1);
	resmat 	 = yrep - X .* brep;
	sigmahat = sqrt( sum(resmat .* resmat) / n);
	se 	 = (SX2').^(-0.5) .* sigmahat';
	betahat  = betahat(:);
	se 	 = se(:);
end

