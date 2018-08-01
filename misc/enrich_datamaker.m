% main function: enrich_data_maker
% local functions: XY_maker, effectsize_maker, usepve, single_linreg 

function [true_para, individual_data, summary_data] = enrich_data_maker(X, thetatype, pve, myseed, path_snps)
% USAGE: simulate individual-level gwas data and compute the single-snp summary-level data
% INPUT:
%	X: n by p genotype matrix
%	thetatype: number of causal effects, [theta0, theta]
%	pve: scalar, user-defined pve
%	myseed: integer, random seed used in data generation
%       path_snps: |path_snps| by 1, index of snps that are in a pathway	
% OUTPUT:
%	true_para: true pve, beta, gamma and sigma 
%	individual_data: centered phenotype y and genotype X
%	summary_data: single-snp analysis betahat and se 

	%% simulate the individual-level data of gwas
	[y, X, beta, gamma, Nsnp, sigma] = XY_maker(X, thetatype, pve, myseed, path_snps);

	%% compute the single-snp summary stats
	[betahat, se] = single_linreg(y, X);

	%% output
	true_para 	= {pve, beta, gamma, sigma};
	individual_data = {y, X};
	summary_data 	= {betahat, se, Nsnp};
end

function [y, X, beta, gamma, Nsnp, sigma] = XY_maker(X, thetatype, pve, myseed, path_snps)
% USAGE: generate continuous phenotype under additive model 
% INPUT:
%	X: n by p genotype matrix
%	thetatype: [theta0, theta]
%	pve: scalar, user-defined pve
%	myseed: integer, random seed used in data generation
%       path_snps: |path_snps| by 1, index of snps that are in a pathway
% OUTPUT:
%	y: n by 1 centered trait vector
%	X: n by p column-centered genotype matrix
%	beta: p by 1, true multi-snp genetic effect
%	gamma: p by 1, causal indicator for each snp
%	Nsnp: p by 1, sample size for each snp
%	sigma: residual SD to obtain the given pve under model Y=XB+E

	rng(myseed, 'twister');
	[n, p] 	= size(X);
	Nsnp 	= n * ones(p, 1);
	
	[beta, gamma] = effectsize_maker(p, thetatype, path_snps);
 
	X 	= X - repmat(mean(X),n,1); 		% center columns of genotypes
	sigma 	= usepve(X, beta, n, pve);		% decide residual sd based on pve
	y 	= X * beta + sigma * randn(n,1);		
	y 	= y - mean(y); 				% center trait
end

function [beta, gamma] = effectsize_maker(p, thetatype, path_snps)
% USAGE: generate effect sizes under the null and enrichment hypothesis
%	 the variance of causal effects is fixed as 1
% INPUT:
%	p: scalar, the total number of snps analyzed
%	thetatype: [theta0, theta]
%	path_snps: |path_snps| by 1, index of snps that are in a pathway
% OUTPUT:
%       beta: p by 1, true multi-snp genetic effect
%       gamma: p by 1, causal indicator for each snp (0/1)

	theta0 = thetatype(1);
	theta  = thetatype(2);

	beta  = zeros(p, 1);
	gamma = zeros(p, 1);

	if theta == 0
		%% effect sizes under the null hypothesis
		pival = 1 / (1+10^(-theta0));
		gamma = (rand(p, 1) < pival);
		beta  = randn(p, 1) .* gamma;
	
	else
		%% effect sizes under the enrichment hypothesis
		% annotate snps based on pathway assignment
		inpath 	= path_snps; 
		outpath = setdiff((1:p)', path_snps);

		% simulate effect sizes of snps outside pathway
		p_out 		= length(outpath);
		out_pival 	= 1 / (1+10^(-theta0)); 
		gamma(outpath) 	= (rand(p_out, 1) < out_pival);
		beta(outpath) 	= randn(p_out, 1) .* gamma(outpath); 

		% simulate effect sizes of snps inside pathway
		p_in 		= length(inpath);
		in_pival 	= 1 / (1+10^(-theta0-theta));
		gamma(inpath) 	= (rand(p_in, 1) < in_pival);
		beta(inpath) 	= randn(p_in, 1) .* gamma(inpath);

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

	part_1 	= (1-pve) / pve;
	xb 	= X * beta;  
	part_2 	= dot(xb, xb) / n;
	sigma 	= sqrt(part_1 * part_2);
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

