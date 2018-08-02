clear;

% add codes to search path
addpath('../../src_vb/');

% load summary-level data
example_data = matfile('example5_simulated_data.mat');

betahat = example_data.betahat;
se      = example_data.se;
R       = example_data.R;
snps    = example_data.snps;

p     = length(betahat);
Si    = 1 ./ se(:);
SiRiS = sparse(repmat(Si, 1, p) .* R .* repmat(Si', p, 1));
clear Si R;

% specify hyper-parameters
theta0 = (-4.5:0.05:-3.5)';      % grid for the genome-wide log-odds (base 10)
theta  = (1.5:0.05:2.5)';        % grid for the log-fold enrichment (base 10)
sigb   = 1;                      % prior SD of genetic effects

% initialize the variational parameters
myseed = 200;

rng(myseed, 'twister');
alpha0 = rand(p,1); 
alpha0 = alpha0 ./ sum(alpha0); 

rng(myseed+1, 'twister');
mu0 = randn(p,1);

n0         = length(theta0);
alpha0_rss = repmat(alpha0, [1 n0]);
mu0_rss    = repmat(mu0, [1 n0]);

% fit baseline model 
tic;
[b_logw,b_alpha,b_mu,b_s] = null_wrapper_fixsb('squarem',betahat,se,SiRiS,sigb,theta0,alpha0_rss,mu0_rss);
fprintf('Baseline model analysis is finished ...\n');

% fit enrichment model
[log10_bf,e_logw,e_alpha,e_mu,e_s] = gsea_wrapper_fixsb('squarem',betahat,se,SiRiS,snps,sigb,theta0,theta,b_logw,b_alpha,b_mu);
fprintf('Enrichment model analysis is finished ...\n');
rss_time = toc;

% save the output
file_name = 'example5_simulated_results.mat'; 
save(file_name,'log10_bf','e_logw','e_alpha','e_mu','e_s','b_logw','b_alpha','b_mu','b_s','rss_time');
