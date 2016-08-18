clear;

% add search paths
addpath(genpath('../src'));

% load summary-level data
fprintf('Loading data.\n');
example_data = matfile('/tmp/pcarbo/example1.mat');

R   = example_data.R; 	% population LD matrix
bwd = example_data.bwd; % bandwidth of R 
BR  = example_data.BR; 	% banded storage of R

betahat = example_data.betahat; % single-SNP effect size estimates
se 	= example_data.se; 	% standard error of betahat
Nsnp 	= example_data.Nsnp; 	% sample size of each SNP

fprintf('Data set is loaded ... \n');

% Check model assumption.
chat = sqrt((betahat(:).^2) ./ (Nsnp(:).*(se(:).^2) + betahat(:).^2));
fprintf('Look at the five-number summary of log10 sample\n');
fprintf('phenotype-genotype correlations: \n')
disp(percentile(log10(chat),0:0.25:1));

% Fit rss-bvsr model.
fprintf('Start RSS-BVSR analysis ... \n');
Ndraw = 2e6;
Nburn = 2e5;
Nthin = 9e1;
tic;

% 1. Simulate posterior samples via mcmc.
[betasam, gammasam, hsam, logpisam, Naccept] = ...
    rss_bvsr(betahat, se, R, Nsnp, Ndraw, Nburn, Nthin);

% 2. Compute the posterior samples of PVE.
matrix_type 	= 1;
M 		= length(hsam);
pvesam 		= zeros(M,1);
progress_bar 	= progress('init','start PVE calculation');
for i = 1:M
  pvesam(i) = compute_pve(betasam(i,:), betahat, se, Nsnp, bwd, BR,...
                          matrix_type);
  progress_bar	= progress(progress_bar, i/M);
end
runtime = toc;

% 3. Save output.
save('/tmp/pcarbo/example1_rssbvsr.mat', 'betasam', 'gammasam', 'hsam',...
     'logpisam','pvesam', 'Naccept', 'runtime', '-v7.3');
clearvars betasam gammasam hsam logpisam pvesam Naccept runtime;
fprintf('RSS-BVSR analysis is done ... \n');

% fit rss-bslmm model
fprintf('Start RSS-BSLMM analysis ... \n');
Ndraw = 2e6;
Nburn = 2e5;
Nthin = 9e1;
tic;   

% 1. simulate posterior samples via mcmc
[bsam, zsam, lpsam, hsam, rsam, Naccept] = ...
    rss_bslmm(betahat, se, R, Nsnp, Ndraw, Nburn, Nthin);

% 2. compute the posterior samples of pve
matrix_type 	= 1;
M 		= length(hsam);
pvesam 		= zeros(M,1);
progress_bar 	= progress('init','start PVE calculation');
for i = 1:M 
  pvesam(i) = compute_pve(bsam(i,:), betahat, se, Nsnp, bwd, BR, matrix_type);
  progress_bar = progress(progress_bar, i/M);
end
runtime = toc;

% 3. save output
save('/tmp/pcarbo/example1_rssbslmm.mat', 'bsam', 'zsam', 'lpsam', 'hsam',...
     'rsam', 'pvesam', 'Naccept', 'runtime', '-v7.3');
clearvars bsam zsam lpsam hsam rsam pvesam Naccept runtime;
fprintf('RSS-BSLMM analysis is done ... \n');

% fit rss-ash model
fprintf('Start RSS-ASH analysis ... \n');
Ndraw = 5e7;
Nburn = 1e7;
Nthin = 1e3;
sigma_beta = [0 0.001 0.003 0.01 0.03 0.1 0.3 1 3];
tic;

% 1. simulate posterior samples via mcmc
[bsam, zsam, wsam, lsam, Naccept] =...
    rss_ash(betahat, se, R, Nsnp, sigma_beta, Ndraw, Nburn, Nthin);

% 2. compute the posterior samples of pve
matrix_type 	= 1;
M 		= length(lsam);
pvesam 		= zeros(M,1);
progress_bar 	= progress('init','start PVE calculation');
for i = 1:M
  pvesam(i) = compute_pve(bsam(i,:), betahat, se, Nsnp, bwd, BR, matrix_type);
  progress_bar	= progress(progress_bar, i/M);
end
runtime = toc;

% 3. save output
save('/tmp/pcarbo/example1_rssash.mat', 'bsam', 'zsam', 'wsam', 'lsam', ...
     'pvesam', 'Naccept', 'runtime', '-v7.3');
clearvars bsam zsam wsam lsam pvesam Naccept runtime;
fprintf('RSS-ASH analysis is done ... \n');

exit;
