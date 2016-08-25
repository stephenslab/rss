clear;

% Set this to the directory containing example1.mat and where the output
% files will be stored.
working_dir = '/tmp/pcarbo';

% add search paths
addpath(genpath('../src'));

% load summary-level data
fprintf('Loading data.\n');
example_data = matfile(cat(2,working_dir,'/','example1.mat'));

R   = example_data.R; 	% population LD matrix
bwd = example_data.bwd; % bandwidth of R 
BR  = example_data.BR; 	% banded storage of R

betahat = example_data.betahat; % single-SNP effect size estimates
se 	= example_data.se; 	% standard error of betahat
Nsnp 	= example_data.Nsnp; 	% sample size of each SNP

fprintf('Data set is loaded ... \n');

% compare the two definitions of se
se_1 = se; 					% the simple version
se_2 = sqrt((betahat.^2) ./ Nsnp + se.^2); 	% the rigorous version 

abs_diff = abs(se_1 - se_2);
fprintf('Look at the five-number summary of the absolute difference between two definitions of SE: \n')
disp(percentile(log10(abs_diff),0:0.25:1));

% fit rss-bvsr model with the first definition of se
fprintf('Start RSS-BVSR analysis on the first definition of SE... \n');

% mcmc info
Ndraw = 2e6;
Nburn = 2e5;
Nthin = 9e1;

tic;

% 1. Simulate posterior samples via mcmc.
[betasam, gammasam, hsam, logpisam, Naccept] = ...
    rss_bvsr(betahat, se_1, R, Nsnp, Ndraw, Nburn, Nthin);

% 2. Compute the posterior samples of PVE.
matrix_type 	= 1;
M 		= length(hsam);
pvesam 		= zeros(M,1);
progress_bar 	= progress('init','start PVE calculation');
for i = 1:M
  pvesam(i) 	= compute_pve(betasam(i,:), betahat, se_1, Nsnp, bwd, BR, matrix_type);
  progress_bar 	= progress(progress_bar, i/M);
end
runtime = toc;

% 3. Save output.
save(cat(2,working_dir,'/','example3_se1.mat'),...
     'betasam', 'gammasam', 'hsam', 'logpisam', 'pvesam',...
     'Naccept', 'runtime', '-v7.3');
clearvars betasam gammasam hsam logpisam pvesam Naccept runtime;

% fit rss-bvsr model with the second definition of se
fprintf('Start RSS-BVSR analysis on the second definition of SE... \n');

tic;

% 1. Simulate posterior samples via mcmc.
[betasam, gammasam, hsam, logpisam, Naccept] = ...
    rss_bvsr(betahat, se_2, R, Nsnp, Ndraw, Nburn, Nthin);

% 2. Compute the posterior samples of PVE.
matrix_type     = 1;
M               = length(hsam);
pvesam          = zeros(M,1);
progress_bar    = progress('init','start PVE calculation');
for i = 1:M
        pvesam(i)       = compute_pve(betasam(i,:), betahat, se_2, Nsnp, bwd, BR, matrix_type);
        progress_bar    = progress(progress_bar, i/M);
end
runtime = toc;

% 3. Save output.
save(cat(2,working_dir,'/','example3_se2.mat'),...
     'betasam', 'gammasam', 'hsam', 'logpisam', 'pvesam',...
     'Naccept', 'runtime', '-v7.3');
clearvars betasam gammasam hsam logpisam pvesam Naccept runtime;
