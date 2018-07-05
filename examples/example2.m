clear;

% Set this to the directory containing example1.mat and where the output
% files will be stored.
working_dir = '.';

% add search paths
addpath(genpath('../src'));

% load summary-level data
fprintf('Loading data.\n');
example_data = matfile(cat(2,working_dir,'/','example2.mat'));

betahat = example_data.betahat; % single-SNP effect size estimates
se 	= example_data.se; 	% standard error of betahat
Nsnp    = example_data.Nsnp; 	% sample size of each SNP

% mcmc length for rss-bvsr
Ndraw = 1e6;
Nburn = 1e5;
Nthin = 9e1;

% load the three types of ld matrix
genotype_data = matfile(cat(2,working_dir,'/','genotype2.mat'));

% fit rss-bvsr model when R is the sample correlation in the panel
cohort_R = genotype_data.cohort_R;
tic;
[betasam, gammasam, hsam, logpisam, Naccept] = ...
  rss_bvsr(betahat, se, cohort_R, Nsnp, Ndraw, Nburn, Nthin);
runtime = toc;
save(cat(2,working_dir,'/','example2_rss_c.mat'),...
     'betasam', 'gammasam', 'hsam', 'logpisam', 'Naccept',...
     'runtime', '-v7.3');
clearvars betasam gammasam hsam logpisam Naccept runtime;
fprintf('RSS-C is done ... \n');

% fit rss-bvsr model when R is the sample correlation in the panel
shrink_R = genotype_data.shrink_R;
tic;
[betasam, gammasam, hsam, logpisam, Naccept] = ...
  rss_bvsr(betahat, se, shrink_R, Nsnp, Ndraw, Nburn, Nthin);
runtime = toc;
save(cat(2,working_dir,'/','example2_rss_s.mat'),...
     'betasam', 'gammasam', 'hsam', 'logpisam', 'Naccept',...
     'runtime', '-v7.3');
clearvars betasam gammasam hsam logpisam Naccept runtime;
fprintf('RSS is done ... \n');

% fit rss-bvsr model when R is the sample correlation in the panel
panel_R = genotype_data.panel_R;
tic;
[betasam, gammasam, hsam, logpisam, Naccept] = ...
  rss_bvsr(betahat, se, panel_R, Nsnp, Ndraw, Nburn, Nthin);
runtime = toc;
save(cat(2,working_dir,'/','example2_rss_p.mat'),...
     'betasam', 'gammasam', 'hsam', 'logpisam', ...
     'Naccept', 'runtime', '-v7.3');
clearvars betasam gammasam hsam logpisam Naccept runtime;
fprintf('RSS-P is done ... \n');

exit;

