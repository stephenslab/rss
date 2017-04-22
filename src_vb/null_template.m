% NB: this script is written only for the environment of RCC midway clusters
% modification maybe required if this script is used in other environments

% copy the input data to /srv/scratch/local on midway partitions
% this is to improve code performance when reading whole genome data
input_file = strcat(input_path,'/',trait_name,'_sumstat.mat');
copyfile(input_file,'/srv/scratch/local/');

% use the copy on local scratch as input data
sumstat_file = strcat('/srv/scratch/local/',trait_name,'_sumstat.mat');

% stop the execution temporarily for 1-10 minutes based on the job id
% this is to avoid license issues when running many jobs in the same time
n_time = abs(mod(case_id, 10))*60;
pause(n_time);

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
pc = parcluster('local');
matlabpool(pc, getenv('SLURM_CPUS_ON_NODE')); %#ok<DPOOL>

% get the number of analyzed SNPs and the (max) sample size
sumstat = matfile(sumstat_file);
betahat = sumstat.betahat; 
p 	= length(cell2mat(betahat));
Nsnp 	= sample_size * ones(p, 1);

% check whether the input data are collected from unique SNPs
chr = double(cell2mat(sumstat.chr));
chr = strread(num2str(chr'),'%s');
pos = double(cell2mat(sumstat.pos));
pos = strread(num2str(pos'),'%s');
snp = strcat({'chr'},chr,{':'},pos);

if all(strcmp(sort(snp),unique(snp)))
  disp('SANITY CHECK PASS: there are no duplicated SNPs ...');
else
  error('ERROR: there exist duplicated SNPs ...');
end
clear sumstat betahat chr pos snp;

% set up initial values for model fitting
switch(stage)

  % step 1: choose an initialization common to all hyper-parameter combinations;
  % this is done by setting the same random seed for all combinations
  case 'step1'
    rng(myseed, 'twister');
    alpha0 = rand(p, 1); 
    alpha0 = alpha0 ./ sum(alpha0);
 
    rng(myseed+1, 'twister');
    mu0 = randn(p, 1);

    fprintf('Use the same random start as initial values for all hyper-parameter combinations ... \n');

  % step 2: choose an initialization common to all the runs of the coordinate
  % ascent algorithm; this is chosen from the hyperparameters with the highest
  % variational estimate of the posterior probability in step 1
  % NB: the file saved in step1_output_path should be generated before step 2
  case 'step2'
    step1_output_path = strrep(output_path,'step2','step1');
    step1_output_path = strcat(step1_output_path,trait_name,'_null_seed_',num2str(myseed),'_',method,'_step1.mat');

    step1_output = matfile(step1_output_path);

    logw  = step1_output.logw;
    alpha = step1_output.alpha;
    mu    = step1_output.mu;

    [~, i] = max(logw(:));
    alpha0 = alpha(:,i);
    mu0    = mu(:,i);

    clear logw alpha mu;
    fprintf('Use the optimal output from step 1 as initial values for all hyper-parameter combinations ... \n');

end

% print job-specific information
fprintf('The stage of analysis: %s ...\n', stage);
fprintf('The path of output: %s ...\n', output_path);
fprintf('The value of theta0: %0.2f ...\n', theta0);
fprintf('The value of h: %0.2f ...\n', h);

% run the whole genome null analysis using rss-varbvsr
tic;
[logw, alpha, mu, s] = null_single(method, sumstat_file, Nsnp, h, theta0, alpha0, mu0);
runtime = toc;

% save the output
output_name = strcat(output_path,trait_name);
output_name = strcat(output_name,'_null_h_',num2str(round(h*100)));
output_name = strcat(output_name,'_theta0_',num2str(round(abs(theta0)*100)));
output_name = strcat(output_name,'_seed_',num2str(myseed),'_',method,'_',stage,'.mat');

save(output_name,'stage','method','myseed','h','theta0','logw','alpha','mu','s','runtime','-v7.3');

% delete the copy of input data on local scratch of the node
% this deletion is acceptable b/c each job takes a single node
delete(sumstat_file);

