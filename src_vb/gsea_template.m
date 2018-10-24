% USAGE: a template script to run gene set enrichment analysis on GWAS summary data
% INPUT:
%	allstat_path: the path of mat file that contains GWAS summary data of whole genome, string
%	sumstat_path: the path of mat file that contains subset of GWAS summary data based on gene set, string
%	data_path: the path of mat file that contains null analysis results, string
%	output_path: the path of mat file that stores analysis output, string
%	trait_name: the name of GWAS, character
%	geneset_name: the name of gene set, character
%	sample_size: the sample size of GWAS data, integer
%	theta0: the grid of the genome-wide log-odds (base 10), n0 by 1 
%	theta: the grid of the enrichment (base 10), n1 by 1
%	h: the fixed proportion of phenotypic variance explained by available genotypes, nh by 1
%	stage: the analysis stage, e.g. 'step 1', 'step 2', ...
%	myseed: seed value of generating random numbers, scalar 
%	method: the implementation of rss-varbvsr, character

% set the fixed part
sumstat  = matfile(allstat_path);
se       = sumstat.se; 
p        = length(cell2mat(se));
Nsnp_all = sample_size * ones(p, 1);
se_all   = cell2mat(se);
clear sumstat se p;

% prepare the input data for gsea
sumstat = matfile(sumstat_path);
snps    = sumstat.snps;

switch method

  case {'original', 'squarem'} % serial implementation: arrays as input
    betahat = sumstat.betahat;
    se      = sumstat.se;
    SiRiS   = sumstat.SiRiS;

  case {'parallel', 'pasquarem'} % parallel implementation: cell arrays as input
    sumstat_cell = partition_geneset(sumstat);
    disp('Input data have been transformed from arrays to cell arrays.');

    betahat = sumstat_cell.betahat;
    se      = sumstat_cell.se;
    SiRiS   = sumstat_cell.SiRiS;

    clear sumstat_cell;

    % start the matlabpool with maximum available workers
    % control how many workers by setting ntasks in your sbatch script
    pc = parcluster('local');
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
end 

% make sure that SNPs assigned to multiple genes are only included once
% i.e. there should be NO duplicated elements in the vector snps
if length(snps) ~= length(unique(snps))
  error('There exist duplicated SNPs!');
end

if all(snps == unique(snps))
  disp('PASS CHECK: SNPs are uniquely included.');
end

% initialize the gsea analysis
p  = length(snps);
n0 = length(theta0);
n1 = length(theta);
nh = length(h);
fprintf('Name of gene set: %s ...\n', geneset_name);
fprintf('Number of SNPs in the gene set: %d \n', p);
fprintf('Grid size of theta0: %d \n', n0);
fprintf('Grid size of theta: %d \n', n1);
fprintf('Grid size of h: %d \n', nh);

% load the variational estimates from the null analysis
logw0  = zeros(n0, nh);
alpha0 = zeros(p, n0, nh);
mu0    = zeros(p, n0, nh);

for k=1:nh
  my_h = h(k);

  for i=1:n0
    my_theta0 = theta0(i);

    my_name = strcat(trait_name,'_null_h_',num2str(round(my_h*100)),'_theta0_',num2str(round(abs(my_theta0)*100)));
    my_name = strcat(my_name,'_seed_',num2str(myseed),'_squarem_',stage,'.mat');

    my_data = matfile(strcat(data_path,my_name));

    % some sanity checks
    if my_h ~= my_data.h
      error('The value of h is wrong!');
    end
    if my_theta0 ~= my_data.theta0
      error('The value of theta0 is wrong!');
    end
    if myseed ~= my_data.myseed
      error('The value of myseed is wrong!');
    end
    %if ~strcmp(method, my_data.method)
    %  error('The value of method is wrong!');
    %end
    if ~strcmp(stage, my_data.stage)
      error('The value of stage is wrong!');
    end

    logw0(i,k) = my_data.logw;

    alpha0_total = my_data.alpha;
    mu0_total    = my_data.mu;

    alpha0(:,i,k) = alpha0_total(snps,1); clear alpha0_total;
    mu0(:,i,k)    = mu0_total(snps,1); clear mu0_total;

  end
end

% randomly initialize the variational estimates for the enrichment analysis
rng(myseed, 'twister');
alpha = rand(p,n0,n1,nh);
alpha = alpha ./ repmat(sum(alpha),p,1);
rng(myseed+1, 'twister');
mu    = randn(p,n0,n1,nh);

% run the enrichment analysis using rss-varbvsr
tic;
[log10bf,logw1,alpha,mu,s] = gsea_wrapper(method,betahat,se,SiRiS,se_all,Nsnp_all,snps,h,theta0,theta,logw0,alpha0,mu0,alpha,mu);
runtime = toc;

% remove unused quantities from memory
clear betahat se SiRiS se_all Nsnp_all alpha0 mu0

% save the output
output_name = strcat(output_path,trait_name,'_gsea_seed_',num2str(myseed),'_',geneset_name,'_',method,'.mat');
save(output_name,'method','myseed','snps','h','theta0','theta','log10bf','logw1','logw0','alpha','mu','s','runtime','-v7.3');

