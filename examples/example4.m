clear;

% add search paths
addpath(genpath('../src_vb'));
addpath('../misc');
addpath('/home/xiangzhu/varbvs-master/varbvs-MATLAB/');

% set the number of replications
prompt = 'What is the number of replications? ';
Nrep   = input(prompt);

% set the genome-wide log-odds 
prompt1 = 'What is the genome-wide log-odds? ';
theta0  = input(prompt1);

% set the log-fold enrichment
prompt2 = 'What is the log-fold enrichment? ';
theta   = input(prompt2);

thetatype = [theta0 theta];

% set the true pve
prompt3 = 'What is the pve (between 0 and 1)? ';
pve     = input(prompt3);

% preallocate output
alpha_n_diff = zeros(Nrep, 1); % max absolute difference of estimated alpha from null model
alpha_e_diff = zeros(Nrep, 1); % max absolute difference of estimated alpha from enrichment model
mu_n_diff    = zeros(Nrep, 1); % max absolute difference of estimated mu from null model
mu_e_diff    = zeros(Nrep, 1); % max absolute difference of estimated mu from enrichment model
s_n_diff     = zeros(Nrep, 1); % max absolute difference of estimated s from null model
s_e_diff     = zeros(Nrep, 1); % max absolute difference of estimated s from enrichment model
bf           = zeros(Nrep, 2); % absolute difference of estimated Bayes factors 
bf_reldiff   = zeros(Nrep, 1); % relative difference of estimated Bayes factors

% generate the individual-level and summary-level data
genotype = matfile('genotype.mat');
C 	 = genotype.C;
[n,p] 	 = size(C);

AH    = matfile('AH_chr16.mat');
H     = AH.H;      			% hypothesis matrix 3323x3160
A     = AH.A;      			% annotation matrix 12758x3323
paths = find(H(:,end));			% signal transduction (Biosystem, Reactome)
snps  = find(sum(A(:,paths),2));    	% index of variants inside the pathway
clear AH A H paths;

% start the replications of checking
for i = 1:Nrep
  myseed = 617 + i;

  % generate data under enrichment hypothesis
  [true_para,individual_data,summary_data] = enrich_datamaker(C,thetatype,pve,myseed,snps);
  fprintf('Individual-level and summary-level data are ready ...\n');

  % fix all hyper-parameters as their true values
  sigma   = true_para{4}^2;       % the true residual variance
  sa      = 1/sigma;              % sa*sigma being the true prior variance of causal effects
  sigb    = 1;                    % the true prior SD of causal effects 
  theta0  = thetatype(1);         % the true genome-wide log-odds (base 10)
  theta   = thetatype(2);         % the true log-fold enrichment (base 10)

  % initialize the variational parameters
  rng(myseed, 'twister');
  alpha0 = rand(p,1); 
  alpha0 = alpha0 ./ sum(alpha0); 

  rng(myseed+1, 'twister');
  mu0 = randn(p,1);

  % load the individual-level data
  y = individual_data{1};
  X = individual_data{2};
  clear individual_data;

  % set the summary-level data for a perfect match b/t rss-varbvsr and varbvs
  betahat = summary_data{1}; 
  se 	  = sqrt(sigma ./ diag(X'*X)); 		% condition 1 for perfect matching
  Si 	  = 1 ./ se(:);
  R 	  = corrcoef(X); 			% condition 2 for perfect matching
  SiRiS   = sparse(repmat(Si, 1, p) .* R .* repmat(Si', p, 1)); 
  clear R Si summary_data; 

  % run varbvs on the individual-level data
  options_n = struct('maxiter',1e8,'sa',sa,'logodds',theta0,'alpha',alpha0,...
		     'mu',mu0,'sigma',sigma,'initialize_params',false,'verbose',false);

  fit_null = varbvs(X,[],y,[],'gaussian',options_n);

  ns = length(theta);
  logodds = repmat(theta0,p,ns);
  logodds(snps,:) = repmat(theta0+theta',length(snps),1);
  options_e = struct('maxiter',1e8,'sa',sa,'logodds',logodds,'alpha',alpha0,...
                     'mu',mu0,'sigma',sigma,'initialize_params',false,'verbose',false);

  fit_gsea = varbvs(X,[],y,[],'gaussian',options_e);

  bf_b = exp(fit_gsea.logw - fit_null.logw);
  clear options_n options_e X y;

  fprintf('Analysis of individual-level data via VARBVS is done ...\n');

  % run rss-varbvsr on the summary-level data
  logodds_n 	  = log(10)*theta0;
  logodds_e 	  = repmat(logodds_n,p,1);
  logodds_e(snps) = log(10)*(theta0+theta);

  options = struct('alpha',alpha0,'mu',mu0,'verbose',false);

  [logw_nr,alpha_nr,mu_nr,s_nr] = rss_varbvsr(betahat,se,SiRiS,sigb,logodds_n,options);
  [logw_er,alpha_er,mu_er,s_er] = rss_varbvsr(betahat,se,SiRiS,sigb,logodds_e,options);

  bf_r = exp(logw_er-logw_nr);
  clear options betahat se SiRiS;

  fprintf('Analysis of summary-level data via RSS-VARBVSR is done ...\n');

  % record the differences
  alpha_nb = fit_null.alpha;
  alpha_eb = fit_gsea.alpha;
  mu_nb    = fit_null.mu;
  mu_eb    = fit_gsea.mu;
  s_nb     = fit_null.s;
  s_eb     = fit_gsea.s;

  alpha_n_diff(i) = max(abs(alpha_nb(:)-alpha_nr(:)));
  alpha_e_diff(i) = max(abs(alpha_eb(:)-alpha_er(:)));
  mu_n_diff(i)    = max(abs(mu_nb(:)-mu_nr(:)));
  mu_e_diff(i)    = max(abs(mu_eb(:)-mu_er(:)));
  s_n_diff(i)     = max(abs(s_nb(:)-s_nr(:)));
  s_e_diff(i)     = max(abs(s_eb(:)-s_er(:)));
  bf(i,:)         = [bf_b bf_r];
  bf_reldiff(i)   = (bf_b-bf_r)/bf_b;

  fprintf('Trial %d is done ...\n', i);

end

% Save the output
output_name = strcat('example4-theta0',num2str(-theta0),'-theta',num2str(theta),'-',num2str(Nrep), 'trials.mat');
save(output_name,'alpha_n_diff','alpha_e_diff','mu_n_diff','mu_e_diff','s_n_diff','s_e_diff','bf','bf_reldiff');
