% Note that the working directory here is assumed to be `rss/examples/example5`.
% Please modify the following code accordingly if a different directory is used.

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

% initialize variational parameters
myseed = 200;

rng(myseed, 'twister');
alpha0 = rand(p,1); 
alpha0 = alpha0 ./ sum(alpha0); 

rng(myseed+1, 'twister');
mu0 = randn(p,1);

n0         = length(theta0);
alpha0_rss = repmat(alpha0, [1 n0]);
mu0_rss    = repmat(mu0, [1 n0]);

% fit the baseline model 
tic;
[b_logw,b_alpha,b_mu,b_s] = null_wrapper_fixsb('squarem',betahat,se,SiRiS,sigb,theta0,alpha0_rss,mu0_rss);
fprintf('Baseline model analysis is finished ...\n');

% fit the enrichment model
[log10_bf,e_logw,e_alpha,e_mu,e_s] = gsea_wrapper_fixsb('squarem',betahat,se,SiRiS,snps,sigb,theta0,theta,b_logw,b_alpha,b_mu);
fprintf('Enrichment model analysis is finished ...\n');

% specify pre-defined genomic segments (genes)
segs_file = 'Aseg_chr16.mat';

% generate gene-level results under baseline model
[b_p1,b_p2] = compute_pip(segs_file,b_logw,b_alpha);

% generate gene-level results under enrichment model

% note that e_logw is 21x21 and e_alpha is 12758x21x21
% but compute_pip.m only allows logw vector and alpha matrix
num_para    = length(theta0) * length(theta);
e_alpha_mat = zeros(p, num_para);
e_logw_vec  = zeros(num_para, 1);

col_count = 0;
for i = 1:length(theta0)
  for j = 1:length(theta)
    col_count = col_count + 1;

    e_logw_vec(col_count)    = e_logw(i,j);
    e_alpha_mat(:,col_count) = e_alpha(:,i,j);  
  end
end
 
[e_p1,e_p2] = compute_pip(segs_file,e_logw_vec,e_alpha_mat);
clear e_logw_vec e_alpha_mat;

fprintf('Gene-level results are generated ...\n');
rss_time = toc;

% save results
file_name = 'example5_simulated_results.mat'; 
save(file_name,'log10_bf','e_logw','e_alpha','e_mu','e_s','e_p1','e_p2','b_logw','b_alpha','b_mu','b_s','b_p1','b_p2','rss_time');

fprintf('Full analysis results are saved ...\n');

% look at some global results

fprintf('Log 10 enrichment Bayes factor: %.4f ...\n', log10_bf);

fprintf('Total number of SNPs with non-zero effect: %d ...\n', sum(example_data.gamma));

b_w   = exp(b_logw - max(b_logw(:)));
b_w   = b_w / sum(b_w(:));
b_ens = sum(b_alpha);
b_ens = dot(b_ens(:), b_w(:));
fprintf('Estimated number of SNPs with non-zero effect under baseline: %.4f ...\n', b_ens);

e_w   = exp(e_logw - max(e_logw(:)));
e_w   = e_w / sum(e_w(:));
e_ens = sum(e_alpha);
e_ens = dot(e_ens(:), e_w(:));
fprintf('Estimated number of SNPs with non-zero effect under enrichment: %.4f ...\n', e_ens);

% look at gene-level results

gene_info  = matfile(segs_file);
gene_start = gene_info.gene_start;
gene_stop  = gene_info.gene_stop;

% load SNP-level causal status
snp_pos = gene_info.pos;
snp_cau = example_data.gamma;

d = 100e3;

% annotate whether a gene is causal
gene_cau = zeros(length(gene_start),1);

for i=1:length(gene_start)
  snp_index   = find(snp_pos > gene_start(i)-d & snp_pos < gene_stop(i)+d);
  gene_cau(i) = max(snp_cau(snp_index));
end

disp([b_p1(1:10) e_p1(1:10) gene_cau(1:10)]);
disp([b_p1(gene_cau==1) e_p1(gene_cau==1)]);
