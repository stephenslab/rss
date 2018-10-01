clear;

% add rss-varbvsr to the search path
addpath('../../src_vb');

% specify model fitting method and random seed
stage  = 'step1';
myseed = 459;
method = 'squarem';

% specify the gene set name
% IL23-mediated signaling events (Pathway Commons 2, PID, 37 genes)
geneset_name = 'path2641'; 

% specify trait name
trait_name  = 'ibd2015';

% specify parameter grids
h      = 0.3;
theta0 = (-2.9:0.025:-2.85)';
theta  = (0:(3/100):3)';

% specify the path of baseline and enrichment results
data_path = './';

% specify the output path
output_path = './'; 

%===== MODIFICATION STOPS =====

% specify pre-defined genomic segments (genes)
segs_file = strcat(trait_name,'_',geneset_name,'_genes.mat');

% get the total number of genome-wide SNPs
load(segs_file);
p = length(pos);

% extract baseline model fitting results
b_logw  = zeros(length(theta0), 1);
b_alpha = zeros(p, length(theta0)); 

for i = 1:length(theta0)
  my_name = strcat(trait_name,'_null_h_',num2str(round(h*100)),'_theta0_',num2str(round(abs(theta0(i))*100)));
  my_name = strcat(my_name,'_seed_',num2str(myseed),'_squarem_',stage,'.mat');
  my_data = matfile(strcat(data_path,my_name));

  b_logw(i,1)  = my_data.logw;
  b_alpha(:,i) = my_data.alpha;
end

% generate gene-level results under baseline model
[b_p1,b_p2] = compute_pip(segs_file,b_logw,b_alpha);

% extract enrichment model fitting results
num_para = length(theta0) * length(theta);
e_alpha  = zeros(p, num_para);
e_logw   = zeros(num_para, 1);

gsea_path = strcat(data_path,trait_name,'_gsea_seed_',num2str(myseed),'_',geneset_name,'_',method,'.mat'); 

enrichment_results = matfile(gsea_path);
alpha_enrichment   = enrichment_results.alpha;
logw_enrichment    = enrichment_results.logw1;

snps = enrichment_results.snps;

col_count = 0;
for i = 1:length(theta0)
  alpha_baseline = b_alpha(:,i);

  for j = 1:length(theta)
    col_count = col_count + 1;

    e_logw(col_count,1) = logw_enrichment(i,j);
 
    % retain the baseline posterior estimates outside the gene set
    e_alpha(:,col_count) = alpha_baseline;
    % use the enrichment posterior estimates inside the gene set
    e_alpha(snps,col_count) = alpha_enrichment(:,i,j); 
  end

end

% generate gene-level results under enrichment model
[e_p1,e_p2] = compute_pip(segs_file,e_logw,e_alpha);

% save gene priortization results
output_file = strcat(output_path,'ibd2015_path2641_genes_results.mat');
save(output_file,'b_p1','b_p2','e_p1','e_p2','gene_chr','gene_id','gene_start','gene_stop');
