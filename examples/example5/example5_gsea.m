% Note that the working directory here is assumed to be `rss/examples/example5`.
% Please modify the following code accordingly if a different directory is used.

clear;

% add rss-varbvsr to the search path
addpath('../../src_vb');

% fitting method and random seed
stage  = 'step1';
myseed = 459;
method = 'squarem';

% specify the gene set name
% IL23-mediated signaling events (Pathway Commons 2, PID, 37 genes)
geneset_name = 'path2641'; 

% specify the path of summary data across whole genome
allstat_path = 'ibd2015_sumstat.mat';

% trait name and sample size
trait_name  = 'ibd2015';
sample_size = (12882+21770); % Table 1; cases: 12,882; controls: 21,770

% parameter grids
h      = 0.3;
theta0 = (-2.9:0.025:-2.85)';
theta  = (0:(3/100):3)';

% specify the path of the results of null analysis
data_path = './';

% specify the path of subset of summary data corresponding to the gene set 
sumstat_path = 'ibd2015_sumstat_path2641.mat';

% specify the output path
output_path = './'; 

%===== MODIFICATION STOPS =====

% run gsea on the gene set
run('gsea_template.m') 
