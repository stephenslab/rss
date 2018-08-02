% input data file: ibd2015_sumstat.mat
% source: Liu, van Sommeren et al, 2015 (PUBMED: 26192919)

% please name the summary data file as
% "trait_name"_sumstat.mat
% in order to run the following code

% specify trait-specific information
trait_name  = 'ibd2015';
sample_size = (12882+21770); % cases: 12,882; controls: 21,770

% specify grid of hyper-parameters
h_rv      = 0.3;
theta0_rv = (-2.9:0.025:-2.85);

% specify stage of analysis
stage = 'step1';

% specify random start and algorithm 
myseed = 459;
method = 'squarem';

% specify input and output paths
input_path  = './'; 
output_path = './';

%%%%%%%%%%%%% DO NOT MODIFY CODES BELOW %%%%%%%%%%%%%

% add search path
addpath('../../src_vb/');

% create all combinations of two vectors
% https://www.mathworks.com/help/nnet/ref/combvec.html
hyper_para = combvec(h_rv, theta0_rv);

% specify hyper-parameters based on job id 
h      = hyper_para(1, case_id);
theta0 = hyper_para(2, case_id);

% run analyses
run('null_template.m');

exit;
