function [pve_mean,pve_sample] = estimate_pve(betahat,se,R,nsam,logw,alpha,mu,s,nsim)
% USAGE: estimate PVE from summary data 
% INPUT:
%       betahat: effect size estimates under single-SNP model, C by 1 cell array
%       se: standard errors of betahat, C by 1 cell array
%       R: estimated LD matrices, C by 1 cell array
%       nsam: sample size of each SNP, integer or p by 1 array
%       logw: k by 1, variational lower bound of the marginal log likelihood (up to some constant)
%       alpha: p by k, variational estimates of the posterior inclusion probabilities 
%       mu: p by k, posterior means of the additive effects (if the SNP is included)
%       s: p by k, posterior variances of the additive effects (if the SNP is included)
%       nsim: size of simulated samples, integer
% OUTPUT:
%	pve_mean: scalar, posterior mean for PVE
%	pve_sample: vector, simulated posterior sample for PVE

% NB: as a sanity check, pve_mean should be close to mean(pve_sample).

  % set the default size of simulated sample 
  if ~exist('nsim','var')
    nsim = 1e4;
  end

  % ensure input data are all stored as C by 1 cell arrays
  cell_check = prod([iscell(betahat),iscell(se),iscell(R)]);
  if cell_check ~= 1
    error('Input data must be stored as cell arrays ...');
  end

  % partition whole genome to chromosomes
  chr_par = partition_genome(betahat);

  % convert sample size vector to a cell array
  if isscalar(nsam)
    nsam_vec = nsam*ones(length(cell2mat(betahat)),1);
  else
    nsam_vec = nsam;
  end
  clear nsam;
  nsam = vec2cell(nsam_vec, chr_par);
  clear nsam_vec;
  
  % obtain the posterior probability distribution for hyper-parameter
  w = normalizelogweights(logw);

  % compute (once) a quantity used later
  betahat_vec = cell2mat(betahat);
  se_vec      = cell2mat(se);
  nsam_vec    = cell2mat(nsam);
  ns2_vec     = nsam_vec .* (se_vec.^2) + betahat_vec.^2;
  clear betahat_vec se_vec nsam_vec;  

  % compute posterior mean PVE for each hyper-parameter
  num_para = length(w);
  pve_para = zeros(num_para, 1);

  for k=1:num_para

    % compute the first part: trace(Cov(b)*A)
    beta_var  = alpha(:,k).*(mu(:,k).^2+s(:,k)) - (alpha(:,k).*mu(:,k)).^2;
    pve_part1 = sum(beta_var ./ ns2_vec); 

    % compute the second part: E(b)'*A*E(b)
    beta_mean = vec2cell(alpha(:,k).*mu(:,k), chr_par);
    pve_part2 = compute_pve(beta_mean, betahat, se, R, nsam);

    % use quadratic expectation formula, that is,
    % E(b'*A*b) = trace(Cov(b)*A) + E(b)'*A*E(b)
    pve_para(k) = pve_part1 + pve_part2;
    clear beta_var beta_mean pve_part1 pve_part2;

  end

  % average posterior mean PVE over all hyper-parameter settings
  pve_mean = sum(w .* pve_para);
  clear pve_para ns2_vec num_para; 

  % simulate posterior sample of PVE
  if nargout > 1

    % get the total number of SNPs
    num_snp = length(cell2mat(betahat));

    % initialize storage for simulated posterior sample
    pve_sample = zeros(nsim, 1);

    % compute PVE for each simulated value
    for i=1:nsim

      % first draw a hyper-parameter setting
      j = randtable(w);

      % next sample regression coefficients ('beta')
      beta_vec = mu(:,j) + sqrt(s(:,j)) .* randn(num_snp,1);
      beta_vec = beta_vec .* (rand(num_snp,1) < alpha(:,j));

      % convert sampled 'beta' to a cell array
      beta = vec2cell(beta_vec, chr_par);
      clear beta_vec;

      % compute PVE for the sampled 'beta' value
      pve_sample(i) = compute_pve(beta, betahat, se, R, nsam);
      clear beta;
  
    end

  end

end

function out_cell = vec2cell(in_vec, par)
% USAGE: convert a vector to a cell array
% INPUT:
%	in_vec: input vector
%	par: [par_start par_end], specifying the cell partition
% OUTPUT:
%	out_cell: output cell array

  out_cell = cell(length(par), 1);

  for c=1:length(par)
    par_start     = par(c,1);
    par_end       = par(c,2);
    out_cell{c,1} = in_vec(par_start:par_end);
  end

  % confirm that out_cell == in_vec
  if ~all(cell2mat(out_cell) == in_vec)
    error('Error occurred in vector to cell conversion ...');
  end

end

function pve = compute_pve(beta, betahat, se, R, nsam)
% USAGE: compute PVE for a given 'beta' value from summary-level data
% INPUT:
%	beta: regression coefficients drawn from the variational posterior, C by 1 cell array
%	betahat: effect size estimates under single-SNP model, C by 1 cell array
%	se: standard errors of betahat, C by 1 cell array
%	R: estimated LD matrices, C by 1 cell array
%	nsam: sample size of each SNP, C by 1 cell array
% OUTPUT:
%	pve: PVE value for the given 'beta' value, scalar

  % preallocate output
  num_cell = length(beta);
  pve_cell = cell(num_cell,1);

  % compute PVE for each cell
  for c=1:num_cell

    % use Equation 3.7 in Zhu and Stephens (2017)
    yxsqr = sqrt(nsam{c,1} .* (se{c,1}.^2) + betahat{c,1}.^2);
    comp1 = beta{c,1} ./ yxsqr;
    comp2 = R{c,1} * comp1;

    pve_cell{c} = sum(comp1 .* comp2);
    clear yxsqr comp1 comp2

  end

  % aggregate results
  pve = sum(cell2mat(pve_cell));

end

function w = normalizelogweights(logw)
% Source: https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/normalizelogweights.m
% Author: Peter Carbonetto 

  % Guard against underflow or overflow by adjusting the
  % log-probabilities so that the largest probability is 1.
  c = max(logw(:));
  w = exp(logw - c);

  % Normalize the probabilities.
  w = w(:) / sum(w(:));

end

function x = randtable(p, n)
% Source: https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/randtable.m
% Author: Peter Carbonetto

  % Get the number of random draws to generate.
  if ~exist('n','var')
    n = 1;
  end

  % Get the number of classes.
  m = numel(p);

  % Normalize the probability table.
  p = p(:)';
  p = p/sum(p);
  
  % Create the CDF.
  ub = cumsum(p);
  lb = [0 ub(1:m-1)];
  
  % Generate the discrete random deviates.
  x = zeros(n,1);
  for i = 1:n
    u    = repmat(rand,1,m);
    x(i) = sum((lb <= u & u < ub) .* (1:m));
  end

end
