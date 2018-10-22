function [lnZ, alpha, mu, s, info] = rss_varbvsr_parallel(betahat, se, SiRiS, sigb, logodds, options)
% USAGE: mean-field variational inference of RSS-BVSR model for a given set of hyperparameters
%        use parfor function (MATLAB Parallel Computing Toolbox) for parallel calculations
% INPUT:
%       betahat: effect size estimates under single-SNP model, C by 1 cell array
%       se: standard errors of betahat, C by 1 cell array
%       SiRiS: inv(S)*R*inv(S), double precision sparse matrix (CCS format), C by 1 cell array
%       sigb: prior SDs of regression coefficients (if included), p by 1 or scalar
%       logodds: log(prior PIP/(1-prior PIP)) of inclusion for each SNP, p by 1 or scalar
%       options: user-specified behaviour of the algorithm, structure
%               - max_walltime: scalar, the maximum wall time (unit: seconds) for this program
%               - tolerance: scalar, convergence tolerance
%               - alpha & mu: p by 1 vectors, initial values of variational parameters
%               - verbose: logical, print program progress if true
% OUTPUT:
%       lnZ: scalar, variational lower bound of the marginal log likelihood (up to some constant)
%       alpha: p by 1, variational estimates of the posterior inclusion probabilities 
%       mu: p by 1, posterior means of the additive effects (if the SNP is included)
%       s: p by 1, posterior variances of the additive effects (if the SNP is included)
%       info: structure with following fields 
%               - iter: integer, number of iterations till convergence
%               - maxerr: maximum relative difference between the parameters at the last 2 iterations
%               - loglik: iter by 1, variational lower bound at each iteration

  % Get the time when the program starts.
  start_time = clock;

  if ~exist('options', 'var')
    options = [];
  end

  % Set the maximum wall time for this program (unit: seconds).
  if isfield(options,'max_walltime')
    max_walltime = double(options.max_walltime);
  else
    max_walltime = (1*24+11)*3600; % 1 day 11 hours
  end

  % Set tolerance for convergence, which is reached when the maximum relative distance
  % between successive updates of the variational parameters is less than this quantity.
  if isfield(options,'tolerance')
    tolerance = double(options.tolerance);
  else
    tolerance = 1e-4;
  end
  fprintf('Tolerance for convergence in this program: %0.2e \n', tolerance);

  % Get the number of analyzed SNPs in the whole genome (p).
  p = length(cell2mat(betahat));

  % Set the hyper-parameters (sigb and logodds).
  if isscalar(sigb)
    disp('Prior SDs (sigb) of all SNPs are the same.');
    sigb = repmat(sigb,p,1);
  end
  if isscalar(logodds)
    disp('Prior log-odds (logodds) of all SNPs are the same.');
    logodds = repmat(logodds,p,1);
  end

  % SiRiS must be a sparse matrix.
  if ~issparse(SiRiS)
    SiRiS = sparse(double(SiRiS));
  end

  % Set initial estimates of variational parameters.
  if isfield(options,'alpha')
    alpha = double(options.alpha(:));
  else
    alpha = rand(p,1);
    alpha = alpha / sum(alpha);
  end
  if isfield(options,'mu')
    mu = double(options.mu(:));
  else
    mu = randn(p,1);
  end
  if length(alpha) ~= p || length(mu) ~= p
    error('options.alpha and options.mu must be vectors of the same length');
  end

  % Determine whether to display the algorithm's progress.
  if isfield(options,'verbose')
    verbose = options.verbose;
  else
    verbose = true;
  end

  clear options;

  % Partition the whole genome data by chromosomes.
  C      = length(betahat);
  chrpar = partition_genome(betahat); 

  lnZ_cell        = zeros(C, 1);
  alpha_cell      = cell(C, 1);
  mu_cell         = cell(C, 1);
  params_cell     = cell(C, 1); 
  s_cell          = cell(C, 1);
  SiRiSr_cell     = cell(C, 1);
  q_cell          = cell(C, 1);
  sesquare_cell   = cell(C, 1);
  logodds_cell    = cell(C, 1);
  sigb_cell       = cell(C, 1);
  sigbsquare_cell = cell(C, 1);

  for c = 1:C
    chr_start             = chrpar(c, 1); 
    chr_end               = chrpar(c, 2);
    alpha_cell{c, 1}      = alpha(chr_start:chr_end); 
    mu_cell{c, 1}         = mu(chr_start:chr_end);
    logodds_cell{c, 1}    = logodds(chr_start:chr_end);
    sigb_cell{c, 1}       = sigb(chr_start:chr_end);
    sigbsquare_cell{c, 1} = sigb(chr_start:chr_end).^2; 
  end

  % Compute a few useful quantities for the main loop.
  parfor c = 1:C
    params_cell{c,1}   = [alpha_cell{c,1}; alpha_cell{c,1} .* mu_cell{c,1}];
    SiRiSr_cell{c,1}   = full(SiRiS{c,1} * (alpha_cell{c,1} .* mu_cell{c,1}));
    sesquare_cell{c,1} = se{c,1} .* se{c,1};
    q_cell{c,1}        = betahat{c,1} ./ sesquare_cell{c,1};
    s_cell{c,1}        = (sesquare_cell{c,1} .* sigbsquare_cell{c,1}) ./ (sesquare_cell{c,1} + sigbsquare_cell{c,1});
  end

  % Aggregate the variational estimates of variance.
  s = cell2mat(s_cell);
  
  % Initialize the fields of the structure info.
  iter   = 0;
  loglik = [];

  % Calculate the variational lower bound based on the initial values.
  parfor c = 1:C
    r = alpha_cell{c,1} .* mu_cell{c,1};

    lnZ_cell(c) = (q_cell{c,1})'*r - 0.5*r'*SiRiSr_cell{c,1} + intgamma(logodds_cell{c,1},alpha_cell{c,1});
    lnZ_cell(c) = lnZ_cell(c) - 0.5*(1./sesquare_cell{c,1})'*betavar(alpha_cell{c,1},mu_cell{c,1},s_cell{c,1});
    lnZ_cell(c) = lnZ_cell(c) + intklbeta_rssbvsr(alpha_cell{c,1},mu_cell{c,1},s_cell{c,1},sigbsquare_cell{c,1});
  end
  lnZ = sum(lnZ_cell);
  fprintf('Calculate the variational lower bound based on the initial values: %+13.6e ...\n', lnZ);

  loglik = [loglik; lnZ]; 

  if verbose
    fprintf('       variational    max. incl max.\n');
    fprintf('iter   lower bound  change vars E[b]\n');
  end

  % Repeat until convergence criterion is met.
  while true

    % Go to the next iteration.
    iter = iter + 1;
    
    % Save the current variational parameters and lower bound.
    lnZ0         = lnZ;
    alpha0_cell  = alpha_cell; 
    mu0_cell     = mu_cell;
    params0_cell = params_cell;
    SiRiSr0_cell = SiRiSr_cell;

    % Parallel for loop over each chromosome.
    maxerr_uni = zeros(C, 1);
    absr_uni   = zeros(C, 1);
    asum_uni   = zeros(C, 1);
 
    parfor c = 1:C

      % Get the number of SNPs on Chr. c
      pchr = length(betahat{c,1});

      % All SNPs are included in the forward/backward updates.
      if mod(iter,2)
        I = (1:pchr);
      else
        I = (pchr:-1:1);
      end

      % Run a forward or backward pass of the coordinate ascent updates on Chr. c.
      [alpha_tmp,mu_tmp,SiRiSr_tmp] = rss_varbvsr_update(SiRiS{c,1},sigb_cell{c,1},...
                                      logodds_cell{c,1},betahat{c,1},se{c,1},...
                                      alpha0_cell{c,1},mu0_cell{c,1},SiRiSr0_cell{c,1},I);

      alpha_cell{c,1}  = alpha_tmp;
      mu_cell{c,1}     = mu_tmp;
      r                = alpha_tmp .* mu_tmp;
      params_tmp       = [alpha_tmp; r];
      params_cell{c,1} = params_tmp;
      SiRiSr_cell{c,1} = SiRiSr_tmp;

      % Compute the lower bound to the marginal log-likelihood of Chr. c.
      lnZ_cell(c) = (q_cell{c,1})'*r - 0.5*r'*SiRiSr_cell{c,1} + intgamma(logodds_cell{c,1},alpha_cell{c,1});
      lnZ_cell(c) = lnZ_cell(c) - 0.5*(1./sesquare_cell{c,1})'*betavar(alpha_cell{c,1},mu_cell{c,1},s_cell{c,1});
      lnZ_cell(c) = lnZ_cell(c) + intklbeta_rssbvsr(alpha_cell{c,1},mu_cell{c,1},s_cell{c,1},sigbsquare_cell{c,1});

      % Compute some quantities related to the stopping rule.
      J_tmp         = find(abs(params_tmp) > 1e-6);
      params0_tmp   = params0_cell{c,1}
      err_tmp       = relerr(params_tmp(J_tmp),params0_tmp(J_tmp));
      maxerr_uni(c) = max(err_tmp);
      absr_uni(c)   = max(abs(r));
      asum_uni(c)   = sum(alpha_tmp);

    end

    % Aggregate per-chromosome results to create stopping rule.
    lnZ    = sum(lnZ_cell);
    maxerr = max(maxerr_uni);
    absr   = max(absr_uni);
    asum   = round(sum(asum_uni));

    % Record the variational lower bound at each iteration.
    loglik = [loglik; lnZ]; %#ok<AGROW>

    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small.
    
    if verbose
      status = sprintf('%4d %+13.6e %0.1e %4d %0.2f',iter,lnZ,maxerr,asum,absr);
      fprintf(status);
      fprintf(repmat('\b',1,length(status)));
    end

    if lnZ < lnZ0

      if verbose
        fprintf('\n');
        fprintf('WARNING: the log variational lower bound decreased by %+0.2e\n',lnZ0-lnZ);
      end
      alpha = cell2mat(alpha0_cell);
      mu    = cell2mat(mu0_cell);
      lnZ   = lnZ0;
      break

    elseif maxerr < tolerance

      alpha = cell2mat(alpha_cell);
      mu    = cell2mat(mu_cell);
      if verbose
        fprintf('\n');
        fprintf('Convergence reached: maximum relative error %+0.2e\n',maxerr);
        fprintf('The log variational lower bound of the last step increased by %+0.2e\n',lnZ-lnZ0);
      end
      break

    end

    % Terminate the for loop after the given maximum wall time.
    exetime = etime(clock, start_time);
    if exetime >= max_walltime

      alpha = cell2mat(alpha_cell);
      mu    = cell2mat(mu_cell);
      if verbose
        fprintf('\n');
        fprintf('Maximum wall time reached: %+0.2e seconds\n',exetime);
        fprintf('The log variational lower bound of the last step increased by %+0.2e\n',lnZ-lnZ0);
      end
      break

    end

  end

  % Save info as a structure array.
  info = struct('iter',iter,'maxerr',maxerr,'loglik',loglik);
end

