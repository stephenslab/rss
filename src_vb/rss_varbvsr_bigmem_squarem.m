function [lnZ, alpha, mu, s, info] = rss_varbvsr_bigmem_squarem(file, sigb, logodds, options)
% USAGE: mean-field variational inference of RSS-BVSR model for a given set of hyperparameters
%        use parfor function (MATLAB Parallel Computing Toolbox) for parallel calculations
%        use SQUAREM as an accelerator (with step length modification to ensure monotonicity)
%        this function is more efficient than rss_varbvsr_pasquarem.m in terms of memory usage
% INPUT:
%       file: the mat file that contains cell arrays of betahat, se and SiRiS, string
%       sigb: prior SDs of the regression coefficients (if included), p by 1 or scalar
%       logodds: log(prior PIP/(1-prior PIP)) of inclusion for each SNP, p by 1 or scalar
%       options: user-specified behaviour of the algorithm, structure
%		- max_walltime: scalar, the maximum wall time (unit: seconds) for this program
%		- tolerance: scalar, convergence tolerance
%		- alpha & mu: p by 1 vectors, initial values of variational parameters
%		- verbose: logical, print program progress if true 
%               - modify_step: logical, modify the step length in SQUAREM if true
%               - elbo_tol: scalar, stop iteration when ELBO increase is below this value
% OUTPUT:
%       lnZ: scalar, variational lower bound of the marginal log likelihood (up to some constant)
%       alpha: p by 1, variational estimates of the posterior inclusion probabilities 
%       mu: p by 1, posterior means of the additive effects (if the SNP is included)
%       s: p by 1, posterior variances of the additive effects (if the SNP is included)
%       info: structure with following fields 
%               - iter: integer, number of iterations till convergence
%               - maxerr: maximum relative difference between the parameters at the last 2 iterations
%               - loglik: iter by 1, variational lower bound at each iteration
% NOTE:
%       The pseudocode for SQUAREM is available at Table 1 of Varadhan and Roland (2008):
%       https://doi.org/10.1111/j.1467-9469.2007.00585.x

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
 
  % Optional: stop the variational updates when
  % the increase in lower bound (ELBO) is small.
  if isfield(options,'elbo_tol')
    elbo_tol = double(options.elbo_tol);
  else
    elbo_tol = NaN;
  end

  % Load summary statistics from the mat file.
  sumstat = matfile(file);
  betahat = sumstat.betahat;
  se      = sumstat.se;

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
  
  % Partition the whole genome data into chromosomes.
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

    % This part is to deal with the big memory requirement issue.
    SiRiS_tmpcell      = sumstat.SiRiS(c,1); %#ok<PFBNS>
    SiRiS_tmp          = SiRiS_tmpcell{1};
    SiRiSr_cell{c,1}   = full(SiRiS_tmp * (alpha_cell{c,1} .* mu_cell{c,1}));

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

    alpha_r_cell = cell(C, 1);
    mu_r_cell    = cell(C, 1);
    alpha_v_cell = cell(C, 1);
    mu_v_cell    = cell(C, 1); 

    alpha_tmp_cell  = cell(C, 1);
    mu_tmp_cell     = cell(C, 1);
    SiRiSr_tmp_cell = cell(C, 1);

    alpha_r_norm2 = zeros(C, 1);
    mu_r_norm2    = zeros(C, 1);
    alpha_v_norm2 = zeros(C, 1);
    mu_v_norm2    = zeros(C, 1);
 
    parfor c = 1:C

      % Get the number of SNPs on Chr. c.
      pchr = length(betahat{c,1});

      % This part is to deal with the big memory requirement issue.
      SiRiS_tmpcell = sumstat.SiRiS(c,1); %#ok<PFBNS>
      SiRiS_tmp     = SiRiS_tmpcell{1};

      % All SNPs are included in the forward/backward updates.
      if mod(iter,2)
        I = (1:pchr);
      else
        I = (pchr:-1:1);
      end

      % Run the first fixed-point mapping step (line 1 of Table 1).
      [alpha_tmp1,mu_tmp1,SiRiSr_tmp1] = rss_varbvsr_update(SiRiS_tmp,sigb_cell{c,1},...
                                         logodds_cell{c,1},betahat{c,1},se{c,1},...
                                         alpha0_cell{c,1},mu0_cell{c,1},SiRiSr0_cell{c,1},I);

      % Run the second fixed-point mapping step (line 2 of Table 1).
      [alpha_tmp2,mu_tmp2,SiRiSr_tmp2] = rss_varbvsr_update(SiRiS_tmp,sigb_cell{c,1},...
                                         logodds_cell{c,1},betahat{c,1},se{c,1},...
                                         alpha_tmp1, mu_tmp1, SiRiSr_tmp1, I);

      % Compute the step length (line 3-4 of Table 1).
      alpha_r_cell{c,1} = alpha_tmp1 - alpha0_cell{c,1};
      mu_r_cell{c,1}    = mu_tmp1 - mu0_cell{c,1};
      alpha_v_cell{c,1} = (alpha_tmp2 - alpha_tmp1) - alpha_r_cell{c,1};
      mu_v_cell{c,1}    = (mu_tmp2 - mu_tmp1) - mu_r_cell{c,1};

      alpha_r_norm2(c) = sum(alpha_r_cell{c,1}.^2);
      mu_r_norm2(c)    = sum(mu_r_cell{c,1}.^2);
      alpha_v_norm2(c) = sum(alpha_v_cell{c,1}.^2);
      mu_v_norm2(c)    = sum(mu_v_cell{c,1}.^2);

      % Save the output of the second fix-point mapping output (for mtp >= -1 case).
      alpha_tmp_cell{c,1}  = alpha_tmp2;
      mu_tmp_cell{c,1}     = mu_tmp2;
      SiRiSr_tmp_cell{c,1} = SiRiSr_tmp2;

    end
    
    % Compute the step length (line 5 of Table 1).
    % The eps term is added to avoid the 0/0 case.
    mtp = - sqrt(sum(alpha_r_norm2)+sum(mu_r_norm2)) / sqrt(sum(alpha_v_norm2)+sum(mu_v_norm2)+eps);

    % Modifiy the step length under three different scenarios.
    alpha_tmp1_cell  = cell(C, 1);
    mu_tmp1_cell     = cell(C, 1);
    SiRiSr_tmp1_cell = cell(C, 1);

    % Scenario 1: use the output of the second fix-point mapping output.
    if mtp >= -1
      for c = 1:C
        alpha_tmp1_cell{c,1}   = alpha_tmp_cell{c,1};
        mu_tmp1_cell{c,1}      = mu_tmp_cell{c,1};
        SiRiSr_tmp1_cell{c,1}  = SiRiSr_tmp_cell{c,1};
      end
    % Scenario 2: use the modified step length mtp computed above (line 7 of Table 1).
    else
      parfor c = 1:C
        % This part is to deal with the big memory requirement issue.
        SiRiS_tmpcell         = sumstat.SiRiS(c,1); %#ok<PFBNS>
        SiRiS_tmp             = SiRiS_tmpcell{1};
        alpha_tmp1_cell{c,1}  = alpha0_cell{c,1} - 2*mtp*alpha_r_cell{c,1} + (mtp^2)*alpha_v_cell{c,1};
        mu_tmp1_cell{c,1}     = mu0_cell{c,1} - 2*mtp*mu_r_cell{c,1} + (mtp^2)*mu_v_cell{c,1};
        SiRiSr_tmp1_cell{c,1} = full(SiRiS_tmp * (alpha_tmp1_cell{c,1} .* mu_tmp1_cell{c,1}));
      end 
    end

    % Run the last fix-point mapping step for Scenario 1 or 2 (line 8 of Table 1).  
    parfor c = 1:C
      pchr = length(betahat{c,1});

      % This part is to deal with the big memory requirement issue.
      SiRiS_tmpcell = sumstat.SiRiS(c,1); %#ok<PFBNS>
      SiRiS_tmp     = SiRiS_tmpcell{1};

      if mod(iter,2)
        I = (1:pchr);
      else
        I = (pchr:-1:1);
      end

      alpha_tmp3  = alpha_tmp1_cell{c,1};
      mu_tmp3     = mu_tmp1_cell{c,1};
      SiRiSr_tmp3 = SiRiSr_tmp1_cell{c,1};

      [alpha_tmp,mu_tmp,SiRiSr_tmp] = rss_varbvsr_update(SiRiS_tmp,sigb_cell{c,1},...
                                      logodds_cell{c,1},betahat{c,1},se{c,1},...
                                      alpha_tmp3,mu_tmp3,SiRiSr_tmp3,I);

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

      absr_uni(c) = max(abs(r));
      asum_uni(c) = sum(alpha_tmp);
    end

    % Compute the lower bound to the marginal log-likelihood for whole genome.
    lnZ = sum(lnZ_cell);

    % Scenario 3: use a simple back-tracking to modify the steplength iteratively.
    if (mtp < -1) && (lnZ < lnZ0)
      num_bt = 0;
      while (lnZ < lnZ0) && (num_bt < 10)
        mtp = 0.5*(mtp-1); % back-tracking

        parfor c = 1:C
          pchr = length(betahat{c,1});

          % This part is to deal with the big memory requirement issue.
          SiRiS_tmpcell = sumstat.SiRiS(c,1); %#ok<PFBNS>
          SiRiS_tmp     = SiRiS_tmpcell{1};

          if mod(iter,2)
            I = (1:pchr);
          else
            I = (pchr:-1:1);
          end

          alpha_tmp3  = alpha0_cell{c,1} - 2*mtp*alpha_r_cell{c,1} + (mtp^2)*alpha_v_cell{c,1};
          mu_tmp3     = mu0_cell{c,1} - 2*mtp*mu_r_cell{c,1} + (mtp^2)*mu_v_cell{c,1};
          SiRiSr_tmp3 = full(SiRiS_tmp * (alpha_tmp3 .* mu_tmp3));

          [alpha_tmp,mu_tmp,SiRiSr_tmp] = rss_varbvsr_update(SiRiS_tmp,sigb_cell{c,1},...
                                          logodds_cell{c,1},betahat{c,1},se{c,1},...
                                          alpha_tmp3, mu_tmp3, SiRiSr_tmp3, I);

          alpha_cell{c,1}  = alpha_tmp;
          mu_cell{c,1}     = mu_tmp;
          r                = alpha_tmp .* mu_tmp;
          SiRiSr_cell{c,1} = SiRiSr_tmp;

          lnZ_cell(c) = (q_cell{c,1})'*r - 0.5*r'*SiRiSr_cell{c,1} + intgamma(logodds_cell{c,1},alpha_cell{c,1});
          lnZ_cell(c) = lnZ_cell(c) - 0.5*(1./sesquare_cell{c,1})'*betavar(alpha_cell{c,1},mu_cell{c,1},s_cell{c,1});
          lnZ_cell(c) = lnZ_cell(c) + intklbeta_rssbvsr(alpha_cell{c,1},mu_cell{c,1},s_cell{c,1},sigbsquare_cell{c,1});
        end

        lnZ    = sum(lnZ_cell);
        num_bt = num_bt + 1; % stop back-tracking after ten steps
      end

      % Aggregate per-chromosome results to create stopping rule.
      parfor c = 1:C
        r                = alpha_cell{c,1} .* mu_cell{c,1};
        params_tmp       = [alpha_cell{c,1}; r];
        params_cell{c,1} = params_tmp;
        J_tmp            = find(abs(params_tmp) > 1e-6);
        params0_tmp      = params0_cell{c,1}
        err_tmp          = relerr(params_tmp(J_tmp),params0_tmp(J_tmp));
        maxerr_uni(c)    = max(err_tmp);
 
        absr_uni(c) = max(abs(r));
        asum_uni(c) = sum(alpha_cell{c,1});
      end

    end

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
      status = sprintf('%4d %+13.6e %0.1e %4d %0.2f\n',iter,lnZ,maxerr,asum,absr);
      fprintf(status);
      fprintf(repmat('\b',1,length(status)));
    end

    % Terminate the for loop if ELBO decreases.
    if lnZ < lnZ0

      fprintf('\nWarning: log variational lower bound decreased by %+0.2e\n',lnZ0-lnZ);
      alpha = cell2mat(alpha0_cell);
      mu    = cell2mat(mu0_cell);
      lnZ   = lnZ0;
      break

    % Terminate the for loop if variational parameters converge.
    elseif maxerr < tolerance

      alpha = cell2mat(alpha_cell);
      mu    = cell2mat(mu_cell);
      fprintf('\nConverged: maximum relative error %+0.2e\n',maxerr);
      fprintf('Log variational lower bound of the last step increased by %+0.2e\n',lnZ-lnZ0);
      break

    end

    % Optional: terminate the for loop if ELBO increases by a small value.
    elbo_change = lnZ - lnZ0;
    if ~isnan(elbo_tol) && elbo_change <= elbo_tol

      alpha = cell2mat(alpha_cell);
      mu    = cell2mat(mu_cell);
      fprintf('\nConverged: log variational lower bound increased by %+0.2e\n',lnZ-lnZ0);
      break

    end

    % Terminate the for loop after the given maximum wall time.
    exe_time = etime(clock, start_time);
    if exe_time >= max_walltime

      alpha = cell2mat(alpha_cell);
      mu    = cell2mat(mu_cell);
      fprintf('\nMaximum wall time reached: %+0.2e seconds\n',exe_time);
      fprintf('Log variational lower bound of the last step increased by %+0.2e\n',lnZ-lnZ0);
      break

    end
    
  end

  % Show total execution time in seconds.
  if exe_time < max_walltime
    fprintf('\nTotal execution time: %+0.2e seconds\n',exe_time);
  end

  % Save info as a structure array.
  info = struct('iter',iter,'maxerr',maxerr,'loglik',loglik,'exe_time',exe_time);

end

