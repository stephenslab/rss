function [lnZ, alpha, mu, s, info] = rss_varbvsr_squarem(betahat, se, SiRiS, sigb, logodds, options)
% USAGE: mean-field variational approximation of the RSS-BVSR model given the hyperparameters
%	 with SQUAREM added as an accelerator (with steplength modification to ensure monotonicity)
% INPUT:
%       betahat: the effect size estimates under single-SNP model, p by 1
%       se: standard errors of betahat, p by 1
%       SiRiS: inv(S)*R*inv(S), double precision sparse matrix (ccs format), p by p
%       sigb: the prior SD of the regression coefficients (if included), scalar
%       logodds: the prior log-odds (i.e. log(prior PIP/(1-prior PIP))) of inclusion for each SNP, p by 1
%       options: user-specified behaviour of the algorithm, structure
%               - max_walltime: scalar, the maximum wall time (unit: seconds) for this program
%               - tolerance: scalar, convergence tolerance
%               - alpha & mu: p by 1 vectors, initial values of variational parameters
%               - verbose: logical, print program progress if true
% OUTPUT:
%	lnZ: scalar, the variational lower bound of the marginal log likelihood (up to some constant)
%	alpha: p by 1, variational estimates of the posterior inclusion probabilities 
%	mu: p by 1, posterior means of the additive effects (given snp included)
%	s: p by 1, posterior variances of the additive effects (given snp included)
%	info: structure with following fields 
%		- iter: integer, number of iterations
%       	- maxerr: the maximum relative difference between the parameters at the last two iterations
%		- sigb: scalar, the maximum likelihood estimate of sigma_beta
%		- loglik: iter by 1, the variational lower bound at each iteration
%               - exetime: scalar, the execution time of the program

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

  % Get the number of variables (p).
  p = length(betahat);

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

  % Determine whether to modify the step length in SQUAREM (step 6 in Table 1).
  if isfield(options,'modify_step')
    modify_step = options.modify_step;
  else
    modify_step = true;
  end

  clear options;
  
  % Compute a few useful quantities for the main loop.
  SiRiSr = full(SiRiS * (alpha .* mu));
  q 	 = betahat ./ (se .^2);

  % Calculate the variance of the coefficients.
  se_square 	= se .* se;
  sigb_square 	= sigb * sigb;
  s 		= (se_square .* sigb_square) ./ (se_square + sigb_square);

  % Initialize the fields of the structure info.
  iter   = 0;
  loglik = [];

  % Calculate the variational lower bound based on the initial values.
  r   = alpha .* mu;
  lnZ = q'*r - 0.5*r'*SiRiSr - 0.5*(1./se_square)'*betavar(alpha, mu, s);
  lnZ = lnZ + intgamma(logodds, alpha) + intklbeta_rssbvsr(alpha, mu, s, sigb_square);
  fprintf('Calculate the variational lower bound based on the initial values: %+13.6e ...\n', lnZ);
  
  loglik = [loglik; lnZ]; 

  if verbose
    fprintf('       variational    max. incl max.       \n');
    fprintf('iter   lower bound  change vars E[b] sigma2\n');
  end

  % Repeat until convergence criterion is met.
  while true

    % Go to the next iteration.
    iter = iter + 1;
    
    % Save the current variational parameters and lower bound.
    alpha0  = alpha; 
    mu0     = mu;
    lnZ0    = lnZ;
    params0 = [alpha; alpha .* mu];

    % Run a forward or backward pass of the coordinate ascent updates.
    if mod(iter,2)
      I = (1:p);
    else
      I = (p:-1:1);
    end

    % Run the first fix-point mapping step (line 1 of Table 1).
    [alpha1, mu1, SiRiSr] = rss_varbvsr_update(SiRiS, sigb, logodds, betahat, se, alpha0, mu0, SiRiSr, I); 
    % Run the second fix-point mapping step (line 2 of Table 1).
    [alpha2, mu2, SiRiSr] = rss_varbvsr_update(SiRiS, sigb, logodds, betahat, se, alpha1, mu1, SiRiSr, I);

    % Compute the step length (line 3-5 of Table 1).
    alpha_r 	= alpha1 - alpha0;
    mu_r 	= mu1 - mu0;
    alpha_v 	= (alpha2 - alpha1) - alpha_r;
    mu_v 	= (mu2 - mu1) - mu_r;
    mtp 	= - sqrt(norm(alpha_r)^2+norm(mu_r)^2) / sqrt(norm(alpha_v)^2+norm(mu_v)^2);

    % Modifiy the step length (optional): three scenarios (Section 6).
    if modify_step
      % scenario 1: use the output of the second fix-point mapping output
      if mtp >= -1 					 % set mtp = -1
        alpha3  = alpha2;
        mu3     = mu2;
        SiRiSr3 = SiRiSr;
      % scenario 2: no need to modify the step length (line 7 of Table 1)  
      else 						% i.e. mtp < -1
        alpha3  = alpha0 - 2*mtp*alpha_r + (mtp^2)*alpha_v;
        mu3     = mu0 - 2*mtp*mu_r + (mtp^2)*mu_v;
        SiRiSr3 = full(SiRiS * (alpha3 .* mu3));
      end
    else % skip the modification of step length
      alpha3  = alpha0 - 2*mtp*alpha_r + (mtp^2)*alpha_v;
      mu3     = mu0 - 2*mtp*mu_r + (mtp^2)*mu_v;
      SiRiSr3 = full(SiRiS * (alpha3 .* mu3));
    end

    % Run the last fix-point mapping step for scenario 1 and 2 (line 8 of Table 1).
    [alpha, mu, SiRiSr] = rss_varbvsr_update(SiRiS, sigb, logodds, betahat, se, alpha3, mu3, SiRiSr3, I);
    r 		  	= alpha .* mu;
    lnZ 		= q'*r - 0.5*r'*SiRiSr - 0.5*(1./se_square)'*betavar(alpha, mu, s);
    lnZ 		= lnZ + intgamma(logodds, alpha) + intklbeta_rssbvsr(alpha, mu, s, sigb_square);
   
    % scenario 3: use a simple back-tracking to modify the steplength iteratively
    if modify_step && (mtp < -1) && (lnZ < lnZ0)
      num_bt = 0;
      while (lnZ < lnZ0) && (num_bt < 10)
        mtp 		    = 0.5*(mtp-1); % back-tracking
    	alpha3              = alpha0 - 2*mtp*alpha_r + (mtp^2)*alpha_v;
    	mu3                 = mu0 - 2*mtp*mu_r + (mtp^2)*mu_v;
    	SiRiSr3             = full(SiRiS * (alpha3 .* mu3));
    	[alpha, mu, SiRiSr] = rss_varbvsr_update(SiRiS, sigb, logodds, betahat, se, alpha3, mu3, SiRiSr3, I);
    	r                   = alpha .* mu;
    	lnZ                 = q'*r - 0.5*r'*SiRiSr - 0.5*(1./se_square)'*betavar(alpha, mu, s);
    	lnZ                 = lnZ + intgamma(logodds, alpha) + intklbeta_rssbvsr(alpha, mu, s, sigb_square);
        num_bt              = num_bt + 1; % stop back-tracking after ten steps
      end
    end

    % Record the variational lower bound at each iteration.
    loglik = [loglik; lnZ]; %#ok<AGROW>

    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small.
    params = [alpha; r];
    I      = find(abs(params) > 1e-6);
    err    = relerr(params(I),params0(I));
    maxerr = max(err);

    if verbose
      status = sprintf('%4d %+13.6e %0.1e %4d %0.2f %5.2f',...
		       iter,lnZ,maxerr,round(sum(alpha)),max(abs(r)),sigb_square);
      fprintf(status);
      fprintf(repmat('\b',1,length(status)));
    end

    if lnZ < lnZ0

      if verbose
        fprintf('\n');
        fprintf('WARNING: the log variational lower bound decreased by %+0.2e\n',lnZ0-lnZ);
      end
      alpha  = alpha0;
      mu     = mu0;
      lnZ    = lnZ0;
      sigb   = sqrt(sigb_square);
      break

    elseif maxerr < tolerance

      sigb = sqrt(sigb_square);
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

      sigb = sqrt(sigb_square);
      if verbose
        fprintf('\n');
        fprintf('Maximum wall time reached: %+0.2e seconds\n',exetime);
        fprintf('The log variational lower bound of the last step increased by %+0.2e\n',lnZ-lnZ0);
      end
      break

    end

  end

  % Save info as a structure array.
  info = struct('iter',iter,'maxerr',maxerr,'sigb',sigb,'loglik',loglik,'exetime',exetime);

end
