function [R, BR] = get_corr(m, Ne, cummap, Hpanel, cutoff, isgeno, negcm)
% USAGE: compute LD matrix using the shrinkage estimator in Wen and Stephens (2010)
% INPUT:
%	m: the number of individuals in the reference panel, integer
%	Ne: the effective population size (diploid), integer
%	cummap: cumulative genetic map in cM, numSNP by 1
%	Hpanel: the (phased) haplotypes from a reference panel, numIND by numSNP
%	cutoff: the hard threshold for small entries being zero, scalar 
%	isgeno: true if Hpanel is an unphased genotype matrix, logical
%	negcm: true if there is negative genetic distance, logical
% OUTPUT:
%	R: the estimated LD matrix, numSNP by numSNP, sparse matrix
%	BR: the banded storage of R, dense matrix

  % decide whether Hpanel is haplotype or genotype
  if ~exist('isgeno', 'var')
    isgeno = false;
  end

  if isgeno
    disp('Input genetic data are genotypes.');
  else
    disp('Input genetic data are haplotypes.');
  end

  % decide whether there exists negative genetic distance or not
  if ~exist('negcm', 'var')
    negcm = false;
  end

  if negcm
    disp('Negative genetic distance is allowed here.');
  else
    disp('Negative genetic distance is not allowed here.');
  end

  % compute the shrinkage estimator of covariance matrix
  disp('Compute Wen-Stephens shrinkage LD estimator ...'); 
  SigHat = shrink_cov(m, Ne, cummap, Hpanel, cutoff, isgeno, negcm);

  % convert covariance to correlation
  R = corrcov(SigHat);
  clear SigHat;

  if nargout > 1
    disp('Convert LD matrix to a banded storage ...');
 
    % get the bandwidth of R
    bwd = find_bandwidth(R);
  
    % get the banded storage of R
    BR = band_storage(R, bwd);
  end

  % store R as a sparse matrix
  R = sparse(R);
end

function SigHat = shrink_cov(m, Ne, cummap, Hpanel, cutoff, isgeno, negcm)
% USAGE: compute the shrinkage estimator of covariance matrix in Wen and Stephens (2010)
% INPUT:
%	m: number of individuals in a reference panel, integer
%	Ne: effective population size (diploid), integer
%	cummap: cumulative genetic map in cM, numSNP by 1
%	Hpanel: phased haplotypes or unphased genotypes in a reference panel, numIND by numSNP
%	cutoff: hard threshold for small entries being zero, scalar
%	isgeno: true if Hpanel is an unphased genotype matrix, logical
%	negcm: true if there is negative genetic distance, logical 
% OUTPUT:
%	SigHat: estimated covariance matrix of haplotype, numSNP by numSNP

  % theta is related to mutation
  nmsum = sum(1 ./ (1:(2*m-1)));
  theta = (1/nmsum) / (2*m + 1/nmsum);
	
  % S is obtained from Sigma_panel by shrinking off-diagonal entries toward 0
  
  % Scenario 1: Hpanel is a phased haplotype matrix
  S = cov(Hpanel);

  % Scenario 2: Hpanel is an unphased genotype matrix
  if isgeno
    disp('Hpanel is an unphased genotype matrix.');
    S = 0.5*S;
  end

  S = triu(S);

  % compute the values on and above the 1st diagonal
  numSNP = size(S,1);

  for i = 1:numSNP
    for j = (i+1):numSNP

      % compute genetic distance between SNPs i and j
      genetic_dist = cummap(j) - cummap(i);

      % deal with potential negative distance value
      % which may occur after lifting genome builds
      if (genetic_dist<0) & (~negcm)
        error('Negative recombination rate is not allowed.')
      end

      if (genetic_dist<0) & (negcm)
        genetic_dist = abs(genetic_dist);
      end

      % compute scaled population recombination rate
      % see https://stephenslab.github.io/rss/Recombination
      rho = 4 * Ne * genetic_dist / 100;

      % compute the shrinkage factor based on recombination rate
      shrinkage = exp(-rho/(2*m));

      % apply hard thresholding to obtain sparse and banded estimator
      if shrinkage < cutoff
        shrinkage = 0;
      end

      % shrink the sample covariance in the reference panel
      S(i,j) = shrinkage * S(i,j);

    end
  end

  % copy the upper half (NO diagonal) to the lower half
  S = S + triu(S,1)';
	
  % SigHat is derived from Li and Stephens model (2003)
  SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);

end

