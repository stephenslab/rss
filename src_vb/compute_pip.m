function [P1,P2] = compute_pip(segs_file, logw, alpha)
% USAGE: compute the posterior statistics P1 and P2 given a segment file and variational results
% INPUT:
%	segs_file: the mat file that contains genomic segment info, string
%	logw: log-importance weight for each combination of hyper-parameters, num_para by 1
%	alpha: variational estimates of the posterior inclusion probabilities, num_snps by num_para
% OUTPUT:
%       P1: the posterior inclusion probability of at least 1 SNP in the segment, num_segs by 1
%       P2: the posterior inclusion probability of at least 2 SNP in the segment, num_segs by 1

  % Load pre-defined genomic segments.
  fprintf('Loading pre-defined genomic segments.\n');
  seginfo = matfile(segs_file);
  Aseg    = seginfo.Aseg;  

  % For each segment, compute the posterior statistics P1
  % and P2 averaged over settings of the hyperparameters.
  fprintf('Compiling posterior results for each segment.\n'); 
  w        = normalizelogweights(logw);
  [P1, P2] = compilesegstats(Aseg,alpha,w);

end 

% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/normalizelogweights.m
% NORMALIZELOGWEIGHTS(LOGW) takes as input an array of unnormalized
% log-importance weights LOGW and returns normalized importance weights such
% that the sum of the normalized importance weights is equal to one.
function w = normalizelogweights (logw)

  % We guard against underflow or overflow by adjusting the log-importance
  % weights so that the largest importance weight is one.
  c = max(logw(:));
  w = exp(logw - c);

  % Normalize the importance weights.
  w = w / sum(w(:));  
end
  
% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/bmapathway/blob/master/MATLAB/results/compilesegstats.m
% For each segment, compute the posterior statistics P1 and P2 averaged over
% settings of the hyperparameters.
function [P1, P2] = compilesegstats (A, PIP, w)

  % Get the number of SNPs (p) and the number of segments (n).
  [~, n] = size(A);

  % These vectors will contain the relevant posterior statistics for each
  % segment.
  P1 = zeros(n,1);
  P2 = zeros(n,1);

  % Repeat for each segment.
  for i = 1:n

    % Get the SNPs in the segment.
    snps = find(A(:,i));
      
    % Compute, for each setting of the hyperparameters, the posterior
    % probability that no SNPs in the segment are included in the model,
    % and the posterior probability that exactly one SNP in the segment is
    % included in the model.
    p0 = computeprob0(PIP(snps,:));
    p1 = computeprob1(PIP(snps,:));
      
    % Compute the posterior probability that at least 1 SNP in the segment is
    % included in the model (P1), and the posterior probability that at
    % least 2 SNPs in the segment are included in the model (P2), using a
    % simple numerical approximation to average over settings of the
    % hyperparameters.
    P1(i) = dot(w,1 - p0);
    P2(i) = dot(w,1 - p1 - p0);
  end
end

% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/bmapathway/blob/master/MATLAB/results/computeprob0.m
% COMPUTEPROB0(PIP) returns the posterior probability that no variable is
% included in the model, assuming that the regression coefficients are
% independent of each other a posteriori. Each row of PIP corresponds to a
% variable, and each column corresponds to a different setting. The return
% value is a column vector equal to the number of columns of PIP.
function p0 = computeprob0 (PIP)
  p0 = prod(1 - PIP,1)';
end

% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/bmapathway/blob/master/MATLAB/results/computeprob1.m
% COMPUTEPROB1(PIP) returns the posterior probability that exactly one SNP
% is included in the model, assuming that the regression coefficients are
% independent of each other a posteriori. Each row of PIP corresponds to a
% variable, and each column corresponds to a different setting. The return
% value is a column vector equal to the number of columns of PIP.
function p1 = computeprob1 (PIP)
  p1 = computeprob0(PIP) .* sumrows(PIP./(1-PIP+eps))';
end

% Author: Peter Carbonetto 
% Source: https://github.com/pcarbo/bmapathway/blob/master/MATLAB/misc/sumrows.m
% SUMROWS(X) returns a column vector in which each entry is the sum over
% each column of X. This is the same as SUM(X,1).
function y = sumrows (x)
  y = sum(x,1);
end
