function chrpar = partition_genome(betahat)
% USAGE: get the starting and ending positions of each chromosome
% INPUT:
%	betahat: single-SNP effect estimates, 22 by 1 cell
% OUTPUT:
%	chrpar: starting and ending positions, 22 by 2 array

  if length(betahat) ~= 1
    % get the length of each chromosome
    chr_length = cellfun(@length, betahat);

    % get the ending positions of each chromosome
    chr_end = cumsum(chr_length);

    % get the starting positions of each chromosome
    chr_start = [1; chr_end(1:end-1)+1];

    % combine results
    chrpar = [chr_start chr_end];
 
  else
    chrpar = [1 length(cell2mat(betahat))];
  end

end
