function B = band_storage(A, p)
% USAGE: this function stores symmetric banded matrix A in a 
% compact form bA in such a way that only the main diagonal,
% and the nonzero superdiagonals are stored. The first column 
% of bA corresponds to the main diagonal of A and the subsequent 
% columns of bA correspond to superdiagonals of A.
% INPUT: 
%	p: upper or lower bandwidth
%	A: symmetric matrix
% OUTPUT: 
%	B: the banded storage of A
% SOURCE: http://www.mathworks.com/matlabcentral/fileexchange/31131-efficient-cholesky-decomposition-of-symmetric-banded-matrix/content/SymmetricBandedCholesky.m

  dim=size(A);
  n = dim(1);
  B = zeros(p+1, n);
  for i=1:n % column
    if i<=n-p
      for j=i:p+i
      B(j-i+1, i) = A(j, i);
      end
    else
      for j=i:n
      B(j-i+1, i) = A(j, i);
      end
    end
  end
end
