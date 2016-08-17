% RANDINT   Generate a matrix of random integers.
%    RANDINT(A,B,M,N) generates an M x N matrix of random numbers
%    uniformly distributed in the set of integers between A and B.
function X = randint (a, b, m, n)
  if nargin < 3
    m = 1;
    n = 1;
  end
  U = rand(m,n);
  U = U .* (U < 1);
  X = floor(a + (b-a+1)*U);
