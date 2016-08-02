function bandwidth = find_bandwidth(A)
% USAGE: find the bandwidth of a symmetric band matrix
% SOURCE: Numerical Computing with MATLAB: Revised Reprint pp 72
  [i, j]    = find(A);
  bandwidth = max(abs(i-j));
end
