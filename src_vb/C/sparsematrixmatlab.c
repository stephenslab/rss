#include "sparsematrixmatlab.h"

// Function definitions.
// -----------------------------------------------------------------

// Get a sparse matrix from a MATLAB array.
SparseMatrix getSparseMatrix (const mxArray* ptr){
  SparseMatrix result;
  result.nr    = mxGetM(ptr);
  result.nc    = mxGetN(ptr);
  result.nz    = mxGetNzmax(ptr);
  result.elems = (double*) mxGetPr(ptr);
  result.ir    = (Index*) mxGetIr(ptr);
  result.jc    = (Index*) mxGetJc(ptr);
  return result;
}

// Copy column k of sparse matrix using the triplet (X, Ir, Jc)
void copySparseColumn (double* y, Index k, const double* X, Index* Ir, Index* Jc, Size n) {
  for (Index i = 0; i != n; ++i) {
	y[i] = 0.0;
  }
  for (Index l = Jc[k]; l != Jc[k+1]; ++l) {
	y[Ir[l]] = X[l];
  } 
}

