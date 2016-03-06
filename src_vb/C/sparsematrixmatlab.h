#ifndef INCLUDE_DOUBLEMATRIXMATLAB
#define INCLUDE_DOUBLEMATRIXMATLAB

#include "types.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Type definitions.
// -----------------------------------------------------------------
// A sparse matrix with double-precision floating-point entries.
typedef struct {
  Size    nr;     // Number of rows.
  Size    nc;     // Number of columns.
  Size    nz;     // Number of non-zero entries.
  double* elems;  // Entries of matrix.
  Index*  ir;     // Row indicators.
  Index*  jc;     // Column indicators.
} SparseMatrix;

// Function declarations.
// -----------------------------------------------------------------
// Get a sparse matrix from a MATLAB array. 
SparseMatrix getSparseMatrix (const mxArray* ptr);

// Copy column j of matrix X. The entries in a single column of the
// matrix are assumed to be stored consecutively in memory. Input n
// is the number of rows in the matrix.
void copySparseColumn (double* y, Index k, const double* X, Index* Ir, Index* Jc, Size n);

#endif
