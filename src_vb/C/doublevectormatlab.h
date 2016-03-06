// Author: Peter Carbonetto 
// Source: https://github.com/pcarbo/varbvs/blob/master/MATLAB/C/doublevectormatlab.h
#ifndef INCLUDE_DOUBLEVECTORMATLAB
#define INCLUDE_DOUBLEVECTORMATLAB

#include "types.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Type definitions.
// -----------------------------------------------------------------
// A vector with double-precision floating-point entries.
typedef struct {
  Size    n;      // Size of vector.
  double* elems;  // Vector entries.
} DoubleVector;

// Function declarations.
// -----------------------------------------------------------------
// Get a double-precision floating-point vector from a MATLAB array.
DoubleVector getDoubleVector (const mxArray* ptr);

// Create a column vector in MATLAB array.
DoubleVector createMatlabVector (Size n, mxArray** ptr);

// Copy entries of one vector to another vector.
void copyDoubleVector (const DoubleVector source, DoubleVector dest);

#endif
