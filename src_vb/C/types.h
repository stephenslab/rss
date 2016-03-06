// Author: Peter Carbonetto
// Source: https://github.com/pcarbo/varbvs/blob/master/R/varbvs/src/types.h
#ifndef INCLUDE_TYPES
#define INCLUDE_TYPES

#ifdef MATLAB_MEX_FILE

// These type definitions are used to build the C routines for MATLAB.
#include "matrix.h"

typedef mwSize  Size;
typedef mwIndex Index;
typedef float   MatrixElem;

#else

// These type definitions are used to build the C routines for R. 
typedef int    Size;
typedef int    Index;
typedef double MatrixElem;

#endif

#endif
