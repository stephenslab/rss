// Author: Peter Carbonetto
// Source: https://github.com/pcarbo/varbvs/blob/master/R/varbvs/src/sigmoid.h
#ifndef INCLUDE_SIGMOID
#define INCLUDE_SIGMOID

// Function declarations.
// -----------------------------------------------------------------
// Compute log(1 + exp(x)) in a numerically stable manner.
double logpexp (double x);

// Return the sigmoid function at x.
double sigmoid (double x);

// Return the logarithm of the sigmoid function at x. Computation is
// performed in a numerically stable manner.
double logsigmoid (double x);

#endif
