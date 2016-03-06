#include "sigmoid.h"
#include <math.h>

// Function definitions.
// -----------------------------------------------------------------
// Computes log(1 + exp(x)).
double logpexp (double x) {
  return (x >= 16) * x + (x < 16)  * log(1 + exp(x));
}

// Return the sigmoid function at x.
double sigmoid (double x) {
  return 1/(1 + exp(-x));
}

// Return the logarithm of the sigmoid function at x.
double logsigmoid (double x) {
  return -logpexp(-x);
}
