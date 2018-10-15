#include "rssvarbvsr.h"
#include "sigmoid.h"
#include "vectorops.h"
#include <math.h>

// This function runs a single coordinate ascent of variational update to fit RSS-BVSR.

void rss_varbvsr_update (double betahat, double se, double sigma_beta,
                         const double* SiRiS_snp, double* SiRiSr, double SiRiSr_snp,
                         double logodds, double* alpha, double* mu, Size p) {

  double se_square         = se * se;
  double sigma_beta_square = sigma_beta * sigma_beta;

  // Compute the variational estimate of the posterior variance.
  double sigma_square = (se_square * sigma_beta_square) / (se_square + sigma_beta_square);
  
  // Update the variational estimate of the posterior mean.
  double r = (*alpha) * (*mu);
  *mu      = sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp);
  
  // Update the variational estimate of the posterior inclusion probability.
  double SSR = (*mu) * (*mu) / sigma_square;
  *alpha     = sigmoid(logodds + 0.5 * (log(sigma_square/(sigma_beta_square)) + SSR));
  
  // Update SiRiSr = inv(S)*R*inv(S)*r.
  double r_new = (*alpha) * (*mu);
  add(SiRiSr, r_new - r, SiRiS_snp, p);
}
