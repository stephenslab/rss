#ifndef INCLUDE_BVSR
#define INCLUDE_BVSR

#include "types.h"

// Function declarations.
// ------------------------------------------------------------------------------
// USAGE: execute a single coordinate ascent of variational update to fit RSS-BVSR model
// INPUT: betahat and se are the single-SNP effect estimate and its SE of one particular SNP
//	  sigma_beta and logodds are the hyperparameters
//	  SiRiS_snp is the corresponding column of SiRiS
//	  SiRiSr is the vector inv(S)*R*inv(S)*r that will be updated to reflect the change to alpha and mu
//	  SiRiSr_snp is the corresponding entry of SiRiSr
//	  alpha and mu are the variational parameters that will be updated 
//	  p is the total number of SNPs
void rss_varbvsr_update (double betahat, double se, double sigma_beta, const double* SiRiS_snp, double* SiRiSr, double SiRiSr_snp, double logodds, double* alpha, double* mu, Size p);

#endif
