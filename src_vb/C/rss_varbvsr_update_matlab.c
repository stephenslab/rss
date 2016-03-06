// For a description of this C code, see rss_varbvsr_update.m.
#include "types.h"
#include "vectorops.h"
#include "doublevectormatlab.h"
#include "sparsematrixmatlab.h"
#include "rssvarbvsr.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// MEX-file gateway routine. Note that varbvsupdate.m checks the
// inputs, so we do not have to do it here.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // Get inputs.
  const SparseMatrix SiRiS	= getSparseMatrix(prhs[0]);
  const double       sigma_beta	= *mxGetPr(prhs[1]);
  const DoubleVector logodds	= getDoubleVector(prhs[2]);
  const DoubleVector betahat	= getDoubleVector(prhs[3]);
  const DoubleVector se		= getDoubleVector(prhs[4]);
  const DoubleVector alpha0	= getDoubleVector(prhs[5]);
  const DoubleVector mu0	= getDoubleVector(prhs[6]);
  const DoubleVector SiRiSr0	= getDoubleVector(prhs[7]);
  const DoubleVector I		= getDoubleVector(prhs[8]);

  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const Size p = betahat.n;
  const Size m = I.n;
  Index* Ir    = SiRiS.ir;
  Index* Jc    = SiRiS.jc;

  // Initialize outputs.
  DoubleVector alpha  = createMatlabVector(p,&plhs[0]);
  DoubleVector mu     = createMatlabVector(p,&plhs[1]);
  DoubleVector SiRiSr = createMatlabVector(p,&plhs[2]);

  copyDoubleVector(alpha0,alpha);
  copyDoubleVector(mu0,mu);
  copyDoubleVector(SiRiSr0,SiRiSr);

  // Store a single column of matrix inv(S)*R*inv(S).
  double* SiRiS_snp = malloc(sizeof(double)*p);

  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  for (Index j = 0; j < m; j++) {
    Index k = (Index) I.elems[j];

    // Copy the kth column of matrix inv(S)*R*inv(S).
    copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);

    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr.elems[k];

    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat.elems[k], se.elems[k], sigma_beta, 
		       SiRiS_snp, SiRiSr.elems, SiRiSr_snp, 
		       logodds.elems[k], alpha.elems+k, mu.elems+k, p);
  }

  // Free the dynamically allocated memory.
  free(SiRiS_snp);
}
