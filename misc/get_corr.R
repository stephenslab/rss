# This script is an R translation of `get_corr.m`.
# Source: https://github.com/stephenslab/rss/blob/master/misc/get_corr.m

#' compute the shrinkage estimator of LD matrix in Wen and Stephens (2010)
#' source: https://pubmed.ncbi.nlm.nih.gov/21479081
#'
#' @param m number of individuals in the reference panel, integer
#' @param Ne effective population size (diploid), integer
#' @param cummap cumulative genetic map in cM, numSNP by 1
#' @param Hpanel phased haplotypes or unphased genotypes from a reference panel, numIND by numSNP
#' @param cutoff hard threshold (e.g. 1e-3) that forces small entries to zero, scalar
#' @param isgeno TRUE if Hpanel is an unphased genotype matrix, logical
#' @param negcm TRUE if there is negative genetic distance, logical
#'
#' @return a list with RHat being estimated LD matrix and SigHat being estimated covariance matrix
#'
get_corr <- function(m, Ne, cummap, Hpanel, cutoff, isgeno=FALSE, negcm=FALSE) {

  # theta is related to mutation suggested by Li and Stephens (2003)
  # see Equation 2.8 of Wen and Stephens (2010)
  nmsum <- sum(1 / (1:(2*m-1)))
  theta <- (1/nmsum) / (2*m + 1/nmsum)

  # S is obtained from Sigma_panel by shrinking off-diagonal entries toward 0
  # see Equation 2.7 of Wen and Stephens (2010)

  # Scenario 1: Hpanel is a phased haplotype matrix
  S <- var(Hpanel)

  # Scenario 2: Hpanel is an unphased genotype matrix
  # see Section 2.4 for justification of "0.5"
  if (isgeno) S <- 0.5*S

  # obtain the upper triangular portion of S
  S[lower.tri(S, diag=FALSE)] <- 0

  # compute shrinkage factors on and above the 1st diagonal
  numSNP <- dim(S)[1]

  for (i in 1:(numSNP-1)) {
    for (j in (i+1):numSNP) {

      # compute genetic distance between SNPs i and j
      genetic_dist <- cummap[j] - cummap[i]

      # deal with potential negative distance value
      # which may occur after lifting genome builds
      if (genetic_dist<0 & !negcm) stop("Negative recombination rate is not allowed.")
      if (genetic_dist<0 & negcm) genetic_dist <- abs(genetic_dist)

      # compute scaled population recombination rate
      # see https://stephenslab.github.io/rss/recombination.html
      rho <- 4 * Ne * genetic_dist / 100

      # compute the shrinkage factor based on recombination rate
      # see Equation 2.7 of Wen and Stephens (2010)
      shrinkage <- exp(-rho/(2*m))

      # apply hard thresholding to obtain sparse and banded estimator
      if (shrinkage < cutoff) shrinkage <- 0

      # shrink the sample covariance in the reference panel
      # see Equation 2.7 of Wen and Stephens (2010)
      S[i,j] <- shrinkage * S[i,j]

    }
  }

  # copy the upper half (NO diagonal) to the lower half
  S_upper <- S
  S_lower <- t(S); diag(S_lower) <- 0

  S <- S_upper + S_lower
  rm(S_upper, S_lower)

  # SigHat is derived from Li and Stephens (2003)
  # see Equation 2.6 of Wen and Stephens (2010)
  SigHat <- (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * diag(x=1,nrow=numSNP,ncol=numSNP,names=FALSE)

  # convert covariance to correlation
  RHat <- cov2cor(SigHat)

  # output results
  return(list(RHat=RHat,SigHat=SigHat))

}
