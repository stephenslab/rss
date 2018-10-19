# This file computes the enrichment likelihood ratios
# for a given list of genomic regions (e.g. gene set)
# using GWAS summary statistics of 31 complex traits.
# We developed this simple, fast calculation as a sanity
# check for our Bayesian model-based enrichment method.
# For more details, please refer to the following article:
#
# Zhu and Stephens (2018), https://www.nature.com/articles/s41467-018-06805-x.
#
# If you find this method and/or script useful for your
# research work, please kindly cite the article above.
#
# Before using this script, you need to download two folders.
# For people with access to Uchicago Midway, please find
# gwas_sumstat <- "/scratch/midway2/xiangzhu/gwas_sumstat/"
# baseline_res <- "/project/mstephens/test_rss/results/ashlrt_baseline/"
# For people with access to Stanford Sherlock, please find
# gwas_sumstat <- "/scratch/PI/whwong/xiangzhu/gwas_sumstat/"
# baseline_res <- "/scratch/PI/whwong/xiangzhu/baseline_results/"
#
# When preparing for the list of genomic regions, please use
# the BED format without header; see below for an example.
#
# $ head -n 3 kidney_P0_enhancr_U6.bed
# chr8	67364174	67364622
# chr6	52447909	52450568
# chr6	72121430	72121779
#
# In addition, please ensure that the BED file is based on hg19.
#
# After all input data are downloaded and correctly specified,
# please type the following command in a R console:
# > source('ash_lrt_31traits.R')
# The analysis results will be written to a plain text file.

# Most users only need to modify the following lines.

# Specify paths to input data files.
gwas_sumstat <- "/Users/xiangzhu/Dropbox/rss/Data/gwas_gsea/viz_gwas_zscore/input_data/"
baseline_res <- "/Users/xiangzhu/Dropbox/rss/Data/gwas_gsea/compute_ash_lrt/baseline_results/"
region_file <- "/Users/xiangzhu/Downloads/mouse_development_enhancer_hg19/kidney_P0_enhancr_U6.bed"

# Sometimes one may want to assign nearby SNPs to a region.
# Please adjust the following variable to set window size (bp).
window_size <- 0

# Save the analysis results as a plain text file.
out_path <- "test.txt"

# PLEASE DO NOT MODIFY ANY CODE BELOW THIS LINE.

# Load R packages.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.matlab))
suppressPackageStartupMessages(library(ashr))

#' Given a list of genomic regions (BED), map SNPs to regions.
#'
#' @param region_file path to the given BED file
#' @param chr a vector of chromosome ID for each SNP
#' @param pos a vector of physical position (base pair) for each SNP
#' @param window_size: length of the "window" for assigning SNPs to regions
#' @return Locations of SNPs mapped to the BED file
#'
snp_to_region <- function(region_file, chr, pos, window_size) {

  # Load and clean the BED file.
  bed_df <- data.table::fread(region_file)
  names(bed_df) <- c("chr","start_pos","end_pos")
  bed_df <- bed_df[bed_df$chr %in% paste0("chr",1:22), ]
  bed_df <- bed_df[ ,c("chr","start_pos","end_pos")]
  bed_df$chr <- as.numeric(gsub("chr","",bed_df$chr))

  # Create an empty list that each region will take one entry.
  num_region <- dim(bed_df)[1]
  snps_list <- list()

  # repeat for each region
  for (i in 1:num_region) {
    valid_region <- (bed_df$start_pos[i] > 0)
    valid_region <- valid_region & (bed_df$end_pos[i] >= bed_df$start_pos[i])
    if (valid_region){
      match_chr <- (chr == bed_df$chr[i])
      match_start <- (pos > bed_df$start_pos[i] - window_size)
      match_end <- (pos < bed_df$end_pos[i] + window_size)
      snps_list[[i]] <- which(match_chr & match_start & match_end)
    }
  }

  # Note that one SNP may be assigned to multiple nearby regions.
  # The following line ensures that each SNP occurs at most once.
  snps <- sort(unique(unlist(snps_list)))

  return(snps)
}

#' Given a list of genomic regions (BED), compute the enrichment
#' likelihood ratio using GWAS summary statistics of a given trait.
#'
#' @param region_file path to the given BED file
#' @param sumstat_file path to the GWAS summary statistics
#' @param baseline_file path to the fitted baseline model
#' @param window_size: length of the "window" for assigning SNPs to regions
#' @return log 10 likelihood ratio
#'
compute_ash_lrt <- function(region_file, sumstat_file, baseline_file, window_size=0) {

  # Load the whole genome GWAS summary statistics.
  sumstat <- R.matlab::readMat(sumstat_file)
  betahat <- c(sumstat$betahat)
  se <- c(sumstat$se)
  chr <- c(sumstat$chr)
  pos <- c(sumstat$pos)

  # Map whole genome SNPs to the BED file.
  snps <- snp_to_region(region_file, chr, pos, window_size)

  # Load fitted baseline model based on whole genome data.
  fitted_g <- readRDS(baseline_file)

  # Fit the baseline model only on SNPs mapped to BED file.
  a2 <- ashr::ash(betahat=betahat[snps], sebetahat=se[snps],
                  mixcompdist="halfuniform", method="shrink")

  # Fit the enrichment model only on SNPs mapped to BED file.
  a3 <- ashr::ash(betahat=betahat[snps], sebetahat=se[snps],
                  mixcompdist="halfuniform", method="shrink", fixg=T, g=fitted_g)

  # Compute log likelihood ratio statistics.
  logLR <- a2$logLR - a3$logLR
  log10_LR <- logLR / log(10)

  if (is.na(log10_LR)) {
    log10_LR <- Inf
  }

  return(log10_LR)
}

# Through the following for loop, we can obtain the enrichment
# likelihood ratios of 31 complex traits analyzed in our article.

sumstat_files <- list.files(path=gwas_sumstat, pattern="*.mat", full.names=TRUE)
num_trait <- length(sumstat_files)

trait_name = log10_lr = NULL

for (j in 1:num_trait) {

  # Specify the GWAS summary statistics file.
  sumstat_file <- sumstat_files[j]

  # Specify the baseline model fitting file.
  sumstat <- R.matlab::readMat(sumstat_file)
  trait_name[j] <- paste0(sumstat$name, sumstat$time)
  baseline_file <- paste0(baseline_res, trait_name[j], ".rds")

  # Perform the simple likelihood ratio calculation.
  log10_lr[j] <- compute_ash_lrt(region_file, sumstat_file, baseline_file, window_size=0)

  cat("Finish analyzing ", trait_name[j], "... \n")
}

# Save the analysis results in a plain text file.
out_df <- data.frame(trait_name=trait_name, log10_lr=log10_lr)
write.table(out_df, file=out_path, col.names=T, row.names=F, sep='\t', quote=F)


