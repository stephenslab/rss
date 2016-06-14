#!/bin/bash

chrid="$1"
dataname="$2"

vcf1000g_path="$HOME/Data/1000g/20130502"
vcfchr="${vcf1000g_path}/ALL.chr${chrid}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# the list of individuals to include in analysis
# each individual ID (as defined in the VCF headerline) should be included on a separate line
# no header line is expected
eursam="${vcf1000g_path}/EUR.samples.1000g.20130502"

# the list of SNP IDs, with one SNP per line
# no header line is expected
snplist_path="$HOME/Data/${dataname}_analysis/analyzed_snp"
snplist="${snplist_path}/chr${chrid}.snp.txt"

analyzed_1000g_path="$HOME/Data/${dataname}_analysis/1000genomes_vcf"
analyzed1000g="${analyzed_1000g_path}/chr${chrid}.analyzed"

# create the sub-folder for vcf files if it does not exist
if [ ! -d "$analyzed_1000g_path" ]; then
  mkdir "$analyzed_1000g_path"
fi

# generate a new vcf file after applying the following filtering options:
# --keep: only include the EUR individuals
# --snps: only include SNPs given in ~/Data/scz2_analysis/analyzed_snp/chr*.snp.txt
# --min/max-alleles: only include bi-allelic sites
vcftools --gzvcf "$vcfchr" --keep "$eursam" --snps "$snplist" --min-alleles 2 --max-alleles 2 --recode --out "$analyzed1000g" 

# output phased haplotypes in IMPUTE reference-panel format ".impute.hap"
# IMPUTE format: one row per SNP and one column per haplotype (0 or 1)
# output files:
# - ".impute.hap": IMPUTE haplotype file
# - ".impute.hap.legend": IMPUTE legend file
# - ".impute.hap.indv": individuals included in the haplotype file
# NB: the last 2 columns of ".impute.hap.legend" (i.e. allele0 and allele1) 
#     specify the alleles underlying the 0/1 coding in the corresponding
# more about IMPUTE format:
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#input_options
vcftools --vcf "$analyzed1000g".recode.vcf --IMPUTE --out "$analyzed1000g"
