### Input Data Formats

#### GWAS summary statistics

This section is modified from Box 1 of [Winkler *et al.* (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24762786).   

The following columns are **required** for any RSS analysis:

- `snp`: identifier of genetic variant, character string such as `rs12498742`; 
- `chr`: chromosome number of genetic variant, a character string from `chr1`, `chr2` ... `chr22`, `chrX`, `chrY`;
- `pos`: physical position, in base pair, of genetic variant, an integer;
- `a1`: allele associated with the trait (corresponding to the regression coefficient), a single upper case character `A`, `C`, `G` or `T`;
- `a2`: the other (non-effect) allele, a single upper case character `A`, `C`, `G` or `T`;
- `betahat`: estimated effect size of genetic variant under the single-marker model, numeric;
- `se`: estimated standard error of `betahat`, numeric.    

The following columns are optional, but they can be very helpful for sanity checks:

- `strand`: strand on which the alleles are reported, a single character `-` or `+`;
- `n`: number of individuals analyzed (a.k.a. sample size) for the genetic variant, an integer;
- `maf`: minor allele frequency, numeric between 0 and 1;
- `p`: p-value of genetic variant association, numeric between 0 and 1;
- `info`: other information about genetic variants.

It is very **important** to make sure that `[a1, betahat, se]` are perfectly aligned. Here is a toy example.

> Consider two SNPs (`rs1`, `rs2`) and four individuals (`i1`, `i2`, `i3`, `i4`):
> ```
> IND, i1, i2, i3, i4
> rs1, AT, TT, AT, AA
> rs2, CG, CC, GG, GC 
> ```
> If the effect alleles (`a1`) of these two SNPs are `A` and `G` respectively, then the genotype data of `rs1` are `X[, 1]=[1, 0, 1, 2]`, and the genotype data of `rs2` are `X[, 2]=[1, 0, 2, 1]`. Further, the single-SNP summary statistics of `rs1` and `rs2` are given by:
> ```
> (betahat[1], se[1]) <- single.SNP.model(y, X[, 1])
> (betahat[2], se[2]) <- single.SNP.model(y, X[, 2])
> ```

Finally, when providing `chr` and `pos` columns, please explicitly specify the [assembly releases and versions](https://genome.ucsc.edu/FAQ/FAQreleases.html) of human genome. If 1000 Genomes Project Phase 3 data are used to estimate LD, please provide `chr` and `pos` columns based on UCSC hg19/GRCh37.   

#### LD matrix estimates

The LD estimates are often derived from the phased haplotype data from [1000 Genomes Project Phase 3 data](http://www.internationalgenome.org/data). Because the 1000 Genomes data are publicly available, the LD estimates only require the list of genetic variants, their physical positions and their effect alleles (i.e. `[snp, chr, pos, a1]` from the summary statistics file).

If there are some internal genotype data that can be used to estimate LD matrix, please organize the genotype data in the same format as 1000 Genomes Phase 3 data:<br>[`ALL.chr${chrid}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)<br>Again, make sure that the physical positions and effect alleles of the internal genotype data are consistent with `[chr, pos, a1]` provided in the GWAS summary statistics file.

#### Genomic annotations

The most statistician-friendly format of genomic annotation data might look like this:

```
 snp   chr   pos ann1 ann2 ann3
 rs1  chr2 52877    0    0    0
 rs2  chr1 50670    0    1    0
 rs3 chr14   854    0    1    1
 rs4  chr4 99620    1    1    1
 rs5 chr16 71537    0    0    0
 rs6 chr22 39741    0    0    0
 rs7  chr6 89331    1    0    0
```
where `ann1`, `ann2` and `ann3` are three types of annotations, `1` indicates that SNP is annotated and `0` otherwise.

Alternatively, a list of annotated SNPs can be saved as a separate file. For example:

```
> cat ann3.txt
 snp   chr   pos
 rs3 chr14   854
 rs4  chr4 99620

> cat ann2.txt
 snp   chr   pos
 rs2  chr1 50670
 rs3 chr14   854
 rs4  chr4 99620

> cat ann1.txt
 snp  chr   pos
 rs4 chr4 99620
 rs7 chr6 89331
```
Sometimes the annotations are based on genes or genomic regions (e.g. biological pathways). For these annotations, it is easier to provide a list of annotated regions:

```
 ensembl_gene_id chromosome_name start_position end_position
 ENSG00000000938               1       27938575     27961788
 ENSG00000008438              19       46522411     46526323
 ENSG00000008516              16        3096682      3110727
 ENSG00000066336              11       47376411     47400127
 ENSG00000077984              20       24929866     24940564
 ENSG00000085265               9      137801431    137809809
```

For all these annotation files, please make sure that the physical positions (`[snp, chr, pos]` or `[chromosome_name, start_position, end_position]`) are **consistent** with `[snp, chr, pos]` in the summary statistics file.                   