---
layout: default
---

## Compute scaled population recombination rate

The following example illustrates how to compute the "scaled population recombination rate" using HapMap genetic map.

We take six SNPs in chromosome 22 from HapMap CEU Phase 2 sample, and we get the effective diploid population size `Ne=11418` from [IMPUTE](https://mathgen.stats.ox.ac.uk/impute/impute_v1.html) software document.

```r
   position             COMBINED_rate(cM/Mb)                                Genetic_Map(cM)
1: 14431347             8.096992                                            0.00000000
2: 14432618             8.131520                                            0.01029128
3: 14433624             8.131967                                            0.01847159
4: 14433659             8.132625                                            0.01875620
5: 14433758             8.129606                                            0.01956133
6: 14434713             8.024772                                            0.02732511
```

To compute recombination rates from any two SNPs, we use a formula from [Li and Stephens (2003)](https://www.ncbi.nlm.nih.gov/pubmed/14704198):

> rho = 4 * effective diploid population size * genetic distance

Eg. 1, for SNP 2 and 3: rho_23 = 4 * 11418 * (0.01847159-0.01029128)/100;

Eg. 2, for SNP 2 and 6: rho_26 = 4 * 11418 * (0.02732511-0.01029128)/100.

Next, as a validation, we compare our calculations with the results provided by [`BLIMP`](http://stephenslab.uchicago.edu/software_pages/blimp/index.html) ([Wen and Stephens, 2010](https://www.ncbi.nlm.nih.gov/pubmed/21479081)).

The file `rmb.ceu.ch22` lists recombination rate between all adjacent markers in the panel. 

```r
library(data.table)
genetic.map.chr22 <- data.table::fread("genetic_map_chr22.txt")
legend.ceu.chr22 <- data.table::fread("legend.ceu.ch22", header = T)
rmb.ceu.chr22 <- data.table::fread("rmb.ceu.ch22")

# locate all the snps specified by the legend file
selected.pos <- intersect(genetic.map.chr22$position, legend.ceu.chr22$position)
map.index <- which(genetic.map.chr22$position %in% selected.pos)
leg.index <- which(legend.ceu.chr22$position %in% selected.pos)

# check the location is correct
location.check = prod(genetic.map.chr22$position[map.index] == legend.ceu.chr22$position[leg.index])
if (location.check != 1) stop('locate the snps wrongly')

# calculate the shrinking coefficient: exp(-rho_{ij}/(2*m)) for adjacent (i,j)
# general formula: rho = 4 * effective population size * genetic distance
# hapmap phase 2 ceu: effective population size (Ne) is 11418
#                     haplotype size (m) is 120
num.snp <- length(map.index);
ro.coef <- rep(0, num.snp);
ro.coef[1] = 1;
map.selected <- genetic.map.chr22[map.index, ]
rcb.selected <- map.selected[[3]]
for (i in 2:num.snp){
  ro.coef[i] = exp(-4*11418*(rcb.selected[i]/100-rcb.selected[i-1]/100)/120)
}

# compare the result with the one provided in BLIMP (Wen and Stephens, 2010)
# make sure rss compute rho in the same way as BLIMP
results <- matrix(0, nrow=num.snp, ncol=2)
results[, 1] <- unlist(rmb.ceu.chr22)
results[, 2] <- ro.coef
```

```r
> head(results)
          [,1]      [,2]
[1,] 1.0000000 1.0000000
[2,] 0.9615886 0.9615886
[3,] 0.9693454 0.9693454
[4,] 0.9989173 0.9989173
[5,] 0.9969404 0.9969404
[6,] 0.9708834 0.9708834

> tail(results)
              [,1]      [,2]
[34021,] 0.9963701 0.9963701
[34022,] 0.9909937 0.9909937
[34023,] 0.9986032 0.9986032
[34024,] 0.9879030 0.9879030
[34025,] 0.9962124 0.9962124
[34026,] 0.9847357 0.9847357

> max(abs(results[,1]-results[, 2]))
[1] 4.999806e-08
```
