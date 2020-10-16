# InTACT `v1.0.0`

`InTACT` is a command line tool for testing gene-trait associations by integrating GWAS summary statistics with gene expression levels from multiple tissues. It returns *P*-values for both multi- and single-tissue association analyses.

InTACT has three compelling features that make it potentially useful in identifying novel gene-trait associations:

- It is very fast and maintains higher power than benchmark methods, that is, the burden sum test (i.e., the standard TWAS) and MultiXcan.
- It effectively integrates several tissue-specific and multi-tissue gene-trait association tests via Cauchy-distribution transformation.
- Its *P*-value is accurate by approximating it through the standard Cauchy distribution, especially when the *P*-values are small.

## Getting Started (Requirements)

First, clone this repository of `InTACT`:

```bash
git clone https://github.com/bae-y/InTACT.git
cd InTACT
```

Before running `InTACT`, you should download the following requirements into your `InTACT/` directory :

- `plink` from [PLINK1.9](https://www.cog-genomics.org/plink/1.9/),
- LD reference data (such as 1000 Genomes European ancestry data),
- the summary statistics data (such as from IGAP1) , and
- a file listing all the information for the weights files (such as TWAS/FUSION).

The file listing information for the weights files must have columns `WGT`, `ID`, `CHR`, `GENE`, `P0`, `P1`, `PANEL`, which will look like this:

           PANEL                                     WGT                 ID     GENE CHR       P0       P1
    136327 Liver Liver/Liver.ENSG00000160305.13.wgt.RDat ENSG00000160305.13    DIP2A  21 47878812 47989926
    143844  Lung    Lung/Lung.ENSG00000228314.1.wgt.RDat  ENSG00000228314.1 CYP4F29P  21 15215454 15220685
Also, please make sure that your R already has these libraries:

```R
R
install.packages(c('optparse','RColorBrewer','Rcpp','RcppArmadillo','bigmemory','mvtnorm','data.table'))
devtools::install_github("gabraham/plink2R/plink2R")
```

## Performing the InTACT

Once we have all the requirements, we can run `InTACT` with the following example command:

```bash
Rscript InTACT.R \
--sumstats IGAP1.rds \
--out out \
--weights TWAS_WEIGHTS.rds \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr 1000G.EUR.MultiTissue.CHR \
--chr 21
```

Since `InTACT` analyzes one chromosome at a time, chromosome number is a required input.

## Reading the Results Files

For each chromosome, `InTACT` returns two files:

1. a result table of single-tissue association tests and
2. a result table of multi-tissue association tests.

The first table shows single-tissue association analyses for every single tissue. The columns `TWAS.P`, `SSU.P`, `ACAT.P` indicate the *P*-values of the burden sum test (i.e., the standard TWAS), the burden score test (i.e., sum of squared U test), and the Cauchy combination test, respectively. `Tk.P` represents the *P*-value of tissue *k*, after combining the three tissue-specific tests.

                             PANEL CHR                 ID ...    TWAS.P     SSU.P     ACAT.P      Tk.P
    7949      Adipose_Subcutaneous  21 ENSG00000160305.13 ... 0.4863835 0.7773842  0.8560749 0.7605017
    14172 Adipose_Visceral_Omentum  21 ENSG00000160305.13 ... 0.7925701 1.0000000  0.8213858 1.0000000
    ...
    25074             Artery_Aorta  21  ENSG00000228314.1 ... 0.8551422 0.1066685 0.11642063 0.2283266
    36475            Artery_Tibial  21  ENSG00000228314.1 ... 0.6329966 0.2443968 0.25933623 0.3494701
    ...
The second table shows the results from gene-trait association tests for `NPANEL` multiple tissues. The columns `T1m.P`, `MultiXcan.P`, `InTACT.P` stand for the *P*-values of the Cauchy-based joint test, MultiXcan, and InTACT, respectively. To compare their performance with the single-tissue test, we also provide the results for the standard TWAS (`BEST.TWAS.P`) to which we ought to apply the stringent Bonferroni-corrected significance threshold in further analysis.

        CHR                 ID ... NPANEL ... BEST.TWAS.P     T1m.P MultiXcan.P  InTACT.P ...    
    87   21 ENSG00000160305.13 ...     28 ...   0.1266740 1.0000000   0.4540494 1.0000000 ...
    120  21  ENSG00000228314.1 ...     22 ...   0.1512981 0.5228292   0.1193236 0.2166768 ...
## Citation

If you use the software, please cite

Bae, Y. E., Bradley, J., Wu, L., and Wu, C. (?) An integrative approach for identifying gene-trait associations in multi-tissue transcriptome-wide association studies. Journal. doi.

## Author

Ye Eun Bae (Department of Statistics, Florida State University)