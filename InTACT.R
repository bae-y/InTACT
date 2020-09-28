library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
library(data.table)

sourceCpp("twas_perm.cpp")
source("aSPUwscore.R")
source("dist_support.R")
source("twas_support.R")

suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

options_list <- list(
  make_option(c("--sumstats"), action = "store", default = NA, type = "character",
              help = "Summary statistics (must have SNP and Z column headers) [required]"),
  make_option(c("--out"), action = "store", default = NA, type = "character",
              help = "Name of output file [required]"),
  make_option(c("--weights"), action = "store", default = NA, type = "character",
              help = "File listing molecular weight (must have columns WGT, ID, CHR, GENE, P0, P1, PANEL) [required]"),
  make_option(c("--weights_dir"), action = "store", default = NA, type = "character",
              help = "Path to directory where weight files (WGT column) are stored [required]"),
  make_option(c("--ref_ld_chr"), action = "store", default = NA, type = "character",
              help = "Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option(c("--chr"), action = "store", default = NA, type = "character",
              help = "Chromosome to analyze [required]")
)
opt = parse_args(OptionParser(option_list=option_list))

##################################################
# Load in summary statistics and reference data ##
##################################################
tryCatch({
  sumstat.origin = readRDS(opt$sumstats)
}, warning = function(w){
  sumstat.origin = fread(opt$sumstats, header = T, data.table = F)
}, error = function(e){
  sumstat.origin = fread(opt$sumstats, header = T, data.table = F)
})

genos.input = read_plink(paste(opt$ref_ld_chr, opt$chr, sep = ""), impute = "avg")
genos.snp = genos.input$bim[, 2]

# Match summary data to input, record NA where summary data is missing
m = match(genos.snp, sumstat.origin$SNP)
sumstat.origin = sumstat.origin[m, ]


##################################################
# Load in list of weights                       ##
##################################################
tryCatch({
  wgtlist = readRDS(opt$weights)
}, warning = function(w){
  wgtlist = fread(opt$weights, header = T, data.table = F)
}, error = function(e){
  wgtlist = fread(opt$weights, header = T, data.table = F)
})

wgtlist = as.data.frame(wgtlist)

if ( sum(names(wgtlist) == "WGT") == 0 | 
     sum(names(wgtlist) == "ID") == 0  |
     sum(names(wgtlist) == "CHR") == 0 |
     sum(names(wgtlist) == "GENE") == 0|
     sum(names(wgtlist) == "P0") == 0  |
     sum(names(wgtlist) == "P1") == 0  ) {
  cat("ERROR : One or more columns (WGT, ID, CHR, GENE, P0, P1, PANEL) required in weights file\n")
  q()
}

wgtlist = wgtlist[wgtlist[, "CHR"] == opt$chr, ]
wgtlist = wgtlist[!duplicated(wgtlist), ]

tissue.panel = unique(wgtlist[, "PANEL"])

used.gene = unique(wgtlist[, "ID"])

G = length(used.gene)


##################################################
# Output will look like this:                   ##
##################################################
res.single.tissue = NULL

res.multi.tissue = as.data.frame(matrix(NA, G, 15))
colnames(res.multi.tissue) = c("CHR", "ID", "GENE", "P0", "P1",
                               "NPANEL", "BEST.PANEL", "BEST.NSNP", "BEST.TWAS.Z", "BEST.TWAS.P",
                               "T1m.P", "MultiXcan.P", "InTACT.P", "MOST.SIG.SNP", "RUNTIME")


##################################################
# For each gene in weights file:                ##
##################################################
for (gene.indx in 1:G) { 
  tryCatch({
    
    sumstat = sumstat.origin
    
    wgtlist1 = wgtlist[wgtlist[, "ID"] == used.gene[gene.indx], ]
    
    genos = genos.input
    
    P0.tmp = as.numeric(wgtlist1[1, "P0"]) - 1e6
    P1.tmp = as.numeric(wgtlist1[1, "P1"]) + 1e6
    genos$bim = genos$bim[genos$bim[, 4] > P0.tmp & genos$bim[, 4] < P1.tmp, ]
    genos$bed = genos$bed[ ,colnames(genos$bed) %in% genos$bim[,2], drop=F]
    
    sumstat = sumstat[sumstat[, 1] %in% genos$bim[, 2], ]
    
    # Match summary data to ref input, record NA where summary data is missing
    m <- match(genos$bim[, 2], sumstat$SNP)
    sum.missing <- is.na(m)
    sumstat <- sumstat[m, ]
    sumstat$SNP <- genos$bim[, 2]
    sumstat$A1[ sum.missing ] <- genos$bim[sum.missing, 5]
    sumstat$A2[ sum.missing ] <- genos$bim[sum.missing, 6]
    
    # find the most significant SNPs
    sumstat.tmp <- sumstat[sumstat[, 1] %in% genos$bim[, 2], ]
    Z <- max(abs(sumstat.tmp[, 4]), na.rm = T)
    most.sig.snp <- (1 - pnorm(Z)) * 2
    
    res.multi.tissue[gene.indx, "MOST.SIG.SNP"] <- most.sig.snp
    
    # QC / allele-flip the ref input and output
    qc <- allele.qc(sumstat$A1, sumstat$A2, genos$bim[, 5], genos$bim[, 6])
    
    # Flip Z-scores for mismatching alleles
    sumstat$Z[ qc$flip ] <- -1 * sumstat$Z[ qc$flip ]
    sumstat$A1[ qc$flip ] <- genos$bim[qc$flip, 5]
    sumstat$A2[ qc$flip ] <- genos$bim[qc$flip, 6]
    
    # Remove strand ambiguous SNPs (if any)
    if (sum(!qc$keep) > 0) {
      genos$bim <- genos$bim[qc$keep, ]
      genos$bed <- genos$bed[, qc$keep]
      sumstat <- sumstat[qc$keep, ]
    }
    
    m <- length(unique(wgtlist1[, "PANEL"]))
    
    start.time = proc.time()[3]
    if (m != 0) {
      wgtlist0 <- wgtlist1
      
      ##################################################
      # Step 1: Tissue-specific association test      ##
      ##################################################
      # TWAS, SSU, ACAT, Tk
      res <- SPrediXcan.sum(wgtlist0, opt, genos)
      
      res.save <- cbind(wgtlist0[, c("PANEL", "CHR", "ID", "GENE", "P0", "P1")], t(res$res.save))
      colnames(res.save) <- c("PANEL", "CHR", "ID", "GENE", "P0", "P1", "NSNP", "NON.ZERO.SNP","R2" ,"TWAS.Z", "TWAS.P", "SSU.P","ACAT.P", "Tk.P")
      
      res.save <- res.save[!is.na(res.save[, "TWAS.P"]),,drop=F]
      
      res.single.tissue <- rbind(res.single.tissue, res.save)
      
      wgt.mat <- res$wgt
      
      if (is.null(wgt.mat)) {
        m <- 0
      } else {
        tmp.colnames <- colnames(wgt.mat)
        tmp.indx <- colSums(abs(wgt.mat)) != 0
        tmp.colnames <- tmp.colnames[tmp.indx]
        
        wgt.mat <- wgt.mat[, tmp.indx, drop=F]
        
        colnames(wgt.mat) <- tmp.colnames #panel
        m <- dim(wgt.mat)[2]
      }
      
      ##################################################
      # Step 2: Multi-tissue association test         ##
      ##################################################       
      if (m > 1) {
        res.tmp = res.save[!is.na(res.save$Tk.P), , drop=F]
        res.tmp = res.tmp[res.tmp$PANEL %in% colnames(wgt.mat), , drop=F]
        Z.score = res.tmp[, "TWAS.Z"]
        
        # Cauchy-based joint test (T_{1:m})
        Tk.P = res.tmp$Tk.P
        T1m.P = ACAT(Tk.P, rep(1/length(Tk.P),length(Tk.P)))
        
        # MultiXcan
        used.weight = as.character(res.tmp[, "PANEL"])
        out.multi <- MultiXcan(wgt.mat, Z.score, used.weight) #using cutoff 30
        MultiXcan.P <- out.multi$multi.p
        
        # InTACT
        InTACT.P = ACAT(c(T1m.P, MultiXcan.P), c(1/2, 1/2))
      }
    }
    end.time <- proc.time()[3]
    
    basic.inf <- cbind(wgtlist1[1, c("CHR", "ID", "GENE", "P0", "P1")], m)
    res.multi.tissue[gene.indx, c("CHR", "ID", "GENE", "P0", "P1", "NPANEL")] <- basic.inf
    
    if (m != 0) {
      # best
      tmp <- res.save[res.save[, "TWAS.P"] == min(res.save[, "TWAS.P"]),, drop=F]
      res.multi.tissue[gene.indx, "BEST.PANEL"] <- as.character(tmp[1, "PANEL"])
      res.multi.tissue[gene.indx, c("BEST.NSNP", "BEST.TWAS.Z", "BEST.TWAS.P")] <- tmp[1, c("NON.ZERO.SNP", "TWAS.Z", "TWAS.P")]
      
      if (m > 1) res.multi.tissue[gene.indx, c("T1m.P", "MultiXcan.P", "InTACT.P")] <- c(T1m.P, MultiXcan.P, InTACT.P)
    }
    
    print(gene.indx)
    run.time <- (end.time - start.time) / 60
    res.multi.tissue[gene.indx, "RUNTIME"] <- run.time
    
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
    basic.inf <- wgtlist1[1, c("CHR", "ID", "GENE", "P0", "P1")]
    res.multi.tissue[gene.indx, c("CHR", "ID", "GENE", "P0", "P1", "NPANEL")] <- basic.inf
    res.multi.tissue[gene.indx, "MOST.SIG.SNP"] <- most.sig.snp
  })
}
cat("Analysis completed.\n")

saveRDS(res.multi.tissue, paste0(opt$out, "_multi.tissue_cHR", opt$chr, ".rds"))
saveRDS(res.single.tissue, paste0(opt$out, "_single.tissue_cHR", opt$chr, ".rds"))
