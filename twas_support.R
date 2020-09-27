allele.qc <- function(a1, a2, ref1, ref2) {
  ref <- ref1
  flip <- ref
  flip[ref == "A"] <- "T"
  flip[ref == "T"] <- "A"
  flip[ref == "G"] <- "C"
  flip[ref == "C"] <- "G"
  flip1 <- flip
  
  ref <- ref2
  flip <- ref
  flip[ref == "A"] <- "T"
  flip[ref == "T"] <- "A"
  flip[ref == "G"] <- "C"
  flip[ref == "C"] <- "G"
  flip2 <- flip
  
  snp <- list()
  snp[["keep"]] <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  snp[["flip"]] <- (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  
  return(snp)
}


SPrediXcan.sum <- function(wgtlist0, opt, genos) {
  res.save <- as.data.frame(matrix(NA, nrow=8, ncol=nrow(wgtlist0)))
  colnames(res.save) <- wgtlist0[, "PANEL"]
  
  wgt.out <- list()
  
  for (w in 1:nrow(wgtlist0)) {
    tryCatch({
      
      # Load in weights
      wgt.file <- paste0(opt$weights_dir, wgtlist0$WGT[w])
      load(wgt.file)
      
      # Remove NAs (these should not be here)
      wgt.matrix[is.na(wgt.matrix)] <- 0
      wgt.matrix <- as.matrix(wgt.matrix)
      
      snps <- snps[snps[,2] %in% rownames(wgt.matrix),]
      snps <- snps[!duplicated(snps),]
      if ( sum(duplicated(snps[,2])) >0 && sum(duplicated(rownames(wgt.matrix)))==0 ) {
        snps <- snps[!duplicated(snps[,2]),]
      }
      
      # Match up the SNPs and weights
      m <- match(snps[, 2], genos$bim[, 2])
      m.keep <- !is.na(m)
      snps <- snps[m.keep, ]
      wgt.matrix <- wgt.matrix[m.keep,, drop = F]
      
      cur.genos <- scale(genos$bed[, m[m.keep]])
      cur.genos[is.na(cur.genos)] = 0
      cur.bim <- genos$bim[m[m.keep], ]
      # Flip WEIGHTS for mismatching alleles
      qc <- allele.qc(snps[, 5], snps[, 6], cur.bim[, 5], cur.bim[, 6])
      wgt.matrix[qc$flip, ] <- -1 * wgt.matrix[qc$flip, ]
      #rm(snps)
      
      cur.FAIL <- FALSE
      
      # Match up the SNPs and the summary stats
      m <- match(cur.bim[, 2], sumstat$SNP)
      cur.Z <- sumstat$Z[m]
      
      # Identify the best model
      if (ncol(wgt.matrix) == 1) {
        wgt.matrix <- cbind(wgt.matrix, wgt.matrix)
        colnames(wgt.matrix) <- c("enet","dummy")
        mod.best = which(colnames(wgt.matrix) == "enet")
      } else {
        mod.best = "enet"
      }
      
      if (length(mod.best) == 0) {
        cat(w, "WARNING : --force_model", mod.best, "does not exist for", unlist(wgtlist0[w, ]), "\n")
        cur.FAIL <- TRUE
        max.r2 <- -9
      } else max.r2 <- cv.performance[1, mod.best]
      
      if (sum(wgt.matrix[, mod.best] != 0) == 0) {
        cat(w, "WARNING : ", unlist(wgtlist0[w, ]), names(cv.performance)[ mod.best ], "had", length(cur.Z), "overlapping SNPs, but none with non-zero expression weights, try more SNPS or a different model\n")
        cur.FAIL <- TRUE
      }
      
      # Compute LD matrix
      if (length(cur.Z) == 0 ) {
        cat(w, "WARNING : ", unlist(wgtlist0[w, ]), " had no overlapping SNPs\n")
        cur.FAIL <- TRUE
      } else if ( length(cur.Z) == 1 && is.na(cur.Z) ){
        cat(w, "WARNING : ", unlist(wgtlist0[w, ]), " had no overlapping SNPs\n")
        cur.FAIL <- TRUE
      } else if (!cur.FAIL) {
        cur.LD <- t(cur.genos) %*% cur.genos / (nrow(cur.genos) - 1)
        nSNP = nrow(cur.LD)

        cur.miss <- is.na(cur.Z)
        cur.Z <- cur.Z[!cur.miss]
        cur.LD = cur.LD[!cur.miss,, drop=F]
        cur.LD = cur.LD[,!cur.miss, drop=F]
        wgt.matrix <- wgt.matrix[!cur.miss,, drop=F]
        
        # Compute TWAS Z-score
        non_zero_ind = (wgt.matrix[, mod.best] != 0)
        wgt.save <- wgt.matrix[non_zero_ind, mod.best]
        cur.twasz <- wgt.save %*% cur.Z[non_zero_ind]
        cur.twasr2pred <- wgt.save %*% cur.LD[non_zero_ind, non_zero_ind, drop=F] %*% wgt.save
        
        # standardize weights
        diag_element <- as.vector(wgt.matrix[non_zero_ind, mod.best])
        diag_sd <- diag_element / sum(abs(diag_element))
        weight_diag <- diag(diag_sd, nrow = length(diag_sd))
        
        Zstat.w <- weight_diag %*% cur.Z[non_zero_ind]
        corSNP.w <- weight_diag %*% cur.LD[non_zero_ind, non_zero_ind, drop=F] %*% t(weight_diag)
             
        if (cur.twasr2pred > 0) {
          cur.twas <- cur.twasz / sqrt(cur.twasr2pred)
          twas.p <- 2 * (pnorm(abs(cur.twas), lower.tail = F))
          ssu.p <- SumSqU(U = Zstat.w, CovS = corSNP.w,method ="Pan")
          acat.p = ACAT(2 * (pnorm(abs(cur.Z[non_zero_ind]), lower.tail = F)), abs(wgt.save))
          Tk.p = ACAT(c(twas.p, ssu.p, acat.p), rep(1/3, 3))
          
          res.save[, w] <- c(nSNP, sum(non_zero_ind), max.r2, cur.twas, twas.p, ssu.p, acat.p, Tk.p)
          
          if (sum(non_zero_ind) == 1) {
            tmp.name <- rownames(wgt.matrix)
            tmp.name <- tmp.name[non_zero_ind]
            names(wgt.save) <- tmp.name
          }
          wgt.out[[w]] <- wgt.save
          
          
        } else {
          cur.FAIL <- T
          cat(w, "WARNING : ", unlist(wgtlist0[w, ]), " had zero predictive accuracy, try a different model.\n")
        }
      }
    }, error = function(e) { cat(w, "ERROR :", conditionMessage(e), "\n") })
      
  } # end of : for (w in 1:nrow(wgtlist0) )
  
  keep = !duplicated(wgtlist0[, "PANEL"])
  res.save = res.save[,keep,drop=F]
  wgt.out = wgt.out[keep]
  
  wgt.out.all <- unlist(wgt.out)
  wgt.out.all <- unique(names(wgt.out.all))
  
  if (length(wgt.out.all) != 0) {
    wgt <- matrix(0, length(wgt.out.all), sum(keep))
    rownames(wgt) <- wgt.out.all
    colnames(wgt) <- wgtlist0$PANEL[keep]
    
    if (nrow(wgtlist0) != length(wgt.out)) {
      wgt.out[[sum(keep)]] <- "flag"
    }
    
    for (i in 1:sum(keep)) {
      wgt.tmp <- wgt.out[[i]]
      if (!is.null(wgt.tmp)) {
        if (wgt.tmp[1] != "flag") {
          wgt[names(wgt.tmp), i] <- wgt.tmp
        }
      }
    }
  } else {
    wgt <- NULL
  }
  
  return(list(res.save = res.save, wgt = wgt))
}


MultiXcan <- function(wgt.mat, Z.score, used.weight,cutoff = 30) {
  sumstat.tmp <- sumstat[sumstat[, 1] %in% rownames(wgt.mat), ]
  sumstat.tmp <- sumstat.tmp[!is.na(sumstat.tmp[, "Z"]), ] # should check if Z scores are missing....
  
  wgt.mat <- wgt.mat[rownames(wgt.mat) %in% sumstat.tmp[, 1],, drop=F]

  cur.genos <- scale(genos$bed[, rownames(wgt.mat)])
  cur.LD <- t(cur.genos) %*% cur.genos / (nrow(cur.genos) - 1)
  
  
  n.tissue <- dim(wgt.mat)[2]
  n.snp <- dim(wgt.mat)[1]
  
  
  nonzero.indx <- list()
  for (i in 1:n.tissue) {
    wi <- wgt.mat [, i]
    wi <- cbind(wi, 1:n.snp)
    wi <- wi[wi[, 1] != 0, ]
    
    
    if (is.null(dim(wi))) {
      wi <- as.matrix(wi)
      wi <- t(wi)
    }
    
    indx.tmp <- wi[, 2]
    nonzero.indx[[i]] <- indx.tmp - 1 # for C codes
  }
  
  # sourceCpp("twas_perm.cpp")
  
  var <- calcMultiVar(cur.LD, wgt.mat, nonzero.indx)
  t.mat <- calcMultiCov(cur.LD, wgt.mat, nonzero.indx, var)
  colnames(t.mat) <- rownames(t.mat) <- colnames(wgt.mat)
  
  Z.score <- as.numeric(Z.score)
  
  t.mat.used <- t.mat[used.weight, ]
  t.mat.used <- t.mat.used[, used.weight]
  tmp <- svd(t.mat.used)
  eigenvalue <- tmp$d
  eigenvalue[eigenvalue < max(eigenvalue) / cutoff] <- 0
  eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
  
  multi.df <- sum(eigenvalue != 0)
  D <- diag(eigenvalue)
  t.mat2 <- tmp$u %*% D %*% t(tmp$v)
  
  multi.stat <- Z.score %*% t.mat2 %*% Z.score
  multi.p <- pchisq(multi.stat, df = multi.df, lower.tail = F)
  
  
  return(list(tissue.cov = t.mat, multi.stat = multi.stat, multi.p = multi.p, multi.mat = t.mat2, multi.df = multi.df))
}




MultiXcan_brain <- function(wgt.mat, Z.score, used.weight,t.mat,cutoff = 30) {
  sumstat.tmp <- sumstat[sumstat[, 1] %in% rownames(wgt.mat), ]
  sumstat.tmp <- sumstat.tmp[!is.na(sumstat.tmp[, "Z"]), ] # should check if Z scores are missing....
  wgt.mat <- wgt.mat[rownames(wgt.mat) %in% sumstat.tmp[, 1], ]
  
  
  Z.score <- as.numeric(Z.score)
  
  t.mat.used <- t.mat[used.weight, ]
  t.mat.used <- t.mat.used[, used.weight]
  tmp <- svd(t.mat.used)
  eigenvalue <- tmp$d
  eigenvalue[eigenvalue < max(eigenvalue) / cutoff] <- 0
  eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
  
  multi.df <- sum(eigenvalue != 0)
  D <- diag(eigenvalue)
  t.mat2 <- tmp$u %*% D %*% t(tmp$v)
  
  multi.stat <- Z.score %*% t.mat2 %*% Z.score
  multi.p <- pchisq(multi.stat, df = multi.df, lower.tail = F)
  
  
  return(list(tissue.cov = t.mat, multi.stat = multi.stat, multi.p = multi.p, multi.mat = t.mat2, multi.df = multi.df))
}
