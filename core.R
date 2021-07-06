################################################################################################################################ #
## Load dependencies ----     

load.deps <- function() {
  #library(DeconRNASeq)
  library(dtangle)
  library(e1071)
  library(parallel)
  library(preprocessCore)
  library(Metrics)
  library(DTWBI)
}

################################################################################################################################ #
## Functions ----     

## Function to calculate GoF
write.gof <- function(measuredExp, estimatedComp, signatureUsed, returnPred = FALSE, log = TRUE) {
  # set to common order
    commonGenes <- rownames(measuredExp)[which(rownames(measuredExp) %in% rownames(signatureUsed))]
    measuredExp <- measuredExp[commonGenes,]; signatureUsed <- signatureUsed[commonGenes,]
    
  # predict expression (predExp) from the estimatedComp * signatureUsed
    predExp <- as.data.frame(matrix(nrow = length(commonGenes), ncol = ncol(measuredExp)))
    rownames(predExp) <- commonGenes
    
    for(j in 1:ncol(predExp)) {
      # storage
      a <- list()
      
      # the contribution of each cell-type to predicted expression
      for(k in colnames(signatureUsed)) { a[[k]] <- estimatedComp[j,k] * signatureUsed[,k] }
      
      # sum expression from all cell-types to a single predicted value
      predExp[,j] <- rowSums(do.call("cbind", a))
    }
  
  # statistics
      # (optional) log transform
      if(log) {
        measuredExp <- log2(measuredExp + 0.5)
        predExp <- log2(predExp + 0.5)
      }
    
    
    # setup
    stats <- as.data.frame(matrix(ncol = 2, nrow = ncol(measuredExp)))
    colnames(stats) <- c("rho", "r")
    rownames(stats) <- colnames(measuredExp)
    
    # spearman correlation between predicted and estimated expression
    stats$rho <- diag(cor(measuredExp, predExp, method = "s")) 
    
    # pearson correlation (comments per spearman)
    stats$r <- diag(cor(measuredExp, predExp, method = "p")) 
    
    # # mean absolute error, root mean squared error, reconstruction accuracy, and cosine similarity
    # for(j in 1:ncol(measuredExp)) { 
    #   a <- measuredExp[,j]
    #   b <- predExp[,j]
    #   
    #   stats$mae[j] <- mae(a, b) 
    #   stats$rmse[j] <- rmse(a, b) 
    #   stats$recon[j] <- 1 - (((eNorm(a - b)) ^ 2) / ((eNorm(a)) ^ 2)) # note: this (I believe) is what PsychEncode calls the "reconstruction accuracy"
    #   stats$cosine[j] <- (sum(a * b) / (eNorm(a) * eNorm(b))) # note: this is the cosine similarity of two vectors per Wikipedia: dot product / product of euclidean norms
    # }
    # 
  # return
    if(returnPred) {
      res <- list()
      res$predExp <- predExp
      res$stats <- stats  
    } else {
      res <- stats
    }
  
    return(res)
}

write.gof.v2 <- function(measuredExp, estimatedComp, signatureUsed, returnPred = FALSE) {
  # set to common row order
  commonGenes <- rownames(measuredExp)[which(rownames(measuredExp) %in% rownames(signatureUsed))]
  measuredExp <- measuredExp[commonGenes,]; signatureUsed <- signatureUsed[commonGenes,]  
  
  # quantile normalise    
  qn <- data.frame(signatureUsed, measuredExp)
  qn <- as.data.frame(normalize.quantiles(as.matrix(qn), copy = FALSE))
  signatureUsed <- qn[,1:ncol(signatureUsed)]
  measuredExp <- qn[,-c(1:ncol(signatureUsed))]  
  
  # predict expression (predExp) from the estimatedComp * signatureUsed
  predExp <- as.data.frame(matrix(nrow = length(commonGenes), ncol = ncol(measuredExp)))
  rownames(predExp) <- commonGenes    
  
  for(j in 1:ncol(predExp)) {
    # storage
    a <- list()      
    
    # the contribution of each cell-type to predicted expression
    for(k in colnames(signatureUsed)) { a[[k]] <- estimatedComp[j,k] * signatureUsed[,k] }     
    
    # sum expression from all cell-types to a single predicted value
    predExp[,j] <- rowSums(do.call("cbind", a))
  }
    
  ## Calculate statistics
  stats <- as.data.frame(matrix(ncol = 5, nrow = ncol(measuredExp)))
  colnames(stats) <- c("rho", "r", "mae", "rmse", "nmae")
  rownames(stats) <- colnames(measuredExp)    
  for(j in 1:ncol(measuredExp)) { 
    a <- measuredExp[,j]
    b <- predExp[,j]      
    stats$r[j] <- cor(log2(a+0.5), log2(b+0.5), method = "p")
    stats$rho[j] <- cor(a, b, method = "s")      
    stats$mae[j] <- mae(a, b) 
    stats$rmse[j] <- rmse(a, b)       
    stats$nmae[j] <- compute.nmae(a, b)     
  }  
  # return
  if(returnPred) {
    res <- list()
    res$predExp <- predExp
    res$stats <- stats  
  } else {
    res <- stats
  }    
  
  return(res)
}

## Code to run CIBERSORT
run.CIB <- function(mixture, signature, interruptCallback, progressSet, progressStart, progressScale, currStep, totalSteps) {

  # run CIB
  res <- CIBERSORT(sig_matrix=signature, mixture_file=mixture, interruptCallback=interruptCallback, progressSet=progressSet, progressStart=progressStart, progressScale=progressScale, currStep=currStep, totalSteps=totalSteps)
    
  # reformat
  res <- as.data.frame(res)
  res <- res[,1:(ncol(res) - 3)] # Note1: This removes some of CIBERSORT's output that is extraneous (for our purposes)
  return(res)
}


## Code to run DRS
run.DRS <- function(mixture, signature) {
  res <- as.data.frame(DeconRNASeq(mixture, signature, use.scale = TRUE)$out.all)
  rownames(res) <- colnames(mixture)
  return(res)
}

## Source code for CIBERSORT
    # CIBERSORT R script v1.04 (last updated 10-24-2016)
    # Note: Signature matrix construction is not currently available; use java version for full functionality.
    # Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
    # Requirements:
    #       R v3.0 or later. (dependencies below might not work properly with earlier versions)
    #       install.packages('e1071')
    #       install.pacakges('parallel')
    #       install.packages('preprocessCore')
    #       if preprocessCore is not available in the repositories you have selected, run the following:
    #           source("http://bioconductor.org/biocLite.R")
    #           biocLite("preprocessCore")
    # Windows users using the R GUI may need to Run as Administrator to install or update packages.
    # This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
    # single-threaded in Windows.
    #
    # Usage:
    #       Navigate to directory containing R script
    #
    #   In R:
    #       source('CIBERSORT.R')
    #       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)
    #
    #       Options:
    #       i)   perm = No. permutations; set to >=100 to calculate p-values (default = 0)
    #       ii)  QN = Quantile normalization of input mixture (default = TRUE)
    #       iii) absolute = Run CIBERSORT in absolute mode (default = FALSE)
    #               - note that cell subsets will be scaled by their absolute levels and will not be
    #                 represented as fractions (to derive the default output, normalize absolute
    #                 levels such that they sum to 1 for each mixture sample)
    #               - the sum of all cell subsets in each mixture sample will be added to the ouput
    #                 ('Absolute score'). If LM22 is used, this score will capture total immune content.
    #       iv)  abs_method = if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
    #               - sig.score = for each mixture sample, define S as the median expression
    #                 level of all genes in the signature matrix divided by the median expression
    #                 level of all genes in the mixture. Multiple cell subset fractions by S.
    #               - no.sumto1 = remove sum to 1 constraint
    #
    # Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
    # Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
    # License: http://cibersort.stanford.edu/CIBERSORT_License.txt
    
    
    #dependencies
    #library(e1071)
    #library(parallel)
    #library(preprocessCore)
    
    #Core algorithm
    CoreAlg <- function(X, y, absolute, abs_method){
      
      #try different values of nu
      svn_itor <- 3
      
      res <- function(i){
        if(i==1){nus <- 0.25}
        if(i==2){nus <- 0.5}
        if(i==3){nus <- 0.75}
        model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
        model
      }
      
      if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
        out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
      
      nusvm <- rep(0,svn_itor)
      corrv <- rep(0,svn_itor)
      
      #do cibersort
      t <- 1
      while(t <= svn_itor) {
        weights = t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights<0)]<-0
        w<-weights/sum(weights)
        u <- sweep(X,MARGIN=2,w,'*')
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- cor(k, y)
        t <- t + 1
      }
      
      #pick best model
      rmses <- nusvm
      mn <- which.min(rmses)
      model <- out[[mn]]
      
      #get and normalize coefficients
      q <- t(model$coefs) %*% model$SV
      q[which(q<0)]<-0
      if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
      if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
      
      mix_rmse <- rmses[mn]
      mix_r <- corrv[mn]
      
      newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
      
    }
    
    #do permutations
    doPerm <- function(perm, X, Y, absolute, abs_method){
      itor <- 1
      Ylist <- as.list(data.matrix(Y))
      dist <- matrix()
      
      while(itor <= perm){
        #random mixture
        yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
        
        #standardize mixture
        yr <- (yr - mean(yr)) / sd(yr)
        
        #run CIBERSORT core algorithm
        result <- CoreAlg(X, yr, absolute, abs_method)
        
        mix_r <- result$mix_r
        
        #store correlation
        if(itor == 1) {dist <- mix_r}
        else {dist <- rbind(dist, mix_r)}
        
        itor <- itor + 1
      }
      newList <- list("dist" = dist)
    }
    
    #main function
    ## Note: GJS modified this function to no longer need to read in matrices from a file, but instead using variables loaded into the workspace
    ## Note: KAW modified this function for shiny progress bar and to fix duplicates
    CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE, absolute=FALSE, abs_method='sig.score', interruptCallback=function(){}, progressSet=function(){}, progressStart=0.0, progressScale=1.0, currStep=1, totalSteps=1){
      
      if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score') stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
      
      #read in data
      # X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
      # Y <- read.table(mixture_file, header=T, sep="\t",check.names=F)
      X <- sig_matrix
      Y <- mixture_file
      
      #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
      dups <- dim(Y)[1] - length(unique(rownames(Y)))
      if(dups > 0) {
        warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
        rownames(Y) <- make.names(rownames(Y), unique=TRUE)
      }
      
      X <- data.matrix(X)
      Y <- data.matrix(Y)
      
      #order
      X <- X[order(rownames(X)),]
      Y <- Y[order(rownames(Y)),]
      
      P <- perm #number of permutations
      
      #anti-log if max < 50 in mixture file
      if(max(Y) < 50) {Y <- 2^Y}
      
      #quantile normalization of mixture file
      if(QN == TRUE){
        tmpc <- colnames(Y)
        tmpr <- rownames(Y)
        Y <- normalize.quantiles(Y)
        colnames(Y) <- tmpc
        rownames(Y) <- tmpr
      }
      
      #store original mixtures
      Yorig <- Y
      Ymedian <- max(median(Yorig),1)
      
      #intersect genes
      Xgns <- row.names(X)
      Ygns <- row.names(Y)
      YintX <- Ygns %in% Xgns
      Y <- Y[YintX,]
      XintY <- Xgns %in% row.names(Y)
      X <- X[XintY,]
      
      #standardize sig matrix
      X <- (X - mean(X)) / sd(as.vector(X))
      
      #empirical null distribution of correlation coefficients
      if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}
      
      header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
      if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
      
      output <- matrix()
      itor <- 1
      mixtures <- dim(Y)[2]
      pval <- 9999
      
      #iterate through mixtures
      time.eta <- "This step can take a while"
      while(itor <= mixtures){
        start.time <- Sys.time()
        
        interruptCallback()

        progressSet(value=progressStart + progressScale * (itor-1) / mixtures, message=sprintf("Step %d/%d: Running CIBERSORT", currStep, totalSteps), detail=sprintf("%s, please DO NOT refresh or close the page", time.eta))
        
        y <- Y[,itor]
        
        #standardize mixture
        y <- (y - mean(y)) / sd(y)
        
        #run SVR core algorithm
        result <- CoreAlg(X, y, absolute, abs_method)
        
        #get results
        w <- result$w
        mix_r <- result$mix_r
        mix_rmse <- result$mix_rmse
        
        if(absolute && abs_method == 'sig.score') {
          w <- w * median(Y[,itor]) / Ymedian
        }
        
        #calculate p-value
        if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
        
        #print output
        out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
        if(absolute) out <- c(out, sum(w))
        if(itor == 1) {output <- out}
        else {output <- rbind(output, out)}
        
        itor <- itor + 1

        end.time <- Sys.time()
        time.eta <- sprintf("ETA: %s", format((end.time - start.time) * (mixtures - itor + 1), format="%H:%M:%S"))
      }
      
      #save results
      #write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
      
      #return matrix object containing all results
      obj <- rbind(header,output)
      obj <- obj[,-1]
      obj <- obj[-1,]
      obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
      rownames(obj) <- colnames(Y)
      if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
      else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
      obj
    }

    run.DTA <- function(mixture, signature, alg = "diff", q = 0.01) { # alg %in% c("diff", "ratio", "p.value", "regression")
      # bind dataframes
      common <- rownames(mixture)[which(rownames(mixture) %in% rownames(signature))]
      dat <- cbind(log2(mixture[common,] + 0.5), log2(signature[common,] + 0.5))
      dat <- as.data.frame(t(dat))  # define where signature samples reside in dat
      ps <- list()
      for (j in colnames(signature)) { ps[[j]] <- grep(paste0("^", j, "$"), rownames(dat)) }  # find markers
      markers <- find_markers(dat, marker_method = alg, data_type = "rna-seq", pure_samples = ps)  # deconvolve
      quant <- lapply(markers$V, function(x) { quantile(x, 1-q) } )
      for (j in 1:length(quant)) {
        if (quant[j] == Inf) { quant[j] <- max(grep("Inf", markers$V[[j]]) + 1) } # a hack for when there are more "infs" than q markers
      }
      n <- sapply(1:length(markers$V), function(i) { max(which(markers$V[[i]] > quant[[i]])) } )
      res <- as.data.frame(dtangle(Y = dat, 
                                  pure_samples = ps, 
                                  markers = markers$L, 
                                  n_markers = n, 
                                  marker_method = alg, 
                                  data_type = "rna-seq")$estimates)  # return(res)
      return(res[1:(nrow(res) - ncol(signature)),])
    }

    
