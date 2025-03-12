library(perturbKit)
library(cmapR)
library(cosineDistr)
library(ggrepel)

# The purpose of this script is to illustrate a practical application of the cosine 
# distribution work. 

L1KcosExample <- function(l1kdatapath, l1kmetapath, outdir=".", mycellid="HA1E", 
                          method="smd", pclim=50){
  
  method <- match.arg(method, c("smd", "minVar0"))
  
  l1kmeta <- perturbKit::read_l1k_meta(l1kmetapath, version=2020)
  mysigs <- l1kmeta$siginfo[l1kmeta$siginfo$cell_id == mycellid & l1kmeta$siginfo$pert_type == "trt_cp",]
  
  cpcount <- table(mysigs$pert_iname)
  topcps <- names(sort(cpcount, decreasing=TRUE)[1:50])
  
  ds <- parse_gctx(perturbKit::get_level5_ds(l1kdatapath), rid=l1kmeta$landmarks$pr_gene_id, 
                   cid=mysigs$sig_id)
  
  # Obtain the eigenvalues and means for the whole dataset along those eigenvectors
  # Oldx is the original dataset in the basis of the eigenvectors of the covariance matrix
  remap <- optimizeVar(ds@mat)
  
  xpr <- prcomp(t(ds@mat), center=TRUE, scale=FALSE)
  
  oldx <- t(remap$xpr$x) + as.numeric(remap$myMeans)
  
  # Don't use this
  if (method == "minVar0"){
    
    pcrange <- seq(pclim)
  
    ax <- sample(dim(ds@mat)[2], 1e4)  
    bgsim0 <- cosineDistr::cosine(oldx[pcrange, ax], oldx[pcrange,ax])
    bgsim1 <- cosineDistr::cosine(newx[pcrange, ax], newx[pcrange,ax])
    
    u0 <- mean(bgsim0[upper.tri(bgsim0)])
    u1 <- mean(bgsim1[upper.tri(bgsim1)])
    
    sd0 <- sd(bgsim0[upper.tri(bgsim0)])
    sd1 <- sd(bgsim1[upper.tri(bgsim1)])
    
    q0 <- quantile(bgsim0[upper.tri(bgsim0)], 0.99)
    q1 <- quantile(bgsim1[upper.tri(bgsim1)], 0.99)
    
    L1Kex <- data.frame(cp=character(), tstat0=numeric(), tstat1=numeric(), p0=numeric(), p1=numeric())
    
    # Consider only the top 100 PCs, as the bottom ones are likely noise
    for (icp in topcps){
      print(icp)
      
      cpix <- which(mysigs$pert_iname == icp)
      cp0 <- cosineDistr::cosine(oldx[pcrange,cpix], oldx[pcrange,cpix])
      cp1 <- cosineDistr::cosine(newx[pcrange,cpix], newx[pcrange,cpix])
      
      t0 <- (mean(cp0[upper.tri(cp0)]) - u0)/sd0
      t1 <- (mean(cp1[upper.tri(cp1)]) - u1)/sd1
      
      qcp0 <- mean(cp0[upper.tri(cp0)] > q0)
      qcp1 <- mean(cp1[upper.tri(cp1)] > q1)
      
      L1Kex <- rbind(L1Kex, data.frame(cp=icp, tstat0=t0, tstat1=t1, p0=qcp0, p1=qcp1))
    }
    
    params <- list(pclim=pclim, mycellid=mycellid)
    
    return(list(L1Kex=L1Kex, params=params))
  }
  
  # Alternate approach, explicitly maximizing standardized mean difference.
  # I originally described this as the t-statistic, but this is incorrect because the t-statistic
  # varies with n, whereas SMD does not. see Zhang 2010 "Strictly standardized mean difference" for more.
  if (method == "smd"){
    
    bx <- sample(dim(ds@mat)[2], 2e3)
    
    L1Ktx <- data.frame(cp=character(), tstat0=numeric(), tstat1=numeric(), p0=numeric(), p1=numeric())
   
    for (icp in topcps){
      print(icp)
      cpix <- which(mysigs$pert_iname == icp)
      
      # Calculate the new eigenvalues
      mysigma <- getMaxTstat(rowMeans(oldx[,cpix]), rowMeans(oldx), xpr$sdev, iter=5)
      
      # Rescale to get the new eigenvalues
      # Is this correct?
      optimx <- mysigma$sigma/xpr$sdev * oldx
      
      # Calculate the statistics
      cp0 <- cosineDistr::cosine(oldx[, cpix], oldx[, cpix])
      cpOptim <- cosineDistr::cosine(optimx[, cpix], optimx[, cpix])
      
      bg0 <- cosineDistr::cosine(oldx[, bx], oldx[, bx])
      bgOptim <- cosineDistr::cosine(optimx[, bx], optimx[, bx])
      
      tstat0 <- (mean(cp0[upper.tri(cp0)]) - mean(bg0[upper.tri(bg0)]))/sd(bg0[upper.tri(bg0)])
      tstat1 <- (mean(cpOptim[upper.tri(cpOptim)]) - mean(bgOptim[upper.tri(bgOptim)]))/sd(bgOptim[upper.tri(bgOptim)])
      
      p0 <- mean(cp0[upper.tri(cp0)] > quantile(bg0[upper.tri(bg0)], 0.99))
      p1 <- mean(cpOptim[upper.tri(cpOptim)] > quantile(bgOptim[upper.tri(bgOptim)], 0.99))
      
      L1Ktx <- rbind(L1Ktx, data.frame(cp=icp, tstat0=tstat0, tstat1=tstat1, p0=p0, p1=p1))
  
      if (TRUE){
        pdf(file.path(outdir, sprintf("%s_optim_%s_%d.pdf", icp, mycellid, length(cpix))), 
            width=6, height=6)
        plot(density(bgOptim[upper.tri(bgOptim)], bw=0.01), col="blue", lwd=2, 
             main=sprintf("%s, Transformed space", icp), xlab="Cosine similarity (bw = 0.01)", 
             xlim=c(-0.5, 0.75))
        lines(density(cpOptim[upper.tri(cpOptim)], bw=0.01), lwd=2, col="red")
        legend(x="topright", legend=c("Compound", "Background"), col=c("red", "blue"), lwd=c(4,4))
        dev.off() 
        
        pdf(file.path(outdir, sprintf("%s_base_%s_%d.pdf", icp, mycellid, length(cpix))), 
            width=6, height=6)
        plot(density(bg0[upper.tri(bg0)], bw=0.01), lwd=2, col="navy",  
             main=sprintf("%s, Unmodified space", icp), 
             xlab="Cosine similarity (bw = 0.01)", xlim=c(-0.5, 0.75))
        lines(density(cp0[upper.tri(cp0)], bw=0.01), col="orange", lwd=2)
        legend(x="topright", legend=c("Compound", "Background"), col=c("darkorange", "navy"), 
               lwd=c(4,4))
        dev.off()
        
      }
    }

    # Save twice for safety
    saveRDS(list(L1Ktx=L1Ktx, params=c(cellid=mycellid, iter=5, backgroundelts=2000)), 
            file = file.path(outdir, sprintf("L1KTstatSimulation_%s.rds", mycellid)))
        
    L1Ktx$moa <- l1kmeta$pertinfo[match(L1Ktx$cp, l1kmeta$pertinfo$cmap_name),"moa"]
    
    saveRDS(list(L1Ktx=L1Ktx, params=c(cellid=mycellid, iter=5, backgroundelts=2000)), 
            file = file.path(outdir, sprintf("L1KTstatSimulation_%s.rds", mycellid)))
    return(L1Ktx)
  }
  
  
}


# This function takes as input a dataset and returns a transformed dataset 
# per Smith et al "On the distribution of cosine similarity"
optimizeVar <- function(mymat){
  
  xpr <- prcomp(t(mymat), center=TRUE, scale=FALSE)
  lambda <- xpr$sdev^2
  myMeans <- xpr$center %*% xpr$rotation
  
  eta <- myMeans/xpr$sdev
  
  # Without loss of generality (scale invariance of cosine similarity), set C = 1
  
  minlambda <- (1 + eta^2)/(1+2*eta^2)
  
  scaleFactor <- minlambda/lambda
  
  # Calculate the rescaled basis and the old data in the new basis. Need to adjust for the 
  # centering used in the PCA:
  newx <- as.vector(scaleFactor) * t(xpr$x)
  oldx <- t(xpr$x) + as.numeric(myMeans)
  
  return(xpr=xpr, oldx=oldx, newx=newx, myMeans=myMeans)
}


# A different approach:
# Maximize t-stat given two distribution means and the eigenvalues sigma
getMaxTstat <- function(m1, m0, sigma, iter=5){
  
  # Make dimensionless
  e1 <- m1/sigma
  e0 <- m0/sigma
  
  # Compute base tstat:
  tstatBase <- tstatExpr(e1, e0, sigma, 0, 0)
  objvals <- numeric()
  
  for (jj in seq(iter)){
    print(jj)
    for (ii in seq_along(sigma)){
      x <- optimize(tstatExpr, interval=c(0, max(abs(sigma))), tol=1e-3, 
                    maximum=TRUE, e1=e1, e0=e0, sigma=sigma, index=ii)
      
      sigma[ii] <- x$maximum
      objvals <- c(objvals, x$objective)
    }
  }
  
  return(list(sigma=sigma, objvals=objvals))
}


tstatExpr <- function(e1, e0, sigma, index, val){
  if (index != 0){
    sigma[index] <- val
  }
  return((sum(sigma^2 * e1^2)/sum(sigma^2*(1+e1^2)) - sum(sigma^2* e0^2)/sum(sigma^2* (1+e0^2))) * 
           (sum(sigma^4 * (1 + 2*e0^2))/sum(sigma^2 * (1 + e0^2))^2)^(-1/2))
}


makeL1KexpPlots <- function(L1Ktx, outdir, mycellid){
  
  maxt <- 1.05*max(max(L1Ktx$tstat0), max(L1Ktx$tstat1))
  
  pdf(file.path(outdir, sprintf("l1kcpExperiment_smdstat_%s_%d.pdf", mycellid, dim(L1Ktx)[1])), 
      width=6, height=6)
  ggplot(L1Ktx, aes(x=tstat0, y=tstat1, label=cp)) + geom_point(color="blue") + theme_minimal() + 
    xlab("Unmodified data, Compound SMD") + ylab("Transformed data, SMD") + 
    ggtitle(sprintf("L1000 optimized embedding, %s N=%d", mycellid, dim(L1Ktx)[1])) + 
    geom_abline(slope=1, intercept=0, lty=2) + xlim(c(0, maxt)) + ylim(c(0, maxt)) + 
    geom_text_repel(size=2.5) + guides(color="none")
  dev.off()
  
  pdf(file.path(outdir, sprintf("l1kcpExperiment_plt01_%s_%d.pdf", mycellid, dim(L1Ktx)[1])), 
      width=6, height=6)
  ggplot(L1Ktx, aes(x=p0, y=p1, label=cp)) + geom_point(color="forestgreen") + theme_minimal() + 
    xlab("Unmodified data, rank < 1%") + ylab("Transformed data, rank < 1%") + 
    ggtitle(sprintf("L1000 optimized embedding, %s N=%d", mycellid, dim(L1Ktx)[1])) + 
    geom_abline(slope=1, intercept=0, lty=2) + xlim(c(0, 1)) + ylim(c(0, 1)) + 
    geom_text_repel(size=2.5) + guides(color="none")
  dev.off()
  
  pdf(file.path(outdir, sprintf("l1kcpExperiment_plt01_%s_%d_7x5.pdf", mycellid, dim(L1Ktx)[1])), 
      width=7, height=5)
  ggplot(L1Ktx, aes(x=p0, y=p1, label=cp)) + geom_point(color="forestgreen") + theme_minimal() + 
    xlab("Original data power, alpha = 0.01") + ylab("Transformed data power, alpha = 0.01") + 
    ggtitle(sprintf("L1000 optimized embedding, %s N=%d", mycellid, dim(L1Ktx)[1])) + 
    geom_abline(slope=1, intercept=0, lty=2) + xlim(c(0, 1)) + ylim(c(0, 1)) + 
    geom_text_repel(size=2.5) + guides(color="none")
  dev.off()
  
}
