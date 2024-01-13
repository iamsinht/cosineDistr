library(cosineDistr)
library(ggplot2)

# This assumes you have run the simulations from cosineDistr::cosineDistGeneral, 
# as well as presupposes naming conventions. It's not particularly general.
makeSimulationPlotWrapper <- function(simpath, outpath){
  
  anorm <- readRDS(file.path(simpath, "cosineDist_normal_100n2000.rds"))
  alapl <- readRDS(file.path(simpath, "cosineDist_laplace_100n2000.rds"))
  amix <- readRDS(file.path(simpath, "cosineDist_mixture_100n1000.rds"))
  amixconst <- readRDS(file.path(simpath, "cosineDist_mixtureconst_100n1000.rds"))
  
  makeSimPlot(anorm, outpath, "normalN=2000_", "Normal")
  makeSimPlot(alapl, outpath, "laplaceN=2000_", "Laplace")
  makeSimPlot(amix, outpath, "mixtureN=1000_", "Gaussian Mixture")
  makeSimPlot(amixconst, outpath, "mixconstN=1000_", "Constrained Gaussian Mixture")
}


makeSimPlot <- function(res, outpath, fname, plotname){
  
  for (ii in seq(9)){ # Fix a weird string formatting of the numerical values
    res[,ii] <- as.numeric(res[,ii])
  }

  pdf(file=file.path(outpath, sprintf("%s_simulationplots.pdf", fname)), width=10, height=8)
  print(ggplot(res, aes(x=obsMean, y=theoryMean, color=as.factor(ndim))) + geom_point() + theme_minimal() + 
    xlab("Observed Mean") + ylab("Theoretical Mean") + geom_abline(slope=1, intercept=0, lty=2) + 
    ggtitle(sprintf("%s observed vs theoretical mean, cor = %0.3f", plotname, cor(res$obsMean, res$theoryMean))) +
    guides(color=guide_legend(title="Dimension")))
  
  print(ggplot(res, aes(x=obsVar, y=theoryVar, color=as.factor(ndim))) + geom_point() + theme_minimal() + 
    xlab("Observed Variance") + ylab("Theoretical Variance") + geom_abline(slope=1, intercept=0, lty=2) + 
    ggtitle(sprintf("%s observed vs theoretical variance, cor = %0.3f", plotname, cor(res$obsVar, res$theoryVar))) +
    guides(color=guide_legend(title="Dimension")))
  
  print(ggplot(res, aes(x=-log10(obsVar), y=-log10(theoryVar), color=as.factor(ndim))) + geom_point() + theme_minimal() + 
    xlab("-log10 Observed Variance") + ylab("-log10 Theoretical Variance") + 
    geom_abline(slope=1, intercept=0, lty=2) + guides(color=guide_legend(title="Dimension")) +
    ggtitle(sprintf("%s observed vs theoretical variance, cor = %0.3f", plotname, cor(res$obsVar, res$theoryVar))))
  
  print(ggplot(res, aes(x=approxVar, y=theoryVar, color=as.factor(ndim))) + geom_point() + theme_minimal() + 
    xlab("Approximation Variance") + ylab("Theoretical Variance") + geom_abline(slope=1, intercept=0, lty=2) + 
    ggtitle(sprintf("%s approximation vs theoretical variance, cor = %0.3f", plotname, cor(res$approxVar, res$theoryVar))) +
    guides(color=guide_legend(title="Dimension")))
  
  print(ggplot(res, aes(x=-log10(approxVar), y=-log10(theoryVar), color=as.factor(ndim))) + geom_point() + theme_minimal() + 
    xlab("-log10 Approximation Variance") + ylab("-log10 Theoretical Variance") + 
    geom_abline(slope=1, intercept=0, lty=2) + guides(color=guide_legend(title="Dimension")) +
    ggtitle(sprintf("%s approximation vs theoretical variance, cor = %0.3f", plotname, cor(res$approxVar, res$theoryVar))))
  
  print(ggplot(res, aes(x=obsMean, y=log10(theoryVar/obsVar), color=as.factor(ndim))) + geom_point() + theme_minimal() + 
    xlab("Observed Mean") + ylab("log10 TheoryVar/ObservedVar") + geom_hline(yintercept=0, lty=2) + 
    ggtitle(sprintf("%s theory/obs variance ratio vs Mean, cor = %0.3f", plotname, cor(res$obsMean, res$theoryVar/res$obsVar))) +
    guides(color=guide_legend(title="Dimension")))
  
  dev.off()
}


combineSimulations <- function(simpath, outpath=".", figname=""){
  # This function takes any simulations that have been run with cosineDistGeneral and
  # outputs a number of summary plots and tables.
  # Grab only the elements with up to 2500 dimensions
  myf <- list.files(path=simpath, pattern="cosineDist.*n2500.*rds")
 
  simdata <- lapply(myf, FUN=function(x) readRDS(file.path(simpath, x))) 
  
  
  resdf <- purrr::reduce(simdata, rbind)
  
  resdf$varRatio <- resdf$theoryVar/resdf$obsVar
  resdf$logVarRatio <- log10(resdf$theoryVar/resdf$obsVar)
  
  meanRatioLims <- c(0.9 * min(resdf$meanLen/resdf$theoryLen), 1.01*max(resdf$meanLen/resdf$theoryLen))
  lvrLims <- c(1.05 * min(resdf$logVarRatio), 1.05*max(resdf$logVarRatio))
  
  # The distributions with mean 0, i.e. centered:
  ix <- which(resdf$meanSd == 0)
  
  pdf(file.path(outpath, sprintf("cosDistrSimSummaryPlots%s.pdf", figname)), width=10, height=10)
  
  # Means
  print(ggplot(resdf, aes(x=obsMean, y=theoryMean, shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("Observed mean") + ylab("Theoretical mean") + 
          ggtitle(sprintf("Theoretical vs empirical mean of cosine similarity, simulated distributions, N=%d", dim(resdf)[1])) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  
  # Variance
  print(ggplot(resdf, aes(x=obsVar, y=theoryVar, shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("Observed Variance") + ylab("Theoretical Variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance of cosine similarity, simulated distributions, N=%d", dim(resdf)[1])) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  
  # Log10 variance
  print(ggplot(resdf, aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("-log10 Theoretical variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance of cosine similarity, simulated distributions, N=%d", dim(resdf)[1])) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  
  # Variance ratio vs mean similarity
  print(ggplot(resdf, aes(x=obsMean, y=logVarRatio, shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("Observed mean") + ylab("log10 TheoryVar/ObsVar") + 
          ggtitle(sprintf("Theoretical vs empirical variance of cosine similarity, simulated distributions, N=%d", dim(resdf)[1])) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  
  
  # Plots per distribution - Variance
  # varPlots <- lapply(unique(resdf$distribution), FUN=function(x) ggplot(resdf[resdf$distribution == x, ], aes(x=obsVar, y=theoryVar, shape=distribution, color=as.factor(ndim))) + 
  #   geom_point() + theme_minimal() + xlab("Observed Variance") + ylab("Theoretical Variance") + 
  #   ggtitle(sprintf("Theoretical vs empirical variance, %s, N=%d", x, dim(resdf)[1])) + 
  #   guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  # 
  # print(grid.arrange(varPlots[[1]], varPlots[[2]], varPlots[[3]], varPlots[[4]], nrow=2))
  
  print(ggplot(resdf, aes(x=obsVar, y=theoryVar, shape=distribution, color=as.factor(ndim))) + 
      geom_point() + theme_minimal() + xlab("Observed Variance") + ylab("Theoretical Variance") +
      ggtitle(sprintf("Theoretical vs empirical variance")) +
      guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2) + 
        facet_wrap(vars(distribution)))

  # Plots per distribution: Mean vs theory mean 
  print(ggplot(resdf, aes(x=obsMean, y=theoryMean, shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("Observed mean") + ylab("Theoretical mean") + 
          ggtitle(sprintf("Theoretical vs empirical mean of cosine similarity, simulated distributions, N=%d", dim(resdf)[1])) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2) + 
          facet_wrap(vars(distribution)))
    
  # Plots per distribution - log10 variance
  # LVarPlots <- lapply(unique(resdf$distribution), FUN=function(x) ggplot(resdf[resdf$distribution == x, ], 
  #  aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
  #  geom_point() + theme_minimal() + xlab("-log10 Observed Variance") + ylab("-log10 Theoretical Variance") + 
  #  ggtitle(sprintf("Theoretical vs empirical variance, %s, N=%d", x, dim(resdf)[1])) + 
  #  guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  # print(grid.arrange(LVarPlots[[1]], LVarPlots[[2]], LVarPlots[[3]], LVarPlots[[4]], nrow=2))
  
  print(ggplot(resdf, aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
    geom_point() + theme_minimal() + xlab("-log10 Observed Variance") + ylab("-log10 Theoretical Variance") + 
    ggtitle(sprintf("Theoretical vs empirical variance (-log10)")) + 
    guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2) + 
    facet_wrap(vars(distribution)))
  
  # Plots per distribution - Log Var Ratio vs observed mean
  # mVsLVarPlots <- lapply(unique(resdf$distribution), FUN=function(x) ggplot(resdf[resdf$distribution == x, ],
  #  aes(x=obsMean, y=logVarRatio, shape=distribution, color=as.factor(ndim))) +
  #  geom_point() + theme_minimal() + xlab("Observed Mean") + ylab("log10 TheoryVar/ObsVar") +
  #  ggtitle(sprintf("Variance ratio vs observed mean, %s; N=%d", x, dim(resdf)[1])) +
  #  guides(color=guide_legend(title="Dimension")) + geom_hline(yintercept=0, lty=2) +
  #  ylim(c(lvrLims[1], lvrLims[2])))
  # print(grid.arrange(mVsLVarPlots[[1]], mVsLVarPlots[[2]], mVsLVarPlots[[3]], mVsLVarPlots[[4]], nrow=2))
  
  print(ggplot(resdf, aes(x=obsMean, y=logVarRatio, shape=distribution, color=as.factor(ndim))) +
    geom_point() + theme_minimal() + xlab("Observed Mean") + ylab("log10 TheoryVar/ObsVar") +
    ggtitle(sprintf("Variance ratio vs observed mean")) +
    guides(color=guide_legend(title="Dimension")) + geom_hline(yintercept=0, lty=2) +
    ylim(c(lvrLims[1], lvrLims[2])) + 
    facet_wrap(vars(distribution)))
           
  
  # Plots per distribution - log Var Ratio vs meanLen/theoryLen - does the VarRatio break down when the 
  # mean assumption breaks down?
  # mratioVsLVarPlots <- lapply(unique(resdf$distribution), FUN=function(x) ggplot(resdf[resdf$distribution == x, ], 
  #  aes(x=meanLen/theoryLen, y=logVarRatio, shape=distribution, color=as.factor(ndim))) + 
  #  geom_point() + theme_minimal() + xlab("Observed Mean") + ylab("log10 TheoryVar/ObsVar") + 
  #  ggtitle(sprintf("Variance ratio vs Expected Length, %s; N=%d", x, dim(resdf)[1])) + 
  #  xlim(c(meanRatioLims[1], meanRatioLims[2])) +
  #  guides(color=guide_legend(title="Dimension")) + geom_hline(yintercept=0, lty=2) + 
  #  ylim(c(lvrLims[1], lvrLims[2])))
  # print(grid.arrange(mratioVsLVarPlots[[1]], mratioVsLVarPlots[[2]], mratioVsLVarPlots[[3]], mratioVsLVarPlots[[4]], nrow=2))
  print(ggplot(resdf, aes(x=meanLen/theoryLen, y=logVarRatio, shape=distribution, color=as.factor(ndim))) + 
    geom_point() + theme_minimal() + xlab("Mean Length/Theory Length") + ylab("log10 TheoryVar/ObsVar") + 
    ggtitle(sprintf("Log10 Variance Ratio (Theory/Obs) vs Observed Length/Theory")) + 
    xlim(c(meanRatioLims[1], meanRatioLims[2])) +
    guides(color=guide_legend(title="Dimension")) + geom_hline(yintercept=0, lty=2) + 
    ylim(c(lvrLims[1], lvrLims[2])) + 
      facet_wrap(vars(distribution)))
  
  # Plots per distribution - logVarRatio vs log10 obsVar - over what regimes is the approximation wrong?
  # LVarObsVarPlots <- lapply(unique(resdf$distribution), FUN=function(x) ggplot(resdf[resdf$distribution == x, ],
  #   aes(x=-log10(obsVar), y=logVarRatio, shape=distribution, color=as.factor(ndim))) +
  #   geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("log10 TheoryVar/ObsVar") +
  #   ggtitle(sprintf("Variance ratio vs observed variance, %s; N=%d", x, dim(resdf)[1])) +
  #   guides(color=guide_legend(title="Dimension")) + geom_hline(yintercept=0, lty=2) +
  #   ylim(c(lvrLims[1], lvrLims[2])))
  # print(grid.arrange(LVarObsVarPlots[[1]], LVarObsVarPlots[[2]], LVarObsVarPlots[[3]], LVarObsVarPlots[[4]], nrow=2))
          
  print(ggplot(resdf, aes(x=-log10(obsVar), y=logVarRatio, shape=distribution, color=as.factor(ndim))) +
    geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("log10 TheoryVar/ObsVar") +
    ggtitle(sprintf("Log10 Variance ratio vs observed variance")) +
    guides(color=guide_legend(title="Dimension")) + geom_hline(yintercept=0, lty=2) +
    ylim(c(lvrLims[1], lvrLims[2])) + 
    facet_wrap(vars(distribution)))
  
  dev.off()
  
  
  # For zero mean guys:
  pdf(file.path(outpath, sprintf("cosDistrSimSummaryPlots_ZeroMean%s.pdf", figname)), width=10, height=10)
  #print(ggplot(resdf[ix,], aes(x=obsMean, y=theoryMean, shape=distribution, color=as.factor(ndim))) + 
  #        geom_point() + theme_minimal() + xlab("Observed mean") + ylab("Theoretical mean") + 
  #        ggtitle(sprintf("Theoretical vs empirical mean of cosine similarity, simulated distributions, N=%d", dim(resdf[ix,])[1])) + 
  #        guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  
  # Variance
  print(ggplot(resdf[ix,], aes(x=obsVar, y=theoryVar, shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("Observed Variance") + ylab("Theoretical Variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance of cosine similarity, simulated distributions, N=%d", dim(resdf[ix,])[1])) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  
  # Log10 variance observed vs theory
  print(ggplot(resdf[ix,], aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("-log10 Theoretical variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance of cosine similarity, simulated distributions, N=%d", dim(resdf[ix,])[1])) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  
  # Log10 variance observed vs theory per distribution
  print(ggplot(resdf[ix,], aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("-log10 Theoretical variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance of cosine similarity, simulated distributions, N=%d", dim(resdf[ix,])[1])) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2) + 
          facet_wrap(vars(distribution)))
  
  
  # Variance ratio vs mean similarity
  #print(ggplot(resdf[ix,], aes(x=obsMean, y=logVarRatio, shape=distribution, color=as.factor(ndim))) + 
  #        geom_point() + theme_minimal() + xlab("Observed mean") + ylab("log10 TheoryVar/ObsVar") + 
  #        ggtitle(sprintf("Theoretical vs empirical variance of cosine similarity, simulated distributions, N=%d", dim(resdf[ix,])[1])) + 
  #        guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  dev.off()
  
  # Show eigenvalue plots
  eigdata <- data.frame(x=numeric(), y=numeric(), iter=numeric())
  myx <- seq(0, 100, 0.01)
  
  for (ii in seq(25)){
    a <- stats::rgamma(1, shape=2, rate=2)
    b <- stats::rgamma(1, shape=2, rate=2)
    
    eigdata <- rbind(eigdata, data.frame(x=myx + 1e-6, y=dgamma(myx, shape=a, rate=b), iter=ii))
  }
  
  pdf(file.path(outpath, sprintf("cosDistr_simulationEigenvalues.pdf")), width=8, height=6)
  print(ggplot(eigdata, aes(x=x, y=y, color=as.factor(iter))) + geom_line(lwd=1.5) + 
          coord_cartesian(ylim=c(0, 0.5), xlim=c(0, 20)) + xlab("Eigenvalues") + 
          ylab("Probability distribution") + ggtitle("Example eigenvalue distributions for simulations") + 
        guides(color=FALSE) + theme_minimal())
  dev.off()
  
  
  #### Manuscript figures: ####
  # Log10 variance observed vs theory for centered data
  pdf(file.path(outpath, sprintf("cosDistrSimulation_Log10Var_ZeroMean.pdf")), width=6, height=6)
  print(ggplot(resdf[ix,], aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("-log10 Theoretical variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance, centered distributions")) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  dev.off()
  
  # Log10 Theoretical variance vs observed variance
  pdf(file.path(outpath, sprintf("cosDistrSimulation_Log10Var_General.pdf")), width=6, height=6)
  print(ggplot(resdf, aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("-log10 Theoretical variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance, general distributions")) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  dev.off()
  
  # Observed Mean vs Theoretical Mean
  pdf(file.path(outpath, sprintf("cosDistrSimulation_Means_General.pdf")), width=6, height=6)
  print(ggplot(resdf, aes(x=obsMean, y=theoryMean, shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("Observed mean") + ylab("Theoretical mean") + 
          ggtitle(sprintf("Theoretical vs empirical mean")) +
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  dev.off()
  
  # Supp - Per distribution: Obs Mean vs Log10 Variance ratio
  pdf(file.path(outpath, sprintf("cosDistrSimulations_MeanVsLog10Var_GeneralPerDistr.pdf")), width=6, height=6)
  print(ggplot(resdf, aes(x=obsMean, y=logVarRatio, color=as.factor(ndim))) + 
       geom_point() + theme_minimal() + xlab("Observed Mean") + ylab("log10 TheoryVar/ObsVar") +
         ggtitle("Log10 Variance ratio vs Observed Mean") + 
       guides(color=guide_legend(title="Dimension")) + geom_hline(yintercept=0, lty=2) +
       ylim(c(lvrLims[1], lvrLims[2])) + 
    facet_wrap(vars(distribution)))
  dev.off()
  
  # Supp - Per distribution: Log10 theoretical vs observed variance
  pdf(file.path(outpath, sprintf("cosDistrSimulation_Log10Var_GeneralPerDistr.pdf")), width=6, height=6)
  print(ggplot(resdf, aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("-log10 Observed Variance") + ylab("-log10 Theoretical Variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance (-log10)")) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2) + 
          facet_wrap(vars(distribution)))
  dev.off()

  
  #### Defense figures: ####
  pdf(file.path(outpath, sprintf("DefensecosDistrSimulation_Log10Var_General.pdf")), width=6.5, height=5)
  print(ggplot(resdf, aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("-log10 Theoretical variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance, general distributions")) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  dev.off()
  
  pdf(file.path(outpath, sprintf("DefensecosDistrSimulation_Means_General.pdf")), width=6.5, height=5)
  print(ggplot(resdf, aes(x=obsMean, y=theoryMean, shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("Observed mean") + ylab("Theoretical mean") + 
          ggtitle(sprintf("Theoretical vs empirical mean")) +
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  dev.off()
  
  pdf(file.path(outpath, sprintf("DefensecosDistrSimulation_Log10Var_ZeroMean.pdf")), width=6.5, height=5)
  print(ggplot(resdf[ix,], aes(x=-log10(obsVar), y=-log10(theoryVar), shape=distribution, color=as.factor(ndim))) + 
          geom_point() + theme_minimal() + xlab("-log10 Observed variance") + ylab("-log10 Theoretical variance") + 
          ggtitle(sprintf("Theoretical vs empirical variance, centered distributions")) + 
          guides(color=guide_legend(title="Dimension")) + geom_abline(slope=1, intercept=0, lty=2))
  dev.off()
  
  
  
    
  saveRDS(resdf, file=file.path(outpath, "cosineDistributionSimulations.rds"))
  
  distSum0 <- resdf[ix,] %>% group_by(distribution) %>% 
    summarise(count=n(), modelCor=cor(obsVar, theoryVar), nCor=cor(obsVar, 1/ndim), 
              modelLogCor=cor(-log10(obsVar), -log10(theoryVar)), nLogCor=cor(-log10(obsVar), -log10(1/ndim)))
  saveRDS(distSum0, file=file.path(outpath, "cosineSimulationsMean0Table.rds"))
  
  distSum <- resdf %>% group_by(distribution) %>% 
    summarise(count=n(), modelCor=cor(obsVar, theoryVar), nCor=cor(obsVar, 1/ndim),
    modelLogCor=cor(-log10(obsVar), -log10(theoryVar)), nLogCor=cor(-log10(obsVar), -log10(1/ndim)))
  saveRDS(distSum, file=file.path(outpath, "cosineSimulationSummaryTable.rds"))
  
  pdf(file.path(outpath, sprintf("cosineSimulationCorrelationGeneralBar.pdf")), width=7, height=5)
  print(ggplot(reshape2::melt(res$distSum, id.vars=c(1,2)), aes(x=distribution, y=value, fill=variable)) + 
    geom_bar(stat="identity", position=position_dodge()) + theme_minimal() + 
    ggtitle("Correlation of simulation variance with theoretical variance") + 
    ylab("Pearson's correlation") + guides(fill=guide_legend(title="Model")))
  dev.off()
  
  pdf(file.path(outpath, sprintf("cosineSimulationCorrelationZeroMeanBar.pdf")), width=7, height=5)
  print(ggplot(reshape2::melt(res$distSum0, id.vars=c(1,2)), aes(x=distribution, y=value, fill=variable)) + 
          geom_bar(stat="identity", position=position_dodge()) + theme_minimal() + 
          ggtitle("Correlation of simulation variance with theoretical variance, Centered data") + 
          ylab("Pearson's correlation") + guides(fill=guide_legend(title="Model")))
  dev.off()
  
  return(list(resdf=resdf, distSum0=distSum0, distSum=distSum))
  
}