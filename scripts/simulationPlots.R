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