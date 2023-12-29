library(ggplot2)

source("sampleVarCos.R")

cosdir <- "~/Work/bhk/analysis/cosine/"

if (!file.exists(file.path(cosdir, "datafiles/sampleVarCos.rds"))){
  resdf <- sampleVarCos()
  saveRDS(resdf, file.path(cosdir, "datafiles/sampleVarCos.rds"))
} else {
  resdf <- readRDS(file.path(cosdir, "datafiles/sampleVarCos.rds"))
}



pdf(file=file.path(cosdir, "cosineVarianceApproximation.pdf"), width=8, height=6)
ggplot(resdf, aes(x=theoryVar, y=obsVar, color=as.factor(ndim))) + geom_point() + theme_minimal() + 
  xlab("Predicted Variance") + ylab("Observed variance") + ggtitle(sprintf("Cosine variance for multivariate normal, Pearson = %0.4f", cor(resdf$theoryVar, resdf$obsVar))) + 
  geom_abline(intercept=0, slope=1, col="black", lty=2) + guides(color=guide_legend(title="Dimension")) + 
  theme(text=element_text(size=16))
dev.off()


ggplot(resdf, aes(x=ndim, y=obsVar*ndim, color=as.factor(ndim))) + geom_point() + theme_minimal() + 
  xlab("Dimension") + ylab("Observed Variance/theoretical minimum") + ggtitle("Observed cosine variance from simulation scaled by theoretical minimum") + 
  ylim(c(0, 1.1*max(resdf$obsVar*resdf$ndim))) + xlim(c(0, 1000)) + geom_hline(yintercept=1, col="black", lty=2) + 
  geom_text("Theoretical minimum")




# Case 3

if (!file.exists(file.path(cosdir, "datafiles/sampleVarCosCase3.rds"))){
  
  totaldf <- c()
  for (mysd in seq(0.1, 1, 0.1)){
    tdf <- sampleGenVarCos(ndim=seq(100, 1000, 100), iter=1, distiter=10, 
                           sPar=c(2,2,2,2), nSamples=1000, meanSd=mysd)
    totaldf <- rbind(totaldf, tdf)
  }
  saveRDS(totaldf, file.path(cosdir, "datafiles/sampleVarCosCase3.rds"))
  
} else {
  totaldf <- readRDS(file.path(cosdir, "datafiles/sampleVarCosCase3.rds"))
}

pdf(file=file.path(cosdir, "cosineVarianceApproxCase3.pdf"), width=8, height=6)
ggplot(totaldf, aes(x=theoryVar, y=obsVar, color=as.factor(ndim))) + geom_point() + 
  geom_abline(intercept=0, slope=1, lty=2) + theme_minimal() + xlim(c(0, 0.31)) + ylim(c(0,0.31)) + 
  xlab("Predicted Variance") + ylab("Observed Variance") + 
  ggtitle(sprintf("Cosine variance for multivariate normal, Pearson = %0.4f", cor(totaldf$theoryVar, totaldf$obsVar))) + 
  theme(text=element_text(size=16)) + guides(color=guide_legend(title="Dimension"))
dev.off()

pdf(file=file.path(cosdir, "cosineVarianceApproxCase3Log.pdf"), width=8, height=6)
ggplot(totaldf, aes(x=-log10(theoryVar), y=-log10(obsVar), color=as.factor(ndim))) + geom_point() + 
  geom_abline(intercept=0, slope=1, lty=2) + theme_minimal() + 
  xlab("-Log10 Predicted Variance") + ylab("-Log10 Observed Variance") + 
  ggtitle(sprintf("Cosine variance for multivariate normal, Pearson = %0.4f", cor(-log10(totaldf$theoryVar), -log10(totaldf$obsVar)))) + 
  theme(text=element_text(size=16)) + guides(color=guide_legend(title="Dimension"))
dev.off()


totaldf$varRatio <- totaldf$theoryVar/totaldf$obsVar

pdf(file=file.path(cosdir, "figs", "cosineTheoryVarianceAnomaly.pdf"), width=10, height=8)
ggplot(totaldf, aes(x=obsMean, y=varRatio, color=as.factor(ndim))) + geom_point() + 
  xlab("Empirical Cosine Mean") + ylab("Ratio of theory variance to empirical variance") + 
  theme_minimal() + guides(color=guide_legend(title="Dimension")) + 
  ggtitle("Variation Ratio vs cosine similarity mean")
dev.off()

pdf(file=file.path(cosdir, "figs", "cosineTheoryVarianceAnomalyLog10.pdf"), width=10, height=8)
ggplot(totaldf, aes(x=obsMean, y=log10(varRatio), color=as.factor(ndim))) + geom_point() + 
  xlab("Empirical Cosine Mean") + ylab("Log10 Ratio of theory variance to empirical variance") + 
  theme_minimal() + guides(color=guide_legend(title="Dimension")) + 
  ggtitle("Log10 variance ratio vs cosine similarity mean")
dev.off()

pdf(file=file.path(cosdir, "figs", "cosineTheoryMeanVsEmpirical.pdf"), width=10, height=8)
ggplot(totaldf, aes(x=obsMean, y=theoryMean, color=as.factor(ndim))) + geom_point() + 
  xlab("Empirical Cosine Mean") + ylab("Theoretical Mean") + 
  theme_minimal() + guides(color=guide_legend(title="Dimension")) + 
  ggtitle("Theoretical vs empirical cosine similarity mean")
dev.off()




#### Experiments to identify the source of the anomalous variance
b <- cosineDistGeneral(iter=1, nSamples=500)
b$varRatio <- b$theoryVar/b$obsVar
b$approxVarRatio <- b$theoryVar/b$approxVar

varmax <- 1.1*max(b$theoryVar, b$obsVar, b$approxVar)

pdf(file=file.path(cosdir, "figs", "cosineTheoryVarAnomalyInquiry.pdf"), width=10, height=8)
ggplot(b, aes(x=theoryVar, y=obsVar, color=as.factor(ndim))) + geom_point() + theme_minimal() + 
  geom_abline(slope=1, intercept=0, lty=2) + xlim(c(0, varmax)) + ylim(c(0, varmax)) + 
  xlab("Theoretical Variance") + ylab("Empirical Variance") + ggtitle("Theoretical vs Empirical cosine variance") + 
  guides(color=guide_legend(title="Dimension"))

ggplot(b, aes(x=theoryVar, y=approxVar, color=as.factor(ndim))) + geom_point() + theme_minimal() + 
  geom_abline(slope=1, intercept=0, lty=2) + xlim(c(0, varmax)) + ylim(c(0, varmax)) + 
  xlab("Theoretical Variance") + ylab("Empirical Variance") + ggtitle("Theoretical vs Empirical cosine variance") + 
  guides(color=guide_legend(title="Dimension"))
dev.off()

